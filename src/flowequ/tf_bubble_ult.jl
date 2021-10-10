"""
温度流的Bubble，ultra low temperature
"""

"""
在极低温的情况下，将只有k和k′都在费米面上面的时候才有贡献
"""



"""
获取所有的Bubble(pp, fs, ex)
qpp = k1+k2; qfs = k3-k2; qex = k1-k3
"""
function all_bubble_tf_ult_mt(Γ4::Gamma4{T, P}, cert::Float64) where {T, P}
    #
    brlu_area = area(Γ4.model.brillouin)
    #因为现在所有的能带都必须有同一个lpats，所以只要每一个patch中有哪些tri
    tris_pat = Γ4.ltris_pat
    #
    bubbval_pp = Array{Float64, 7}(
        undef,
        #alpha, beta, b1, b2
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k1(b1), k2(b2)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    Threads.@threads for idxs in CartesianIndices(bubbval_pp)
        alpha, beta, b1, b2, i_n, i1, i2 = Tuple(idxs)
        k1, k2 = Γ4.patches[b1][i1], Γ4.patches[b2][i2]
        q_pp = Γ4.model.kadd(k1, k2)
        bubbres = pi_αβ_minus_tf_ult(
            tris_pat[i_n], brlu_area, cert,
            q_pp, Γ4.model.dispersion[alpha], Γ4.model.dispersion[beta],
            Γ4.model.kadd
        )
        bubbval_pp[alpha, beta, b1, b2, i_n, i1, i2] = bubbres
    end
    bubb_qpp = Bubble{:tf, :minus}(Inf, bubbval_pp)
    #
    bubbval_fs = Array{Float64, 7}(
        undef,
        #alpha, beta, b2, b3
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k2(b2), k3(b3)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    Threads.@threads for idxs in CartesianIndices(bubbval_fs)
        alpha, beta, b2, b3, i_n, i2, i3 = Tuple(idxs)
        k2, k3 = Γ4.patches[b2][i2], Γ4.patches[b3][i3]
        q_fs = Γ4.model.kadd(k3, -k2)
        bubbres = pi_αβ_plus_tf_ult(
            tris_pat[i_n], brlu_area, cert,
            q_fs, Γ4.model.dispersion[alpha],
            Γ4.model.dispersion[beta], Γ4.model.kadd
        )
        bubbval_fs[alpha, beta, b2, b3, i_n, i2, i3] = bubbres
    end
    bubb_qfs = Bubble{:tf, :plus}(Inf, bubbval_fs)
    #
    bubbval_ex = Array{Float64, 7}(
        undef,
        #alpha, beta, b1, b3
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k1(b3), k3(b3)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    Threads.@threads for idxs in CartesianIndices(bubbval_ex)
        alpha, beta, b1, b3, i_n, i1, i3 = Tuple(idxs)
        k1, k3 = Γ4.patches[b1][i1], Γ4.patches[b3][i3]
        q_ex = Γ4.model.kadd(k1, -k3)
        bubbres = pi_αβ_plus_tf_ult(
            tris_pat[i_n], brlu_area, cert,
            q_ex, Γ4.model.dispersion[alpha],
            Γ4.model.dispersion[beta], Γ4.model.kadd
        )
        bubbval_ex[alpha, beta, b1, b3, i_n, i1, i3] = bubbres
    end
    bubb_qex = Bubble{:tf, :plus}(Inf, bubbval_ex)
    return bubb_qpp, bubb_qfs, bubb_qex
end




"""
温度流的Bubble, 这里必须要保证ltris在费米面上面
"""
function pi_αβ_plus_tf_ult(
    ltris::Vector{T},
    area::Float64, cert::Float64,
    qval::Point2D,
    dispα::Function, dispβ::Function,
    kadd::Function) where T <: Basics.AbstractTriangle
    """温度流的+
    这里的lamb就是T，ltris中的所有三角都应该要在同一个patch中,
    tarea是每个小三角形的面积(已经不用了)，
    dispa是和k相关的那个能带，dispb是k-q相关的
    """
    nega_q = -qval
    result = 0.
    for tri in ltris
        #这个小三角形的k值
        kval = tri.center
        #k-q
        kprim = kadd(kval, nega_q)
        #epsilon_k
        eps_k = dispα(kval.x, kval.y)
        if !isapprox(eps_k, 0., atol=cert)
            continue
        end
        #epsilon_{k-q}
        eps_kp = dispβ(kprim.x, kprim.y)
        if !isapprox(eps_kp, 0., atol=cert)
            continue
        end
        #这个小区域的贡献
        result += tri.edges[1].length
    end#对ltris的循环
    result = result / area
    return result
end



"""
温度流的Bubble, 这里必须要保证ltris在费米面上面
"""
function pi_αβ_minus_tf_ult(
    ltris::Vector{T},
    area::Float64, cert::Float64,
    qval::Point2D,
    dispα::Function, dispβ::Function,
    kadd::Function) where T <: Basics.AbstractTriangle
    result = 0.
    for tri in ltris
        #这个小三角形的k值
        kval = tri.center
        nega_k = -kval
        #-k+q
        kprim = kadd(nega_k, qval)
        #epsilon_k
        eps_k = dispα(kval.x, kval.y)
        if !isapprox(eps_k, 0., atol=cert)
            continue
        end
        #-epsilon_{-k+q}
        neps_kp = -dispβ(kprim.x, kprim.y)
        if !isapprox(neps_kp, 0., atol=cert)
            continue
        end
        #这个时候，因为epsilon_{-k+q}前面已经有了负号，分母上还是负号
        #计算这个小区域的贡献
        result += tri.edges[1].length
    end#对ltris的循环
    result = result / area
    return result
end

