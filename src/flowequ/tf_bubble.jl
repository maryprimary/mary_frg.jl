"""
温度流的Bubble
"""



"""
获取所有的Bubble(pp, fs, ex)
qpp = k1+k2; qfs = k3-k2; qex = k1-k3
"""
function all_bubble_tf_mt(Γ4::Gamma4{T, P}, lval) where {T, P}
    #
    lamb = Γ4.λ_0 * exp(-lval)
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
        q_pp = kadd(Γ4.model, k1, k2)
        bubbres = pi_αβ_minus_tf(
            Γ4.model, tris_pat[i_n], brlu_area, lamb,
            q_pp, alpha, beta
        )
        bubbval_pp[alpha, beta, b1, b2, i_n, i1, i2] = bubbres
    end
    bubb_qpp = Bubble{:tf, :minus}(lval, bubbval_pp)
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
        q_fs = kadd(Γ4.model, k3, -k2)
        bubbres = pi_αβ_plus_tf(
            Γ4.model, tris_pat[i_n], brlu_area,
            lamb, q_fs, alpha, beta
        )
        bubbval_fs[alpha, beta, b2, b3, i_n, i2, i3] = bubbres
    end
    bubb_qfs = Bubble{:tf, :plus}(lval, bubbval_fs)
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
        q_ex = kadd(Γ4.model, k1, -k3)
        bubbres = pi_αβ_plus_tf(
            Γ4.model, tris_pat[i_n], brlu_area, lamb,
            q_ex, alpha, beta
        )
        bubbval_ex[alpha, beta, b1, b3, i_n, i1, i3] = bubbres
    end
    bubb_qex = Bubble{:tf, :plus}(lval, bubbval_ex)
    return bubb_qpp, bubb_qfs, bubb_qex
end



"""
温度流的Bubble
"""
function pi_αβ_plus_tf(
    model::P, ltris::Vector{T},
    area::Float64, lamb::Float64,
    qval::Point2D,
    dαidx::Int64, dβidx::Int64) where {
        T <: Basics.AbstractTriangle, P <: Fermi.Abstract2DModel
        }
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
        kprim = kval + nega_q
        #kprim = kadd(kval, nega_q)
        #epsilon_k#dispα(kval.x, kval.y)
        eps_k = dispersion(model, dαidx, kval.x, kval.y)
        #epsilon_{k-q}#dispβ(kprim.x, kprim.y)
        eps_kp = dispersion(model, dβidx, kprim.x, kprim.y)
        #这个小区域的贡献
        if abs(eps_k - eps_kp) < 1.e-10
            #如果特别小，可以利用
            # lim (eps_k -> eps_kp) Pi^{+} =
            # 1/T (e^{eps/T} (-eps/T*e^{eps/T} +
            #eps/T + e^{eps/T} + 1)) / (e^{eps/T} + 1)^3
            bval = eps_kp / lamb
            #如果本身就很大，分母会比较大导致接近0
            if bval > 25
                #@info "数值不稳定"
                num = 0.
                den = 100.
            else
                expb = exp(bval)
                num = expb * (-bval * expb + bval + expb + 1)
                den = (1+expb)^3
            end
            d_val = num / den / lamb
        else#如果本身是正常的数值
            #分子左侧的数值
            if (eps_k / lamb) > 25
                #@info "数值不稳定"
                num_left = 0.
            else
                #exp^{epsilon_k / T}
                exp_k_t = exp(eps_k / lamb)
                num_left = eps_k / lamb * exp_k_t / ((1 + exp_k_t)^2)
            end
            #分子右侧的数值
            if (eps_kp / lamb) > 25
                #@info "数值不稳定"
                num_righ = 0.
            else
                #e^{epsilon_{k-q} / T}
                exp_kp_t = exp(eps_kp / lamb)
                num_righ = eps_kp / lamb * exp_kp_t / ((1 + exp_kp_t)^2)
            end
            d_val = (num_left - num_righ) / (eps_k - eps_kp)
        end# end if 小区域的贡献计算完
        result += d_val * tri.area
    end#对ltris的循环
    result = result / area
    return result
end


function pi_αβ_minus_tf(
    model::P, ltris::Vector{T},
    area::Float64, lamb::Float64,
    qval::Point2D,
    dαidx::Int64, dβidx::Int64) where {
        T <: Basics.AbstractTriangle, P <: Fermi.Abstract2DModel
        }
    result = 0.
    for tri in ltris
        #这个小三角形的k值
        kval = tri.center
        nega_k = -kval
        #-k+q
        kprim = nega_k + qval
        #kprim = kadd(nega_k, qval)
        #epsilon_k, #dispα(kval.x, kval.y)
        eps_k = dispersion(model, dαidx, kval.x, kval.y)
        #-epsilon_{-k+q}#-dispβ(kprim.x, kprim.y)
        neps_kp = -dispersion(model, dβidx, kprim.x, kprim.y)
        #这个时候，因为epsilon_{-k+q}前面已经有了负号，分母上还是负号
        #计算这个小区域的贡献
        if abs(eps_k - neps_kp) < 1.e-10
            #如果两个数值比较接近, Pi^{-}和Pi^{+}的公式完全一样，就是第二个能量要加个负号
            # lim (eps_k -> -eps_kp) Pi^{-} =
            # 1/T (e^{eps/T} (-eps/T*e^{eps/T} + eps/T + e^{eps/T} + 1)) 
            #/ (e^{eps/T} + 1)^3
            bval = eps_k / lamb
            if bval > 25
                #@info "数值不稳定"
                num = 0.
                den = 100.
            else
                expb = exp(bval)
                num = expb * (-bval * expb + bval + expb + 1)
                den = (1+expb)^3
            end
            d_val = num / den / lamb
        else
            if (eps_k / lamb) > 25
                #@info "数值不稳定"
                num_left = 0.
            else
                #e^{epsilon_k / T}
                exp_k_t = exp(eps_k / lamb)
                num_left = eps_k / lamb * exp_k_t / ((1 + exp_k_t)^2)
            end
            if (neps_kp / lamb) > 25
                #@info "数值不稳定"
                num_righ = 0.
            else
                #e^{-epsilon_{-k+q} / T}
                exp_nkp_t = exp(neps_kp / lamb)
                num_righ = neps_kp / lamb * exp_nkp_t / 
                ((1 + exp_nkp_t)^2)
            end
            d_val = (num_left - num_righ) / (eps_k - neps_kp)
        end # end if 小区域的贡献计算完
        result += d_val * tri.area
    end#对ltris的循环
    result = result / area
    return result
end

