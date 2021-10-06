"""
Γ^{4}，有效势
"""

module FlowEquation
    

using ..Basics
using ..Triangulated
using ..Refine
using ..Fermi
using ..Fermi.Surface
using ..Fermi.Patch


export Gamma4
export ECGamma4
export pi_αβ_plus_ec, pi_αβ_minus_ec
export all_bubble_ec_mt

export dl_ec_mt

export TFGamma4
export pi_αβ_plus_tf, pi_αβ_minus_tf
export all_bubble_tf_mt

export dl_tf_mt

export fliter_away_surface
export refine_to_surface
export refine_list_triangle

export TFGamma4_refine_ltris_mt
export TFGamma4_reweight_ltris_mt
export TFGamma4_addition_ltris_mt


"""
最后一个动量可以由动量守恒确定
"""
struct Gamma4{T <: Fermi.Abstract2DModel, P <: Basics.AbstractTriangle}
    model :: T
    λ_0 :: Float64
    V :: Array{Float64, 7}
    k4tab :: Array{Int64, 7}
    #patch允许不同能带取不同的位置，但是为了方便必须有一样的patch数目
    #对相同的布里渊区，patch切分的方法相同，所以相同的patchnum和ltris意味着相同的lpats
    patchnum :: Int64
    patches :: Vector{Vector{Basics.Point2D}}
    ltris :: Vector{P}
    lpats :: Vector{Int64}
    ladjs :: Union{Nothing, TYPE_LADJS{P}}
    #Vector{Tuple{Union{Missing, P}, Union{Missing, P}, Union{Missing, P}}}
    ltris_pat :: Union{Nothing, Vector{Vector{P}}}
end



"""
能量截断的Γ4
"""
function ECGamma4(
    model::Union{QuadrateSystem, TriangularSystem},
    λ_0::Float64, patchnum::Int64, splitnum::Int64
    )
    mpats::Vector{Vector{Basics.Point2D}} = []
    for midx in 1:1:model.bandnum
        push!(mpats, patches_under_vonhove(
            model.brillouin, model.dispersion[midx], patchnum
        ))
    end
    #
    V = zeros(
        model.bandnum, model.bandnum, model.bandnum, model.bandnum,
        patchnum, patchnum, patchnum
    )
    #
    k4tab = Array{Int64, 7}(undef,
    model.bandnum, model.bandnum, model.bandnum, model.bandnum,
    patchnum, patchnum, patchnum
    )
    #
    find_algo = isa(model, QuadrateSystem) ? find_patch_index_squa :
    find_patch_index_hexa
    for idxs in CartesianIndices(k4tab)
        m1, m2, m3, m4, k1, k2, k3 = Tuple(idxs)
        k1p = mpats[m1][k1]
        k2p = mpats[m2][k2]
        k3p = mpats[m3][k3]
        k4p = model.kadd(k1p, k2p)
        k4p = model.kadd(k4p, -k3p)
        k4tab[m1, m2, m3, m4, k1, k2, k3] =
        find_algo(k4p, model.brillouin, patchnum)
    end
    #
    split_algo = isa(model, QuadrateSystem) ? split_square :
    split_hexagon
    ltris, ladjs = split_algo(model.brillouin, splitnum)
    #
    lpats = group_ltris_into_patches_mt(ltris, model.brillouin, patchnum)
    return Gamma4(
        model,
        λ_0,
        V,
        k4tab,
        patchnum,
        mpats,
        ltris, lpats, ladjs, nothing
    )
end



function TFGamma4(
    model::Union{QuadrateSystem, TriangularSystem},
    λ_0::Float64, patchnum::Int64, splitnum::Int64
    )
    mpats::Vector{Vector{Basics.Point2D}} = []
    for midx in 1:1:model.bandnum
        push!(mpats, patches_under_vonhove(
            model.brillouin, model.dispersion[midx], patchnum
        ))
    end
    #
    V = zeros(
        model.bandnum, model.bandnum, model.bandnum, model.bandnum,
        patchnum, patchnum, patchnum
    )
    #
    k4tab = Array{Int64, 7}(undef,
    model.bandnum, model.bandnum, model.bandnum, model.bandnum,
    patchnum, patchnum, patchnum
    )
    #
    find_algo = isa(model, QuadrateSystem) ? find_patch_index_squa :
    find_patch_index_hexa
    for idxs in CartesianIndices(k4tab)
        m1, m2, m3, m4, k1, k2, k3 = Tuple(idxs)
        k1p = mpats[m1][k1]
        k2p = mpats[m2][k2]
        k3p = mpats[m3][k3]
        k4p = model.kadd(k1p, k2p)
        k4p = model.kadd(k4p, -k3p)
        k4tab[m1, m2, m3, m4, k1, k2, k3] =
        find_algo(k4p, model.brillouin, patchnum)
    end
    #
    split_algo = isa(model, QuadrateSystem) ? split_square :
    split_hexagon
    ltris, ladjs = split_algo(model.brillouin, splitnum)
    lpats = group_ltris_into_patches_mt(ltris, model.brillouin, patchnum)
    #
    #因为现在所有的能带都必须有同一个lpats，所以只要每一个patch中有哪些tri
    ltris_pat = Vector{typeof(ltris)}(undef, patchnum)
    for idx in 1:1:patchnum
        ltris_pat[idx] = []
    end
    #
    for (tri, pat) in zip(ltris, lpats)
        push!(ltris_pat[pat], tri)
    end
    #
    return Gamma4(
        model,
        λ_0,
        V,
        k4tab,
        patchnum,
        mpats,
        ltris, lpats, nothing, ltris_pat
    )
end



include("bubble.jl")

include("ec_bubble.jl")



"""
能量截断的导数
"""
function dl_ec_mt(Γ4::Gamma4, bubb_pp::T1, bubb_fs::T2,
    bubb_nfs::T3, bubb_ex::T4, bubb_nex::T5) where {
        T1<:Bubble, T2<:Bubble,
        T3<:Bubble, T4<:Bubble, T5<:Bubble}
    #
    sys = Γ4.model
    dl_val = zeros(size(Γ4.V))
    #所有的自由度
    Threads.@threads for idxs in CartesianIndices(dl_val)
        b1, b2, b3, b4, i1, i2, i3 = Tuple(idxs)
        i4 = Γ4.k4tab[b1, b2, b3, b4, i1, i2, i3]
        #需要进行的积分
        value = 0.
        place_holder = Array{Int8, 3}(undef, sys.bandnum, sys.bandnum, Γ4.patchnum)
        for intidxs in CartesianIndices(place_holder)
            α, β, i_n = Tuple(intidxs)
            pi_min_αβ_n_qpp = bubb_pp.V[α, β, b1, b2, i_n, i1, i2]
            pi_plu_αβ_n_qfs = bubb_fs.V[α, β, b2, b3, i_n, i2, i3]
            pi_plu_αβ_n_nqfs = bubb_nfs.V[α, β, b2, b3, i_n, i2, i3]
            pi_plu_αβ_n_qex = bubb_ex.V[α, β, b1, b3, i_n, i1, i3]
            pi_plu_αβ_n_nqex = bubb_nex.V[α, β, b1, b3, i_n, i1, i3]
            #
            value += Γ4.V[b2, b1, α, β, i2, i1, i_n] *
            Γ4.V[b3, b4, α, β, i3, i4, i_n] * pi_min_αβ_n_qpp
            #
            value += Γ4.V[b1, b2, α, β, i1, i2, i_n] *
            Γ4.V[b4, b3, α, β, i4, i3, i_n] * pi_min_αβ_n_qpp
            #
            value += 2*Γ4.V[α, b4, b1, β, i_n, i4, i1] *
            Γ4.V[α, b2, b3, β, i_n, i2, i3] * pi_plu_αβ_n_qfs
            #
            value += 2*Γ4.V[α, b1, b4, β, i_n, i1, i4] *
            Γ4.V[α, b3, b2, β, i_n, i3, i2] * pi_plu_αβ_n_nqfs
            #
            value -= Γ4.V[b4, α, b1, β, i4, i_n, i1] *
            Γ4.V[α, b2, b3, β, i_n, i2, i3] * pi_plu_αβ_n_qfs
            #
            value -= Γ4.V[b1, α, b4, β, i1, i_n, i4] *
            Γ4.V[α, b3, b2, β, i_n, i3, i2] * pi_plu_αβ_n_nqfs
            #
            value -= Γ4.V[α, b4, b1, β, i_n, i4, i1] *
            Γ4.V[b2, α, b3, β, i2, i_n, i3] * pi_plu_αβ_n_qfs
            #
            value -= Γ4.V[α, b1, b4, β, i_n, i1, i4] *
            Γ4.V[b3, α, b2, β, i3, i_n, i2] * pi_plu_αβ_n_nqfs
            #
            value -= Γ4.V[b3, α, b1, β, i3, i_n, i1] *
            Γ4.V[b2, α, b4, β, i2, i_n, i4] * pi_plu_αβ_n_qex
            #
            value -= Γ4.V[b1, α, b3, β, i1, i_n, i3] *
            Γ4.V[b4, α, b2, β, i4, i_n, i2] * pi_plu_αβ_n_nqex
            #
        end# 
        ##总的负号
        value = -value
        ##太大的时候就线暂停
        if abs(value) > 1e32
            value = 0
        end
        #计算完成
        dl_val[b1, b2, b3, b4, i1, i2, i3] = value
    end
    return dl_val
end


include("tf_bubble.jl")


"""
温度流的导数
"""
function dl_tf_mt(Γ4::Gamma4, bubb_pp::T1, bubb_fs::T2, bubb_ex::T3
    ) where {T1<:Bubble, T2<:Bubble, T3<:Bubble}
    #
    sys = Γ4.model
    dl_val = zeros(size(Γ4.V))
    #所有的自由度
    Threads.@threads for idxs in CartesianIndices(dl_val)
        b1, b2, b3, b4, i1, i2, i3 = Tuple(idxs)
        i4 = Γ4.k4tab[b1, b2, b3, b4, i1, i2, i3]
        #需要进行的积分
        value = 0.
        place_holder = Array{Int8, 3}(undef, sys.bandnum, sys.bandnum, Γ4.patchnum)
        for intidxs in CartesianIndices(place_holder)
            α, β, i_n = Tuple(intidxs)
            pi_min_αβ_n_qpp = bubb_pp.V[α, β, b1, b2, i_n, i1, i2]
            pi_plu_αβ_n_qfs = bubb_fs.V[α, β, b2, b3, i_n, i2, i3]
            pi_plu_αβ_n_qex = bubb_ex.V[α, β, b1, b3, i_n, i1, i3]
            #
            value += Γ4.V[b2, b1, α, β, i2, i1, i_n] *
                Γ4.V[b3, b4, α, β, i3, i4, i_n] * pi_min_αβ_n_qpp
            value += 2*Γ4.V[α, b4, b1, β, i_n, i4, i1] *
                Γ4.V[α, b2, b3, β, i_n, i2, i3] * pi_plu_αβ_n_qfs
            value -= Γ4.V[b4, α, b1, β, i4, i_n, i1] *
                Γ4.V[α, b2, b3, β, i_n, i2, i3] * pi_plu_αβ_n_qfs
            value -= Γ4.V[α, b4, b1, β, i_n, i4, i1] *
                Γ4.V[b2, α, b3, β, i2, i_n, i3] * pi_plu_αβ_n_qfs
            value -= Γ4.V[b3, α, b1, β, i3, i_n, i1] *
                Γ4.V[b2, α, b4, β, i2, i_n, i4] * pi_plu_αβ_n_qex
        end#结束一个微分数值的计算
        #
        value = -value
        if abs(value) > 1.e32
            value = 0.
        end
        #计算完成
        dl_val[b1, b2, b3, b4, i1, i2, i3] = value
    end# 结束对所有能带和动量的循环
    return dl_val
end



"""
取掉离费米面太远的
"""
function fliter_away_surface(Γ4::Gamma4{T, P}, lval) where {T, P}
    lamb = Γ4.λ_0 * exp(-lval)
    resltris::Vector{P} = []
    reslpats::Vector{Int64} = []
    cert = 1#length(Γ4.ltris) > 2000000 ? 5 : 15
    for (tri, pat) in zip(Γ4.ltris, Γ4.lpats)
        engs = [Γ4.model.dispersion[bidx](
            tri.center.x, tri.center.y) / lamb for bidx in 1:1:Γ4.model.bandnum
        ]
        #如果还有可能有贡献
        if minimum(abs.(engs)) < cert #&& minimum(abs.(engs2)) > cert # && pat == 1# && minimum(abs.(engs)) > 0.1
            push!(resltris, tri)
            push!(reslpats, pat)
        end
    end
    ltris_pat = Vector{typeof(resltris)}(undef, Γ4.patchnum)
    for idx in 1:1:Γ4.patchnum
        ltris_pat[idx] = []
    end
    #
    for (tri, pat) in zip(resltris, reslpats)
        push!(ltris_pat[pat], tri)
    end
    #
    return Gamma4(
        Γ4.model,
        Γ4.λ_0,
        Γ4.V,
        Γ4.k4tab,
        Γ4.patchnum,
        Γ4.patches,
        resltris, reslpats, nothing, ltris_pat
    )
end


"""
将Γ4中的ltris进行重新的分割
"""
function TFGamma4_refine_ltris_mt(Γ4::Gamma4{T, P}, lval; maxnum=2500000) where {T, P}
    #
    if !isnothing(Γ4.ladjs)
        @warn "会丢失ladjs的信息"
    end
    if length(Γ4.ltris) >= maxnum
        @warn "已经超过最大数量"
        return fliter_away_surface(Γ4, lval)
    end
    #
    lamb = Γ4.λ_0 * exp(-lval)
    max2suf = Vector{Float64}(undef, Γ4.model.bandnum)
    for bidx in 1:1:length(max2suf)
        disp = Γ4.model.dispersion[bidx]
        suf = const_energy_triangle(Γ4.ltris, 0., disp)
        engs = zeros(length(suf))
        for (idx, tri) in enumerate(suf)
            center = tri.center
            engs[idx] = abs(disp(center.x, center.y)) / lamb
        end
        max2suf[bidx] = maximum(engs)
    end
    #如果已经足够贴近费米面，直接返回原来的
    if maximum(max2suf) < 0.0001
        return fliter_away_surface(Γ4, lval)
    end
    @info string(length(Γ4.ltris))*" -> "*string(4*length(Γ4.ltris))
    #
    #将所有的nltris进行重新的切分，切分一定是切成4个
    spltris = Vector{P}(undef, 4*length(Γ4.ltris))
    for (idx, tri) in enumerate(Γ4.ltris)
        spltris[4*idx-3:4*idx] = split_triangle(tri)
    end
    #去掉离得太远的
    resltris::Vector{P} = []
    cert = length(spltris) > 1000000 ? 5 : 15
    for tri in spltris
        engs = [Γ4.model.dispersion[bidx](
            tri.center.x, tri.center.y) / lamb for bidx in 1:1:Γ4.model.bandnum
        ]
        #如果还有可能有贡献
        if minimum(abs.(engs)) < cert
            push!(resltris, tri)
        end
    end
    #重新计算一下所属的patch，可能会发生变化
    reslpats = group_ltris_into_patches_mt(resltris, Γ4.model.brillouin, Γ4.patchnum)
    #重新计算一下每个patch包含的三角形
    ltris_pat = Vector{typeof(resltris)}(undef, Γ4.patchnum)
    for idx in 1:1:Γ4.patchnum
        ltris_pat[idx] = []
    end
    #
    for (tri, pat) in zip(resltris, reslpats)
        push!(ltris_pat[pat], tri)
    end
    #
    return Gamma4(
        Γ4.model,
        Γ4.λ_0,
        Γ4.V,
        Γ4.k4tab,
        Γ4.patchnum,
        Γ4.patches,
        resltris, reslpats, nothing, ltris_pat
    )
end


"""
将Γ4中的ltris，靠近费米面的切割的更加细致
"""
function TFGamma4_reweight_ltris_mt(Γ4::Gamma4{T, P}, lval
    ; maxnum=5000000) where {T, P}
    #
    truncated = false
    if !isnothing(Γ4.ladjs)
        @warn "会丢失ladjs的信息"
    end
    if length(Γ4.ltris) >= maxnum
        @warn "已经超过最大数量"
        truncated = true
    end
    lamb = Γ4.λ_0 * exp(-lval)
    #切分离费米面近的
    ltris = Γ4.ltris
    newltris::Vector{P} = []
    for tri in ltris
        add_tri = false
        for didx in 1:1:Γ4.model.bandnum
            disp = Γ4.model.dispersion[didx]
            ver1 = tri.vertex[1]
            sgn1 = sign(disp(ver1.x, ver1.y))
            ver2 = tri.vertex[2]
            sgn2 = sign(disp(ver2.x, ver2.y))
            ver3 = tri.vertex[3]
            sgn3 = sign(disp(ver3.x, ver3.y))
            #
            engs = abs(disp(tri.center.x, tri.center.y)) / lamb
            #不在费米面上
            if sgn1 == sgn2 && sgn1 == sgn3
                #如果没有超过最大，或者还有明显贡献
                if !truncated || engs < 3
                    add_tri = true
                end
                #继续循环下一个能带看是否要refine
                continue
            end
            #如果没有超过最大，而且还不够精细
            if !truncated && engs > 0.0001
                ftris = split_triangle(tri)
                push!(newltris, ftris[1])
                push!(newltris, ftris[2])
                push!(newltris, ftris[3])
                push!(newltris, ftris[4])
                #如果在某个能带上已经切掉，就不再计算别的能带了
                add_tri = false
                break
            #费米面上的tri一定会被保留
            else
                add_tri = true
                #继续循环下一个能带看是否要refine
                continue
            end
        end
        if add_tri
            push!(newltris, tri)
        end
    end
    #
    refltris = newltris
    #重新计算一下所属的patch，可能会发生变化
    reflpats = group_ltris_into_patches_mt(refltris, Γ4.model.brillouin, Γ4.patchnum)
    #重新计算一下每个patch包含的三角形
    ltris_pat = Vector{typeof(refltris)}(undef, Γ4.patchnum)
    for idx in 1:1:Γ4.patchnum
        ltris_pat[idx] = []
    end
    #
    for (tri, pat) in zip(refltris, reflpats)
        push!(ltris_pat[pat], tri)
    end
    #
    return Gamma4(
        Γ4.model,
        Γ4.λ_0,
        Γ4.V,
        Γ4.k4tab,
        Γ4.patchnum,
        Γ4.patches,
        refltris, reflpats, nothing, ltris_pat
    )
end



"""
将等能面上面的一串三角形变成更接近的
"""
function refine_to_surface(ltris::Vector{P}, eng, disp) where P
    newltris::Vector{P} = []
    for tri in ltris
        ftris = split_triangle(tri)
        for ftri in ftris
            if @onsurface ftri disp eng
                push!(newltris, ftri)
            end
        end
    end
    return newltris
end


"""
将等能面上面的一串三角形变成更细的
"""
function refine_list_triangle(ltris::Vector{P}, itertime) where P
    nowltris = ltris
    for idx in 1:1:itertime
        newltris = Vector{P}(undef, 4*length(nowltris))
        for (idx, tri) in enumerate(nowltris)
            newltris[4*idx-3:4*idx] = split_triangle(tri)
        end
        nowltris = newltris
    end
    return nowltris
end


function TFGamma4_addition_ltris_mt(Γ4::Gamma4{T, P}, blval
    ; maxnum=5000000) where {T, P}
    if !isnothing(Γ4.ladjs)
        @warn "会丢失ladjs的信息"
    end
    if length(Γ4.ltris) >= maxnum
        @warn "已经超过最大数量"
    end
    lamb = Γ4.λ_0 * exp(-blval)
    #找到合适的有贡献的位置
    #println(eng)
    away::Vector{P} = []
    inner::Vector{P} = []
    for tri in Γ4.ltris
        isaway = true
        for bidx in 1:1:Γ4.model.bandnum
            disp = Γ4.model.dispersion[bidx]
            teng = disp(tri.center.x, tri.center.y)
            if abs(teng) <= lamb
                isaway = false
                break
            end
        end
        if isaway 
            push!(away, tri)
        else
            push!(inner, tri)
        end
    end
    #将有贡献的位置拆分的更细致
    iterf = log(maxnum - length(away)) - log(length(inner))
    iterf = iterf / log(4)
    iteri = floor(Int64, iterf)
    inner = refine_list_triangle(inner, iteri)
    #
    refltris = vcat(inner, away)
    #重新计算一下所属的patch，可能会发生变化
    reflpats = group_ltris_into_patches_mt(refltris, Γ4.model.brillouin, Γ4.patchnum)
    #重新计算一下每个patch包含的三角形
    ltris_pat = Vector{typeof(refltris)}(undef, Γ4.patchnum)
    for idx in 1:1:Γ4.patchnum
        ltris_pat[idx] = []
    end
    #
    for (tri, pat) in zip(refltris, reflpats)
        push!(ltris_pat[pat], tri)
    end
    #
    return Gamma4(
        Γ4.model,
        Γ4.λ_0,
        Γ4.V,
        Γ4.k4tab,
        Γ4.patchnum,
        Γ4.patches,
        refltris, reflpats, nothing, ltris_pat
    )
end


end # end module


