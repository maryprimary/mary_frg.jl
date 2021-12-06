"""
验证一些值是否正确
"""


using JLD
using Plots
using LinearAlgebra
using Printf
using MARY_fRG.Basics
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Patch
using MARY_fRG.FlowEquation


"""
没有动量转移的通道和PI通道
"""
function channel_noQ(model::P, V,
    kpats::Vector{Point2D}, ppats::Vector{Point2D}) where P <: Fermi.Abstract2DModel
    Vssc = zeros(length(kpats), length(ppats))
    Vtsc = zeros(length(kpats), length(ppats))
    Vspi = zeros(length(kpats), length(ppats))
    Vtpi = zeros(length(kpats), length(ppats))
    #寻找patch编号
    findpidx = isa(model, QuadrateSystem) ? find_patch_index_squa :
    find_patch_index_hexa
    #
    for cidx in CartesianIndices(Vssc)
        kidx, pidx = Tuple(cidx)
        nkidx = findpidx(-kpats[kidx], model.brillouin, length(kpats))
        npidx = findpidx(-ppats[pidx], model.brillouin, length(ppats))
        Vssc[kidx, pidx] = V[kidx, nkidx, pidx] + V[kidx, nkidx, npidx] +
        V[nkidx, kidx, pidx] + V[nkidx, kidx, npidx]
        Vtsc[kidx, pidx] = -V[kidx, nkidx, pidx] + V[kidx, nkidx, npidx] +
        V[nkidx, kidx, pidx] - V[nkidx, kidx, npidx]
        Vspi[kidx, pidx] = V[kidx, pidx, pidx] - 0.5*V[kidx, pidx, kidx]
        Vtpi[kidx, pidx] = -0.5*V[kidx, pidx, kidx]
    end
    return Vssc, Vtsc, Vspi, Vtpi
end



"""
验证是否满足旋转对称性
"""
function rotation_check(V, model, pats1, pats2)
    Vssc, Vtsc, Vspi, Vtpi = channel_noQ(model, V, pats1, pats2)
    #验证6度对称性
    gap = Int64(length(pats1) // 6)
    for idx in 1:1:length(pats1)-gap
        if !isapprox(Vssc[idx, idx], Vssc[idx+gap, idx+gap])
            println("Vssc ", idx, " != ", idx+gap)
        #else
            println(Vssc[idx, idx], " ", Vssc[idx+gap, idx+gap])
        end
    end
    #验证6度对称性
    gap = Int64(length(pats1) // 6)
    for idx1 in 1:1:length(pats1)-gap; for idx2 in 1:1:length(pats1)-gap; for idx3 in 1:1:length(pats1)-gap
        if !isapprox(V[idx1, idx2, idx3], V[idx1+gap, idx2+gap, idx3+gap])
            println("V[", idx1, ",", idx2, ",", idx3, "] != V[", idx1+gap, ",", idx2+gap, ",", idx3+gap, "]")
            println(V[idx1, idx2, idx3], " ", V[idx1+gap, idx2+gap, idx3+gap])
        end
    end; end; end
    println("rotation_check finish")
end

