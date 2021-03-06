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
    Vssc = zeros(ComplexF64, length(kpats), length(ppats))
    Vtsc = zeros(ComplexF64, length(kpats), length(ppats))
    Vspi = zeros(ComplexF64, length(kpats), length(ppats))
    Vtpi = zeros(ComplexF64, length(kpats), length(ppats))
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
function rotation_check(V, model, pats1, pats2; sailence=false)
    Vssc, Vtsc, Vspi, Vtpi = channel_noQ(model, V, pats1, pats2)
    symm = true
    #验证
    if isa(model, Fermi.TriangularSystem)
        gap = Int64(length(pats1) // 6)
    else
        gap = Int64(length(pats1) // 4)
    end
    for idx in 1:1:length(pats1)-gap
        if !isapprox(Vssc[idx, idx], Vssc[idx+gap, idx+gap])
            symm = false
            if !sailence
                println("Vssc ", idx, " != ", idx+gap)
                println(Vssc[idx, idx], " ", Vssc[idx+gap, idx+gap])
            end
        end
    end
    #验证6度对称性
    for idx1 in 1:1:length(pats1)-gap; for idx2 in 1:1:length(pats1)-gap; for idx3 in 1:1:length(pats1)-gap
        if !isapprox(V[idx1, idx2, idx3], V[idx1+gap, idx2+gap, idx3+gap])
            symm = false
            if !sailence
                println("V[", idx1, ",", idx2, ",", idx3, "] != V[", idx1+gap, ",", idx2+gap, ",", idx3+gap, "]")
                println(V[idx1, idx2, idx3], " ", V[idx1+gap, idx2+gap, idx3+gap])
            end
        end
    end; end; end
    println("rotation_check finish")
    return symm
end


"""
验证是否满足旋转对称性
"""
function rotation_check2(V, model, pats1, pats2)
    #验证
    if isa(model, Fermi.TriangularSystem)
        gap = Int64(length(pats1) // 6)
    else
        gap = Int64(length(pats1) // 4)
    end
    for idx1 in 1:1:length(pats1)-gap; for idx2 in 1:1:length(pats1)-gap; for idx3 in 1:1:length(pats1)-gap
        gdiff = V[idx1+gap, idx2+gap, idx3+gap] - V[idx1, idx2, idx3]
        if gdiff != gap && gdiff != gap - length(pats1)
            println("V[", idx1, ",", idx2, ",", idx3, "] != V[", idx1+gap, ",", idx2+gap, ",", idx3+gap, "]")
        #else
            println(V[idx1, idx2, idx3], " ", V[idx1+gap, idx2+gap, idx3+gap])
        end
    end; end; end
    println("rotation_check finish")
end


"""
验证反转指标的对称性
"""
function reflection_check(V, k4tab)
    for idxs in CartesianIndices(k4tab)
        k1i, k2i, k3i = Tuple(idxs)
        k4i = k4tab[k1i, k2i, k3i]
        if !isapprox(V[k1i, k2i, k3i], V[k4i, k3i, k2i])
            println(k1i, " ", k2i, " ", k3i, " ", k4i, " not full inverse")
        end
        if !isapprox(V[k1i, k2i, k3i], V[k2i, k1i, k4i])
            println(k1i, " ", k2i, " ", k3i, " ", k4i, " not anti inverse")
            println(k1i, " ", k2i, " ", k3i, " ", V[k1i, k2i, k3i], 
            " ", k2i, " ", k1i, " ", k4i, " ", V[k2i, k1i, k4i])
        end
    end
    println("reflection_check finish")
end


"""
验证k4tab的对称性
"""
function reflection_check2(k4tab)
    npat = size(k4tab)[1]
    for idxs in CartesianIndices(k4tab)
        k1i, k2i, k3i = Tuple(idxs)
        k4i = k4tab[k1i, k2i, k3i]
        if k4i == -1
            #println(k1i, " ", k2i, " ", k3i, " ", k4i)
            continue
        #else
        #    println(k1i, " ", k2i, " ", k3i, " ", k4i)
        end
        #println(k1i, " ", k2i, " ", k3i, " ", k4i)
        k1i2 = k4tab[k4i, k3i, k2i]
        if k1i2 != k1i
            println(k1i, " ", k2i, " ", k3i, " ", k4i, " not full inverse ", k4i, " ", k3i, " ", k2i, " ", k1i2)
        end
        k3i2 = k4tab[k2i, k1i, k4i]
        if k3i2 != k3i
            println(k1i, " ", k2i, " ", k3i, " ", k4i, " not anti inverse ", k2i, " ", k1i, " ", k4i, " ", k3i2)
        end
    end
    println("reflection_check2 finish")
end
