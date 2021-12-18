"""
测试kagome的相关功能
"""

using MARY_fRG.Basics
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Patch
using MARY_fRG.Fermi.Surface
using MARY_fRG.Drawers
using Plots
using MARY_fRG.FlowEquation
using LinearAlgebra

"""
运行
"""
function run()
    k1 = Point2D(-1.2258407235596014, -3.6275987284684357)
    k2 = Point2D(
        k1.x * cos(π/3) - k1.y * sin(π/3),
        k1.x * sin(π/3) + k1.y * cos(π/3)
    )
    println(k1)
    println(k2)
    println(k2.x/4 + √3*k2.y/4)
    println(k1.x/2)

    mat1 = zeros(3, 3)
    mat1[1, 2] = 2cos(k1.x/4 + √3*k1.y/4)
    mat1[1, 3] = 2cos(k1.x/2)
    mat1[2, 3] = 2cos(-k1.x/4 + √3*k1.y/4)
    mat1[2, 1] = 2cos(k1.x/4 + √3*k1.y/4)
    mat1[3, 1] = 2cos(k1.x/2)
    mat1[3, 2] = 2cos(-k1.x/4 + √3*k1.y/4)

    mat2 = zeros(3, 3)
    mat2[1, 2] = 2cos(k2.x/4 + √3*k2.y/4)
    mat2[1, 3] = 2cos(k2.x/2)
    mat2[2, 3] = 2cos(-k2.x/4 + √3*k2.y/4)
    mat2[2, 1] = 2cos(k2.x/4 + √3*k2.y/4)
    mat2[3, 1] = 2cos(k2.x/2)
    mat2[3, 2] = 2cos(-k2.x/4 + √3*k2.y/4)

    println(mat1)
    println(mat2)
    
    println(eigvecs(mat1))
    println(eigvecs(mat2))
    println(eigvecs(mat1)[:, 1])
    println(eigvecs(mat2)[:, 1])
    println(eigvecs(mat2)[:, 2])
    println(eigvecs(mat2)[:, 3])
    n1, n2, n3 = get_kagome_ν(k1.x, k1.y)
    println(n1)
    println(n2)
    println(n3)
    n1, n2, n3 = get_kagome_ν(k2.x, k2.y)
    println(n1)
    println(n2)
    println(n3)

    mat3 = zeros(ComplexF64, 3, 3)
    mat3[1, 2] = 1 + exp(-im*(k1.x/2+√3*k1.y/2))#2cos(2x)
    mat3[1, 3] = 1 + exp(-im*(k1.x))#exp(-im*(4x))#2cos(x+y)
    mat3[2, 3] = 1 + exp(-im*(k1.x/2 - √3*k1.y/2))#exp(im*(-2x+2y))#2cos(-x+y)
    mat3[2, 1] = 1 + exp(im*(k1.x/2+√3*k1.y/2))#2cos(2x)
    mat3[3, 1] = 1 + exp(im*(k1.x))#exp(im*(4x))#2cos(x+y)
    mat3[3, 2] = 1 + exp(im*(k1.x/2 - √3*k1.y/2))#exp(-im*(-2x+2y))#2cos(-x+y)
    println(eigvals(mat1))
    println(eigvals(mat2))
    println(eigvals(mat3))
end


#run()


function run2()
    model = wuxx_av3sb5_kagome_lattice()
    Γ4 = TFGamma4(
        model, 0.1, 144, 10
    )
    sitea = []
    siteb = []
    sitec = []
    for pat in Γ4.patches[1]
        nu1, nu2, nu3 = get_kagome_ν(pat.x, pat.y)
        push!(sitea, nu2[1])
        push!(siteb, nu2[2])
        push!(sitec, nu2[3])
    end
    plt = plot(1:1:length(sitea), sitea)
    plot!(plt, 1:1:length(sitea), siteb)
    plot!(plt, 1:1:length(sitea), sitec)
    png(plt, "kag_band_contrib4")
end

run2()
