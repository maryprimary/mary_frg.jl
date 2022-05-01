"""
测试单带的正方格子
"""

using Plots
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation

include("../../src/helpers/checks.jl")

#阻止画图
ENV["GKSwstype"] = "100"


"""
从头开始运行代码
"""
function run_ec()
    model = common_square_lattice(0.20)
    Γ4 = ECGamma4CMPLX(
        model, 1.0, 16, 200
    )
    rotation_check2(Γ4.k4tab[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
    Γ4.V .+= 1.0 + 0.1im
    lval = 0.
    lstep = 0.01
    for idx in 1:1:701
        bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex = all_bubble_ec_mt(Γ4, lval; usesymm=false)
        bubb_pp2, bubb_fs2, bubb_nfs2, bubb_ex2, bubb_nex2 = all_bubble_ec_mt(Γ4, lval)
        rotation_check(bubb_pp.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        rotation_check(bubb_fs.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        rotation_check(bubb_ex.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        println(sum(abs.(bubb_pp.V - bubb_pp2.V)))
        println(sum(abs.(bubb_fs.V - bubb_fs2.V)))
        println(sum(abs.(bubb_ex.V - bubb_ex2.V)))
        dl = dl_ec_mt(Γ4, bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(real(Γ4.V[1, 1, 1, 1, :, :, 1]))
            png(plt, "Gamma4"*string(idx))
            plt = heatmap(imag(Γ4.V[1, 1, 1, 1, :, :, 1]))
            png(plt, "iGamma4"*string(idx))
        end
        lval += lstep
    end
end


#run_ec()


function run_tf()
    model = common_square_lattice(0.00)
    Γ4 = TFGamma4CMPLX(
        model, 1.0, 16, 100
    )
    rotation_check2(Γ4.k4tab[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
    Γ4.V .+= 1.0 + 0.1im
    lval = 0.
    lstep = 0.01
    for idx in 1:1:701
        bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, lval; usesymm=false)
        bubb_pp2, bubb_fs2, bubb_ex2 = all_bubble_tf_mt(Γ4, lval)
        rotation_check(bubb_pp.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        rotation_check(bubb_fs.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        rotation_check(bubb_ex.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        println(sum(abs.(bubb_pp.V - bubb_pp2.V)))
        println(sum(abs.(bubb_fs.V - bubb_fs2.V)))
        println(sum(abs.(bubb_ex.V - bubb_ex2.V)))
        dl = dl_tf_mt(Γ4, bubb_pp, bubb_fs, bubb_ex)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(real(Γ4.V[1, 1, 1, 1, :, :, 1]))
            png(plt, "Gamma4"*string(idx))
            plt = heatmap(imag(Γ4.V[1, 1, 1, 1, :, :, 1]))
            png(plt, "iGamma4"*string(idx))
        end
        lval += lstep
    end
end


run_tf()
