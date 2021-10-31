"""
测试单带三角格子，tJ
"""

using Plots
using MARY_fRG.Interactions
using MARY_fRG.Fermi, MARY_fRG.Fermi.Surface
using MARY_fRG.FlowEquation


#阻止画图
ENV["GKSwstype"] = "100"


function run()
    model = common_triangle_lattice(-2.0)
    Γ4 = ECGamma4(
        model, 6.85, 18, 50
    )
    Γ4.V[1, 1, 1, 1, :, :, :] .+= triangular_system_heisenberg(
        Γ4.patches[1], Γ4.model.kadd, 1.0
    )
    lval = 0.
    lstep = 0.01
    for idx = 1:1:301
        bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex = all_bubble_ec_mt(Γ4, lval)
        dl = dl_ec_mt(Γ4, bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma4"*string(idx))
        end
        lval += lstep
    end
end


run()
