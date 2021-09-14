"""
测试温度流的refine
"""

using Plots
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation


function run_rftf()
    model = common_square_lattice(0.20)
    Γ4 = TFGamma4(
        model, 8.0, 16, 200
    )
    Γ4.V .+= 2.0
    lval = 0.
    lstep = 0.01
    plt = draw_points([tri.center for tri in Γ4.ltris])
    @savepng plt "brlu"*string(lval)
    for idx in 1:1:2501
        if idx % 10 == 1
            newΓ4 = TFGamma4_refine_ltris_mt(Γ4, lval)
            if length(newΓ4.ltris) != length(Γ4.ltris)
                plt = draw_points([tri.center for tri in newΓ4.ltris])
                @savepng plt "brlu"*string(lval)
            end
            ##重复到不需要再refine
            #while !(newΓ4 === Γ4)
            #    Γ4 = newΓ4
            #    newΓ4 = TFGamma4_refine_ltris_mt(Γ4, lval)
            #    plt = draw_points([tri.center for tri in Γ4.ltris])
            #    @savepng plt "brlu"*string(lval)
            #end
            Γ4 = newΓ4
        end
        bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, lval)
        dl = dl_tf_mt(Γ4, bubb_pp, bubb_fs, bubb_ex)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma4"*string(idx))
        end
        lval += lstep
    end
end


run_rftf()
