"""
测试温度流的refine
"""

using Plots
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation

#阻止画图
ENV["GKSwstype"] = "100"


function run_rftf()
    model = common_square_lattice(0.20)
    Γ4 = TFGamma4(
        model, 8.0, 16, 100
    )
    Γ4.V .+= 1.0
    lval = 0.
    lstep = 0.01
    #
    for idx in 1:1:2501
        #if idx % 10 == 1
        #    blval = lval - 5
        #    blval = max(blval, 0.)
        #    blval = min(blval, 5.)
        #    Γ4b = TFGamma4_addition_ltris_mt(Γ4, blval; maxnum=1000000)
        #    println(length(Γ4.ltris), "->", length(Γ4b.ltris))
        #    #plt = draw_points([tri.center for tri in Γ4b.ltris])
        #    #@savepng plt "brlu"*string(lval)
        #end
        #blval = min(lval, 10.0)
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
