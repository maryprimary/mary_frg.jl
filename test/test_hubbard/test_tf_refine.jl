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
        model, 8.0, 16, 1000
    )
    Γ4.V .+= 1.0
    lval = 0.
    lstep = 0.01
    #plt = draw_points([tri.center for tri in Γ4.ltris])
    #draw_polygon(Γ4.ltris; pcidx=Γ4.lpats)
    #@savepng plt "brlu"*string(lval)
    for idx in 1:1:1501
        if idx % 10 == 1
            oldΓ4 = Γ4
            Γ4 = TFGamma4_refine_ltris_mt(Γ4, lval)
            if !(oldΓ4 === Γ4)
                println("refine")
                println(length(Γ4.ltris))
                #plt = draw_points([tri.center for tri in Γ4.ltris])
                #@savepng plt "brlu"*string(lval)
            end
        end
        #bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, lval)
        #dl = dl_tf_mt(Γ4, bubb_pp, bubb_fs, bubb_ex)
        #Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma4"*string(idx))
        end
        lval += lstep
    end
end


run_rftf()
