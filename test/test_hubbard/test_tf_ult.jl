"""
测试温度流的refine
"""

using Plots
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation

#阻止画图
ENV["GKSwstype"] = "100"

"""
获取极低温的时候的bubble
"""
function get_ult_bubb(Γ4; maxnum=1000000)
    Γ4_ult = fliter_away_surface(Γ4, 2.0)
    Γ4_ult = refine_to_surface(Γ4_ult)
    while length(Γ4_ult.ltris) < maxnum
        Γ4_ult = refine_to_surface(Γ4_ult)
    end
    plt = draw_points([tri.center for tri in Γ4_ult.ltris])
    @savepng plt "ult_surface"
    maxi, mini = engpeak_to_surface(Γ4_ult)
    bubb_pp_u, bubb_fs_u, bubb_ex_u = all_bubble_tf_ult_mt(Γ4_ult, maxi)
    return bubb_pp_u, bubb_fs_u, bubb_ex_u
end



function run_rftf()
    model = common_square_lattice(0.20)
    Γ4 = TFGamma4(
        model, 8.0, 16, 200
    )
    Γ4.V .+= 1.0
    lval = 0.
    lstep = 0.01
    #plt = draw_points([tri.center for tri in Γ4.ltris])
    #@savepng plt "brlu"*string(lval)
    #
    bpu, bfu, beu = get_ult_bubb(Γ4)
    @info "ult结束"
    #
    for idx in 1:1:3351
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
        dl = dl_tf_mix_ult_mt(Γ4, bubb_pp, bubb_fs, bubb_ex, bpu, bfu, beu)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1 || idx > 3301
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma4"*string(idx))
        end
        lval += lstep
    end
end


run_rftf()
