"""
测试上wuxx的 A V_3 Sb_5
"""

using Plots
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Surface
using MARY_fRG.FlowEquation


#阻止画图
ENV["GKSwstype"] = "100"


"""
获取极低温的时候的bubble
"""
function get_ult_bubb(Γ4; maxnum=2000000)
    maxi, mini = engpeak_to_surface(Γ4)
    println(maxi, " ", mini)
    Γ4_ult = fliter_away_surface(Γ4, 2.0)
    Γ4_ult = refine_to_surface(Γ4_ult)
    while length(Γ4_ult.ltris) < maxnum
        Γ4_ult = refine_to_surface(Γ4_ult)
    end
    plt = draw_points([tri.center for tri in Γ4_ult.ltris])
    @savepng plt "ult_surface"
    maxi, mini = engpeak_to_surface(Γ4_ult)
    println(maxi, " ", mini)
    bubb_pp_u, bubb_fs_u, bubb_ex_u = all_bubble_tf_ult_mt(Γ4_ult, maxi)
    return bubb_pp_u, bubb_fs_u, bubb_ex_u
end


function run_tf()
    model = wuxx_av3sb5_kagome_lattice()
    Γ4 = TFGamma4(
        model, 6.0, 24, 233
    )
    Γ4.V .+= get_wuxx_U_mt(4.0, 2.0, 2.0, Γ4.patches[1], Γ4.patches[2])
    println(
        all(isfinite.(Γ4.V))
    )
    plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
    png(plt, "Gamma41")
    #return
    lval = 0.
    lstep = 0.01
    bpu, bfu, beu = get_ult_bubb(Γ4)
    #
    #plt = heatmap(1e5*bpu.V[1, 1, 1, 1, 1, :, :])
    #png(plt, "wuxx_bpu1")
    #plt = heatmap(1e5*bpu.V[2, 2, 2, 2, 1, :, :])
    #png(plt, "wuxx_bpu2")
    #plt = heatmap(1e5*bpu.V[1, 2, 2, 1, 1, :, :])
    #png(plt, "wuxx_bpu3")
    #plt = heatmap(1e5*bpu.V[2, 1, 1, 2, 1, :, :])
    #png(plt, "wuxx_bpu4")
    ##
    #plt = heatmap(1e5*bfu.V[1, 1, 1, 1, 1, :, :])
    #png(plt, "wuxx_bfu1")
    #plt = heatmap(1e5*bfu.V[2, 2, 2, 2, 1, :, :])
    #png(plt, "wuxx_bfu2")
    #plt = heatmap(1e5*bfu.V[1, 2, 2, 1, 1, :, :])
    #png(plt, "wuxx_bfu3")
    #plt = heatmap(1e5*bfu.V[2, 1, 1, 2, 1, :, :])
    #png(plt, "wuxx_bfu4")
    ##
    #plt = heatmap(1e5*beu.V[1, 1, 1, 1, 1, :, :])
    #png(plt, "wuxx_beu1")
    #plt = heatmap(1e5*beu.V[2, 2, 2, 2, 1, :, :])
    #png(plt, "wuxx_beu2")
    #plt = heatmap(1e5*beu.V[1, 2, 2, 1, 1, :, :])
    #png(plt, "wuxx_beu3")
    #plt = heatmap(1e5*beu.V[2, 1, 1, 2, 1, :, :])
    #png(plt, "wuxx_beu4")
    @info "ult结束"
    for idx in 1:1:5001
        bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, lval)
        dl = dl_tf_mix_ult_mt(Γ4, bubb_pp, bubb_fs, bubb_ex, bpu, bfu, beu)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma41111"*string(idx))
            plt = heatmap(Γ4.V[2, 2, 2, 2, :, :, 1])
            png(plt, "Gamma42222"*string(idx))
            plt = heatmap(Γ4.V[1, 2, 2, 1, :, :, 1])
            png(plt, "Gamma41221"*string(idx))
            plt = heatmap(Γ4.V[2, 1, 1, 2, :, :, 1])
            png(plt, "Gamma42112"*string(idx))
        end
        lval += lstep
    end
end


run_tf()


