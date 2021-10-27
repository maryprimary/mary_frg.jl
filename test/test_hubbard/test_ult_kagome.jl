"""
测试上VHS的kagome
"""

using Plots
using MARY_fRG.Drawers
using MARY_fRG.Fermi.Kagome
using MARY_fRG.Fermi.Surface
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


function run_tf()
    model = upperband_kagome_lattice(0.)
    Γ4 = TFGamma4(
        model, 6.0, 24, 223
    )
    Γ4.V[1, 1, 1, 1, :, :, :] .+= get_kagome_U_mt(3.0, Γ4.patches[1])
    println(
        all(isfinite.(Γ4.V[1, 1, 1, 1, :, :, :]))
    )
    #plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
    #png(plt, "Gamma41")
    #return
    lval = 0.
    lstep = 0.01
    bpu, bfu, beu = get_ult_bubb(Γ4)
    @info "ult结束"
    for idx in 1:1:3501
        bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, lval)
        dl = dl_tf_mix_ult_mt(Γ4, bubb_pp, bubb_fs, bubb_ex, bpu, bfu, beu)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1 || idx > 1351
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma4"*string(idx))
        end
        lval += lstep
    end
end


#run_tf()


"""
绘制patch和费米面
"""
function draw_brlu()
    model = upperband_kagome_lattice(0.)
    Γ4 = ECGamma4(
        model, 6.0, 24, 223
    )
    plt = draw_lines(Γ4.model.brillouin.edges)
    suf = const_energy_line(Γ4.ltris, Γ4.ladjs, 0., Γ4.model.dispersion[1])
    draw_lines!(plt, suf; lcidx=ones(Int64, length(suf)))
    draw_points!(plt, Γ4.patches[1], text=1:1:24)
    @savepng plt "kag_pat_surf"
end



draw_brlu()

