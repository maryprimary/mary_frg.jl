"""
测试A V_3 Sb_5的基本功能
"""


using Plots
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Patch
using MARY_fRG.Fermi.Surface
using MARY_fRG.FlowEquation
using MARY_fRG.Drawers

"""
绘制费米面
"""
function draw_surface()
    ws = wuxx_av3sb5_kagome_lattice()
    Γ4 = ECGamma4(ws, 6.0, 24, 100)
    #pats1 = patches_under_vonhove(ws.brillouin, ws.dispersion[1], 24)
    #pats2 = patches_under_vonhove(ws.brillouin, ws.dispersion[2], 24)
    #
    plt = draw_lines(ws.brillouin.edges)
    draw_points!(plt, Γ4.patches[1])
    draw_points!(plt, Γ4.patches[2])
    #
    suf1, cidx1 = const_energy_line_in_patches(
        Γ4.model, Γ4.ltris, Γ4.ladjs, Γ4.lpats, 0., 1
    )
    draw_lines!(plt, suf1; lcidx=cidx1)
    suf2, cidx2 = const_energy_line_in_patches(
        Γ4.model, Γ4.ltris, Γ4.ladjs, Γ4.lpats, 0., 2
    )
    draw_lines!(plt, suf2; lcidx=cidx2)
    draw_polygon!(plt, Γ4.ltris, pcidx=Γ4.lpats)
    @savepng plt "av3sb5_suface"
end

draw_surface()


"""
相互作用
"""
function draw_interaction()
    ws = wuxx_av3sb5_kagome_lattice()
    Γ4 = ECGamma4(ws, 6.0, 24, 100)
    Γ4.V .+= get_wuxx_U_mt(4.0, 2.0, 2.0, Γ4.patches[1], Γ4.patches[2])
    println(Γ4.V[1, 1, 1, 1, 1, 2, 2], " ", Γ4.V[1, 1, 1, 1, 2, 1, 1])
    println(Γ4.V[1, 2, 2, 1, 1, 2, 2], " ", Γ4.V[2, 1, 1, 2, 2, 1, 1])
    println(Γ4.V[1, 2, 2, 1, 11, 7, 11], " ", Γ4.V[2, 1, 1, 2, 7, 11, 7])
    println(Γ4.V[1, 2, 2, 1, 1, 13, 2], " ", Γ4.V[1, 2, 2, 1, 14, 2, 13])
    #println(Γ4.V[1, 2, 2, 1, 15, 7, 11], " ", Γ4.V[1, 2, 2, 1, 15, 11, 7])
    #println(Γ4.V[2, 1, 1, 2, 15, 7, 11], " ", Γ4.V[2, 1, 1, 2, 15, 11, 7])
    plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
    @savepng plt "av3sb5_int1"
    plt = heatmap(Γ4.V[2, 2, 2, 2, :, :, 1])
    @savepng plt "av3sb5_int2"
    plt = heatmap(Γ4.V[1, 2, 2, 1, 1, :, :])
    @savepng plt "av3sb5_int3"
    plt = heatmap(Γ4.V[2, 1, 1, 2, :, :, 1])
    @savepng plt "av3sb5_int4"
end

draw_interaction()
