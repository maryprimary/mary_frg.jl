"""
测试单带三角格子，tJ
"""

using Plots
using MARY_fRG.Drawers
using MARY_fRG.Interactions
using MARY_fRG.Fermi, MARY_fRG.Fermi.Surface
using MARY_fRG.FlowEquation

include("../../src/helpers/checks.jl")

#阻止画图
ENV["GKSwstype"] = "100"


function run()
    model = common_triangle_lattice(-2.0)
    Γ4 = ECGamma4(
        model, 2.85, 24, 500
    )
    rotation_check2(Γ4.k4tab[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
    plt = draw_lines(model.brillouin.edges)
    println(Γ4.patches[1][20], Γ4.patches[1][12], Γ4.patches[1][15])
    k4 = Fermi.hexagon_kadd2(Γ4.patches[1][20], Γ4.patches[1][12])
    k4 = Fermi.hexagon_kadd2(k4, -Γ4.patches[1][15])
    println(k4, Γ4.k4tab[1, 1, 1, 1, 20, 12, 15])
    println(find_patch_index_hexa(k4, model.brillouin, 24))
    scatter!(plt, [Γ4.patches[1][20].x, Γ4.patches[1][12].x, Γ4.patches[1][15].x],
    [Γ4.patches[1][20].y, Γ4.patches[1][12].y, Γ4.patches[1][15].y], markercolor=:red, markersize=10
    )
    scatter!(plt, [k4.x], [k4.y], markercolor=:yellow, markersize=10)
    println(Γ4.patches[1][24], Γ4.patches[1][16], Γ4.patches[1][19])
    k4 = Fermi.hexagon_kadd2(Γ4.patches[1][24], Γ4.patches[1][16])
    k4 = Fermi.hexagon_kadd2(k4, -Γ4.patches[1][19])
    println(k4, Γ4.k4tab[1, 1, 1, 1, 24, 16, 19])
    scatter!(plt, [Γ4.patches[1][24].x, Γ4.patches[1][16].x, Γ4.patches[1][19].x],
    [Γ4.patches[1][24].y, Γ4.patches[1][16].y, Γ4.patches[1][19].y], markercolor=:green,  markersize=5
    )
    scatter!(plt, [k4.x], [k4.y], markercolor=:blue, markersize=5)
    png(plt, "testk4")
    Γ4.V[1, 1, 1, 1, :, :, :] .+= triangular_system_heisenberg(
        Γ4.model, Γ4.patches[1], 1.0
    )
    lval = 0.
    lstep = 0.01
    for idx = 1:1:301
        bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex = all_bubble_ec_mt(Γ4, lval)
        bubb_pp2, bubb_fs2, bubb_nfs2, bubb_ex2, bubb_nex2 = all_bubble_ec_mt(Γ4, lval; usesymm=false)
        println(sum(abs.(bubb_pp.V - bubb_pp2.V)))
        println(sum(abs.(bubb_fs.V - bubb_fs2.V)))
        println(sum(abs.(bubb_ex.V - bubb_ex2.V)))
        rotation_check(bubb_pp2.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        rotation_check(bubb_fs2.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
        rotation_check(bubb_ex2.V[1, 1, 1, 1, :, :, :], model, Γ4.patches[1], Γ4.patches[1])
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
