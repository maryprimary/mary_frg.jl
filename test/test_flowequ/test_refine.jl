using Test

using MARY_fRG.Fermi
using MARY_fRG.FlowEquation
using MARY_fRG.Drawers


@testset "测试精细化" begin
    qs = common_triangle_lattice(0.11)
    Γ4 = TFGamma4(qs, 6.0, 24, 10)
    plt = draw_lines(Γ4.model.brillouin.edges)
    draw_polygon!(plt, Γ4.ltris; pcidx=Γ4.lpats)
    @savepng plt "refine1"
    #
    for idx in 2:1:9
        newΓ4 = TFGamma4_refine_ltris_mt(Γ4, 5.0)
        println(idx, ": ", length(newΓ4.ltris))
        plt = draw_lines(newΓ4.model.brillouin.edges)
        draw_polygon!(plt, newΓ4.ltris; pcidx=newΓ4.lpats)
        @savepng plt "refine"*string(idx)
        Γ4 = newΓ4
    end
end

