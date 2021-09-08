using Test

using MARY_fRG.Fermi
using MARY_fRG.FlowEquation
using MARY_fRG.Drawers


@testset "切分费米面" begin
    qs = common_square_lattice(0.20)
    Γ4 = TFGamma4(qs, 8.0, 16, 4)
    plt = draw_lines(Γ4.model.brillouin.edges)
    draw_polygon!(plt, Γ4.ltris; pcidx=Γ4.lpats)
    @savepng plt "brlu-20"
    newΓ4 = TFGamma4_reweight_ltris_mt(Γ4, 0.0, 1.0, 4)
    plt = draw_lines(newΓ4.model.brillouin.edges)
    draw_polygon!(plt, newΓ4.ltris; pcidx=newΓ4.lpats)
    @savepng plt "brlu-201"
end
