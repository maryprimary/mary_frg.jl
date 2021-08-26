"""
测试Gamma4相关的功能
"""


using Test
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation
using MARY_fRG.Drawers
using MARY_fRG.Fermi.Surface


@testset "普通的正方格子" begin
    qs = common_square_lattice(0.)
    Γ4 = ECGamma4(
        qs, 8.0, 16, 50
    )
    plt = draw_lines(Γ4.model.brillouin.edges)
    sf, sfpidx = const_energy_line_in_patches(Γ4.ltris, Γ4.ladjs,
    Γ4.lpats, 0., Γ4.model.dispersion[1])
    draw_lines!(plt, sf; lcidx=sfpidx)
    patcs = Γ4.patches[1]
    draw_points!(plt, patcs; text=1:1:length(patcs))
    draw_polygon!(plt, Γ4.ltris; pcidx=Γ4.lpats)
    @savepng plt "Gamma4"
    @test Γ4.k4tab[1, 1, 1, 1, 1, 12, 12] == 1
end


@testset "普通的三角格子" begin
    ts = common_triangle_lattice(0.)
    Γ4 = ECGamma4(
        ts, 6.0, 18, 50
    )
    plt = draw_lines(Γ4.model.brillouin.edges)
    sf, sfpidx = const_energy_line_in_patches(Γ4.ltris, Γ4.ladjs,
    Γ4.lpats, 0., Γ4.model.dispersion[1])
    draw_lines!(plt, sf; lcidx=sfpidx)
    patcs = Γ4.patches[1]
    draw_points!(plt, patcs; text=1:1:length(patcs))
    draw_polygon!(plt, Γ4.ltris; pcidx=Γ4.lpats)
    @savepng plt "Gamma42"
    @test Γ4.k4tab[1, 1, 1, 1, 1, 12, 12] == 1
end


@testset "两个能带的正方格子" begin
    disp1(x, y) = -2*(cos(x) + cos(y))
    disp2(x, y) = -2*(cos(x) + cos(y)) + 1.0
    qs = QuadrateSystem{:DOUBLE}([disp1, disp2])
    Γ4 = ECGamma4(
        qs, 8.0, 16, 50
    )
    plt = draw_lines(Γ4.model.brillouin.edges)
    #
    sf, sfpidx = const_energy_line_in_patches(Γ4.ltris, Γ4.ladjs,
    Γ4.lpats, 0., Γ4.model.dispersion[1])
    draw_lines!(plt, sf; lcidx=sfpidx)
    #
    sf2, sf2pidx = const_energy_line_in_patches(Γ4.ltris, Γ4.ladjs,
    Γ4.lpats, 0., Γ4.model.dispersion[2])
    draw_lines!(plt, sf2; lcidx=sf2pidx)
    #
    patcs = Γ4.patches[1]
    draw_points!(plt, patcs; text=1:1:length(patcs))
    patcs2 = Γ4.patches[2]
    draw_points!(plt, patcs2; text=1:1:length(patcs))
    draw_polygon!(plt, Γ4.ltris; pcidx=Γ4.lpats)
    @savepng plt "Gamma43"
    @test Γ4.k4tab[1, 1, 2, 2, 1, 9, 1] == 9
end

