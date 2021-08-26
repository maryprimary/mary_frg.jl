using Test
using MARY_fRG.Triangulated
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Surface
using MARY_fRG.Drawers

@testset "正方格子的费米面" begin
    qs = common_square_lattice(0.)
    ltris, ladjs = split_square(qs.brillouin, 100)
    sf = const_energy_line(ltris, ladjs, 0., qs.dispersion[1])
    plt = draw_lines(qs.brillouin.edges)
    draw_lines!(plt, sf; lcidx=ones(Int64, length(sf)))
    #draw_polygon!(plt, ltris)
    sf2 = const_energy_line(ltris, ladjs, -0.2, qs.dispersion[1])
    draw_lines!(plt, sf2; lcidx=2ones(Int64, length(sf2)))
    @savepng plt "surface1"
end


@testset "三角格子的费米面" begin
    ts = common_triangle_lattice(0.)
    ltris, ladjs = split_hexagon(ts.brillouin, 100)
    sf = const_energy_line(ltris, ladjs, 0., ts.dispersion[1])
    plt = draw_lines(ts.brillouin.edges)
    draw_lines!(plt, sf; lcidx=ones(Int64, length(sf)))
    sf2 = const_energy_line(ltris, ladjs, 2.0, ts.dispersion[1])
    draw_lines!(plt, sf2; lcidx=2ones(Int64, length(sf2)))
    @savepng plt "surface2"
end
