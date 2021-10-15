using Test
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Surface


@testset "正方格子的功能" begin
    qs = common_square_lattice(0.)
    println(qs)
end


@testset "三角格子的功能" begin
    ts = common_triangle_lattice(0.)
    println(ts)
end


@testset "kagome格子的功能" begin
    ks = upperband_kagome_lattice(0.)
    plt = draw_lines(ks.brillouin.edges)
    @savepng plt "test_latt_ks1"
end
