using Test
using MARY_fRG.Fermi


@testset "正方格子的功能" begin
    qs = common_square_lattice(0.)
    println(qs)
end


@testset "三角格子的功能" begin
    ts = common_triangle_lattice(0.)
    println(ts)
end
