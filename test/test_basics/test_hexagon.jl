using Test
using MARY_fRG.Basics


@testset "等边六边型" begin
    hexa1 = EqHexagon(Point2D(0., 0.), 1)
    println(hexa1)
    hexa2 = EqHexagon(Point2D(-2., 3.), 1)
    @test abs(area(hexa2) - 3*sqrt(3) / 2) < 1e-6
end
