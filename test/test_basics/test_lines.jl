using Test
using MARY_fRG.Basics


@testset "线段相关的内容" begin
    seg1 = Segment(Point2D(1.0, 1.0), Point2D(-2, 7.3))
    println(seg1)
    @test seg1.length == 6.977822009767804
end
