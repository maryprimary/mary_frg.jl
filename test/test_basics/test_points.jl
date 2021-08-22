using Test
using MARY_fRG.Basics


@testset "点相关的测试" begin
    pt1 = Point2D(1., 1.)
    println(pt1)
    pt2 = Point2D(2., 3.)
    pt3 = pt1 + pt2
    @test (pt3.y == 4) && (pt3.x == 3)
    pt4 = middle_point(pt1, pt3)
    @test (pt4.x == 2) && (pt4.y == 2.5)
    ang1 = absolute_angle(pt1)
    @test ang1 == pi / 4
    ang1 = absolute_angle(Point2D(-1., 1.))
    @test ang1 == 3pi / 4
    ang1 = absolute_angle(Point2D(1., -1.))
    @test ang1 == 7pi / 4
    pt5 = -pt1
    println(pt5)
    pt6 = pt5 - pt1
    println(pt6)
    pt7 = pt5 - 2pt1
    println(pt7)
end
