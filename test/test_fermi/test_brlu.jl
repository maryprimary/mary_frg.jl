"""
测试正方晶系
"""

using Test
using MARY_fRG.Basics, MARY_fRG.Fermi


@testset "正方晶系" begin
    disp(x, y) = -2 * (cos(x) + cos(y))
    ss = QuadrateSystem{:COM}([disp])
    println(ss)
    println(ss.brillouin.vertex)
    pt1 = Point2D(2., 2.)
    pt2 = ss.kadd(pt1, pt1)
    @test pt2.x == 4 - 2*pi
    pt3 = ss.kadd(pt2, -pt1)
    @test pt3.x == 2
end


@testset "三角晶系" begin
    disp(x, y) = -2*numpy.cos(kxv) - 4*numpy.cos(kyp)*numpy.cos(kxv/2)
    ts = TriangularSystem{:COM}([disp])
    println(ts)
    println(ts.brillouin.vertex)
    pt1 = Point2D(2., 2.)
    pt2 = ts.kadd(pt1, pt1)
    @test pt2.y == 4 - 2pi/sqrt(3)
    @test pt2.x == 4 - 2pi
    pt3 = Point2D(2pi, 0.)
    pt3 = ts.kadd(pt2, pt3)
    @test pt3.y == 4 - 2pi/sqrt(3) - 2pi/sqrt(3)
    @test pt3.x == 4 + 2pi - 4pi
end
