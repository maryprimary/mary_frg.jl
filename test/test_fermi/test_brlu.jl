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
    for idx in 1:1:1000000
        pt1 = Point2D((rand()-0.5)*4π, (rand()-0.5)*4π)
        pt2 = Point2D((rand()-0.5)*4π, (rand()-0.5)*4π)
        pt3 = Fermi.hexagon_kadd(pt1, pt2)
        pt4 = Fermi.hexagon_kadd2(pt1, pt2)
        pt5 = Fermi.hexagon_kadd3(pt1, pt2)
        #println(pt3)
        #println(pt4)
        #if pt3.x != pt4.x || pt3.y != pt4.y
        #    println(pt1)
        #    println(pt2)
        #end
        @test isapprox(pt3.x, pt4.x, atol=1e-12)
        @test isapprox(pt3.y, pt4.y, atol=1e-12)
        @test isapprox(pt3.x, pt5.x, atol=1e-12)
        @test isapprox(pt3.y, pt5.y, atol=1e-12)
    end
    #println(Fermi.hexagon_kadd2(
    #    Point2D(0.8806377184930461, 5.305646990668821),
    #    Point2D(-6.126445471920293, -6.016421208448352)
    #))
end
