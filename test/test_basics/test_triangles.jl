using Test
using MARY_fRG.Basics
using Printf


@testset "测试三角形" begin
    tri1 = Triangle([Point2D(1., 0.), Point2D(0., 1), Point2D(0., 0.)])
    println(tri1)
    println(Point2D(1., 0.))
    @test typeof(tri1) <: AbstractPolygon
    tri2 = Triangle(Point2D(1., 0.), Point2D(0., 1), Point2D(0., 0.))
    println(tri2)
    @test area(tri1) == 0.5
end


@testset "测试直角三角形" begin
    tri1 = RtTriangle([Point2D(0., 0.), Point2D(0., 1), Point2D(1., 0.)])
    @test area(tri1) == 0.5
    tri2 = RtTriangle(Point2D(1., 1.), Point2D(0., 0), Point2D(2., 0.))
    @test area(tri2) == 1.0
end


@testset "测试等边三角形" begin
    tri1 = EqTriangle([Point2D(0., 0.), Point2D(1., sqrt(3)), Point2D(2., 0.)])
    @test area(tri1) == sqrt(3)
    tri2 = EqTriangle(Point2D(1., 1.), Point2D(2., 1+sqrt(3)), Point2D(3., 1.))
    @test area(tri2) == sqrt(3)
end

