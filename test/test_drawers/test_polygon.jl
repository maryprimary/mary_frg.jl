"""测试多边形的所有功能"""

using Test
using MARY_fRG.Basics, MARY_fRG.Drawers


@testset "点的绘制" begin
    plt = draw_points([Point2D(0., 0.), Point2D(1., 2.)])
    @savepng plt "pts1"
    draw_points!(plt, [Point2D(3., 9.)]; text=["1"])
    @savepng plt "pts2"
end


@testset "线的绘制" begin
    plt = draw_lines([Segment(Point2D(0., 0.), Point2D(1., 2.))])
    @savepng plt "lns1"
    draw_lines!(plt, [
        Segment(Point2D(3., 3.), Point2D(-1., -9.)),
        Segment(Point2D(-3., 3.), Point2D(-1., 9.))
    ])
    @savepng plt "lns2"
end


@testset "多边形的绘制" begin
    plt = draw_polygon([
        Square(Point2D(0., 0.), 1.)
    ])
    @savepng plt "pls1"
    draw_polygon!(plt,
    [
        Square(Point2D(1., 1.), 1.),
        EqHexagon(Point2D(3., 3.), 1)
    ]; pcidx = [4, 7])
    @savepng plt "pls2"
end


@testset "图形的绘制" begin
    tri = RtTriangle(Point2D(2., 2.), Point2D(1., 1.), Point2D(3., 1.))
    squ = Square(Point2D(0., 0.), 2.5)
    hexa = EqHexagon(Point2D(-2., -2.), 3.)
    plt1 = draw_figure(hexa.vertex, hexa.edges, [hexa])
    @savepng plt1 "fig1"
    plt2 = draw_figure(
        vcat([tri.vertex; squ.vertex; hexa.vertex]),
        nothing, nothing
    )
    @savepng plt2 "fig2"
    plt3 = draw_figure(
        nothing,
        vcat([tri.edges; squ.edges; hexa.edges]),
        [tri, squ, hexa];
        lcidx = ones(Int64, 13),
        pcidx = 1:1:3
    )
    @savepng plt3 "fig3"
end
