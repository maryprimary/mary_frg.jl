"""
测试重新切分
"""

using Test
using MARY_fRG.Refine, MARY_fRG.Basics, MARY_fRG.Triangulated
using MARY_fRG.Drawers

@testset "三角形的切分" begin
    rtri = RtTriangle(Point2D(0., 0.), Point2D(1., 1.), Point2D(1., -1.))
    newtris = split_triangle(rtri)
    println(newtris)
    plt = draw_polygon(newtris)
    @savepng plt "newtris"
    eqtri = EqTriangle(Point2D(0., 0.), Point2D(1., sqrt(3)), Point2D(2., 0.))
    newtris = split_triangle(eqtri)
    println(newtris)
    plt = draw_polygon(newtris)
    @savepng plt "newtris2"
end


@testset "正方形的细化" begin
    squ1 = Square(Point2D(0., 0.), 1)
    ltris1, ladjs1 = split_square(squ1, 3)
    plt = draw_polygon(ltris1)
    @savepng plt "refine1"
    ltris2 = refine_triangles(ltris1)
    plt = draw_polygon(ltris2)
    @savepng plt "refine2"
end

@testset "六边型的细化" begin
    hexa1 = EqHexagon(Point2D(0., 0.), 1)
    ltris1, ladjs1 = split_hexagon(hexa1, 4)
    plt = draw_polygon(ltris1)
    @savepng plt "refine3"
    ltris2 = refine_triangles(ltris1)
    plt = draw_polygon(ltris2)
    @savepng plt "refine4"
end


@testset "寻找最近邻" begin
    #正方格子
    squ1 = Square(Point2D(0., 0.), 1.)
    ltris, ladjs = split_square(squ1, 11)
    ladjs2 = find_adjs_by_adjoint(ltris)
    save_triangulated(squ1, ltris, ladjs, "squ1")
    save_triangulated(squ1, ltris, ladjs2, "squ2")
    fl1 = open("squ1.adj", "r")
    fl2 = open("squ2.adj", "r")
    @test read(fl1) == read(fl2)
    #六角格子
    hexa1 = EqHexagon(Point2D(0., 0.), 1.)
    ltris, ladjs = split_hexagon(hexa1, 13)
    ladjs2 = find_adjs_by_adjoint(ltris)
    save_triangulated(squ1, ltris, ladjs, "hexa1")
    save_triangulated(squ1, ltris, ladjs2, "hexa2")
    fl1 = open("hexa1.adj", "r")
    fl2 = open("hexa2.adj", "r")
    @test read(fl1) == read(fl2)
    rm("hexa1.adj")
    rm("hexa2.adj")
    rm("hexa1.ply")
    rm("hexa2.ply")
    rm("hexa1.tri")
    rm("hexa2.tri")
    rm("squ1.adj")
    rm("squ2.adj")
    rm("squ1.ply")
    rm("squ2.ply")
    rm("squ1.tri")
    rm("squ2.tri")
end
