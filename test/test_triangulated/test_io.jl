using Test
using MARY_fRG.Triangulated, MARY_fRG.Basics



@testset "保存相关的功能" begin
    squ1 = Square(Point2D(0., 0.), 2.)
    ltris, ladjs = split_square(squ1, 10)

    save_triangulated(squ1, ltris, ladjs, "squ1")

    @test_broken save_triangulated(ltris[1], ltris, ladjs, "tri1")

    hexa1 = EqHexagon(Point2D(0., 0.), 2.)
    ltris, ladjs = split_hexagon(hexa1, 10)
    save_triangulated(hexa1, ltris, ladjs, "hexa1")
end


@testset "读取相关的功能" begin
    #
    squ1, ltris, ladjs = load_triangulated("squ1")
    save_triangulated(squ1, ltris, ladjs, "squ2")
    orifile = open("squ1.tri", "r")
    cpyfile = open("squ2.tri", "r")
    @test read(orifile) == read(cpyfile)
    rm("squ2.tri")
    orifile = open("squ1.ply", "r")
    cpyfile = open("squ2.ply", "r")
    @test read(orifile) == read(cpyfile)
    rm("squ2.ply")
    orifile = open("squ1.adj", "r")
    cpyfile = open("squ2.adj", "r")
    @test read(orifile) == read(cpyfile)
    rm("squ2.adj")
    #
    hexa1, ltris, ladjs = load_triangulated("hexa1")
    save_triangulated(hexa1, ltris, ladjs, "hexa2")
    orifile = open("hexa1.tri", "r")
    cpyfile = open("hexa2.tri", "r")
    @test read(orifile) == read(cpyfile)
    rm("hexa2.tri")
    orifile = open("hexa1.ply", "r")
    cpyfile = open("hexa2.ply", "r")
    @test read(orifile) == read(cpyfile)
    rm("hexa2.ply")
    orifile = open("hexa1.adj", "r")
    cpyfile = open("hexa2.adj", "r")
    @test read(orifile) == read(cpyfile)
    rm("hexa2.adj")
end
