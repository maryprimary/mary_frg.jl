using Test
using MARY_fRG.Triangulated, MARY_fRG.Basics
using MARY_fRG.Drawers


#rts = Triangulated.split_single_square(Square(Point2D(0., 0.), 1.))
#
#plt = draw_figure(
#    [rt.center for rt in rts], nothing, rts;
#    ptext = 1:1:4
#)
#
#@savepng plt "sss"

@testset "正方形的切割" begin
    squ1 = Square(Point2D(0., 0.), 2.)
    ltris, ladjs = split_square(squ1, 10)
    #plt = draw_figure(
    #    nothing, nothing, ltris
    #)
    #@savepng plt "ss1"
    tri2idx = Dict(((tri, idx) for (idx, tri) in enumerate(ltris)))
    #pts = vcat([ltris[1].center; ladjs[1][1].center])#; ladjs[1][2].center; ladjs[1][3].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[1]], tri2idx[ladjs[1][1]]])
    #@savepng plt "ss2"
    @test tri2idx[ladjs[202][1]] == 102
    @test ismissing(ladjs[202][2])
    @test tri2idx[ladjs[202][3]] == 302
    @test tri2idx[ladjs[300][1]] == 200
    @test tri2idx[ladjs[300][2]] == 90
    @test tri2idx[ladjs[300][3]] == 400
    println(ladjs[310])
end
