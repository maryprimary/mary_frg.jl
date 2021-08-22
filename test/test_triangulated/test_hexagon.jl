using Test
using MARY_fRG.Drawers, MARY_fRG.Triangulated, MARY_fRG.Basics



@testset "六边型的切割" begin
    hexa1 = EqHexagon(Point2D(0., 0.), 2)
    ltris, ladjs = split_hexagon(hexa1, 10)
    #plt = draw_figure(
    #    nothing, hexa1.edges, nothing
    #)
    #@savepng plt "sh1"
    #
    tri2idx = Dict(((tri, idx) for (idx, tri) in enumerate(ltris)))
    #pts = vcat([ltris[1].center; ladjs[1][2].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[1]], tri2idx[ladjs[1][2]]])
    #draw_lines!(plt, [ltris[1].edges[2]])
    #
    #pts = vcat([ltris[22].center; ladjs[22][2].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[22]], tri2idx[ladjs[22][2]]])
    #draw_lines!(plt, [ltris[22].edges[2]])
    #
    #pts = vcat([ltris[21].center; ladjs[21][1].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[21]], tri2idx[ladjs[21][1]]])
    #draw_lines!(plt, [ltris[21].edges[1]])
    #
    #pts = vcat([ltris[44].center; ladjs[44][1].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[44]], tri2idx[ladjs[44][1]]])
    #draw_lines!(plt, [ltris[44].edges[1]])
    #
    #
    #pts = vcat([ltris[56].center; ladjs[56][1].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[56]], tri2idx[ladjs[56][1]]])
    #draw_lines!(plt, [ltris[56].edges[1]])
    @test tri2idx[ladjs[56][1]] == 57
    #
    #
    #pts = vcat([ltris[556].center; ladjs[556][1].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[556]], tri2idx[ladjs[556][1]]])
    #draw_lines!(plt, [ltris[556].edges[1]])
    @test ismissing(ladjs[556][2])
    #
    #pts = vcat([ltris[600].center; ladjs[600][1].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[600]], tri2idx[ladjs[600][1]]])
    #draw_lines!(plt, [ltris[600].edges[1]])
    @test ismissing(ladjs[600][2])
    #
    #pts = vcat([ltris[580].center; ladjs[580][2].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[580]], tri2idx[ladjs[580][2]]])
    #draw_lines!(plt, [ltris[580].edges[2]])
    @test ismissing(ladjs[580][1])
    #
    #
    #pts = vcat([ltris[532].center; ladjs[532][2].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[532]], tri2idx[ladjs[532][2]]])
    #draw_lines!(plt, [ltris[532].edges[2]])
    @test ismissing(ladjs[532][1])
    #
    #pts = vcat([ltris[521].center; ladjs[521][2].center])
    #draw_points!(plt, pts; 
    #text=[tri2idx[ltris[521]], tri2idx[ladjs[521][2]]])
    #draw_lines!(plt, [ltris[521].edges[2]])
    @test tri2idx[ladjs[521][2]] == 522
    #
    #@savepng plt "sh2"
    #
end

