using Test
using MARY_fRG.Basics
using MARY_fRG.Triangulated
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Patch


#@testset "三角晶系" begin
#    hexa1 = common_triangle_lattice(0.)
#    #
#    ltris, ladjs = split_hexagon(hexa1.brillouin, 200)
#    pnum = 42
#    lpats = group_ltris_into_patches_mt(ltris, hexa1.brillouin, pnum)
#    patcs = patches_under_vonhove(hexa1.brillouin, hexa1.dispersion[1], pnum)
#    #lpats = [find_patch_index_hexa(tri.center, hexa1.brillouin, pnum) for tri in ltris]
#    #plt = draw_polygon(ltris; pcidx=lpats)
#    #@savepng plt "patch1"
#    for idx in 1:1:pnum
#        plt = draw_lines(hexa1.brillouin.edges)
#        pltris::Vector{Basics.AbstractTriangle} = []
#        for (tri, pid) in zip(ltris, lpats)
#            if pid == idx
#                push!(pltris, tri)
#            end
#        end
#        draw_points!(plt, [patcs[idx]]; text=idx)
#        draw_polygon!(plt, pltris; pcidx=ones(Int64, length(pltris)))
#        @savepng plt "patch"*string(idx)
#    end
#end


@testset "正方晶系" begin
    quad = common_square_lattice(1.0)
    #
    ltris, ladjs = split_square(quad.brillouin, 200)
    #find_patch_index_squa(ltris[613].center, quad.brillouin, 28)
    pnum = 44
    lpats = group_ltris_into_patches_mt(ltris, quad.brillouin, pnum)
    patcs = patches_under_vonhove(quad.brillouin, quad.dispersion[1], pnum)
    #[find_patch_index_squa(tri.center, quad.brillouin, pnum) for tri in ltris]
    for idx in 1:1:pnum
        plt = draw_lines(quad.brillouin.edges)
        pltris::Vector{Basics.AbstractTriangle} = []
        for (tri, pid) in zip(ltris, lpats)
            if pid == idx
                push!(pltris, tri)
            end
        end
        draw_points!(plt, [patcs[idx]]; text=idx)
        draw_polygon!(plt, pltris; pcidx=ones(Int64, length(pltris)))
        @savepng plt "patch"*string(idx)
    end
end

