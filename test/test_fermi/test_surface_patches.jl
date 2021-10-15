"""
测试两个功能
"""

using Test
using MARY_fRG.Drawers
using MARY_fRG.Triangulated
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Surface
using MARY_fRG.Fermi.Patch


#@testset "正方晶系" begin
#    quad = common_square_lattice(0.)
#    ltris, ladjs = split_square(quad.brillouin, 100)
#    pnum = 44
#    lpats = [find_patch_index_squa(tri.center, quad.brillouin, pnum) for tri in ltris]
#    sf, sfidx = const_energy_line_in_patches(ltris, ladjs, lpats, 0., quad.dispersion[1])
#    plt = draw_lines(quad.brillouin.edges)
#    draw_lines!(plt, sf; lcidx=sfidx)
#    sf2, sf2idx = const_energy_line_in_patches(ltris, ladjs, lpats, -0.2, quad.dispersion[1])
#    draw_lines!(plt, sf2; lcidx=sf2idx)
#    @savepng plt "surface1"
#end
#
#
#@testset "三角格子的费米面" begin
#    ts = common_triangle_lattice(0.)
#    ltris, ladjs = split_hexagon(ts.brillouin, 100)
#    pnum = 42
#    lpats = [find_patch_index_hexa(tri.center, ts.brillouin, pnum) for tri in ltris]
#    sf, sfidx = const_energy_line_in_patches(ltris, ladjs, lpats, 0., ts.dispersion[1])
#    plt = draw_lines(ts.brillouin.edges)
#    draw_lines!(plt, sf; lcidx=sfidx)
#    sf2, sf2idx = const_energy_line_in_patches(ltris, ladjs, lpats, 2.0, ts.dispersion[1])
#    draw_lines!(plt, sf2; lcidx=sf2idx)
#    @savepng plt "surface2"
#end



@testset "kagome格子的费米面" begin
    ts = upperband_kagome_lattice(-0.2)
    ltris, ladjs = split_hexagon(ts.brillouin, 100)
    pnum = 42
    lpats = [find_patch_index_hexa(tri.center, ts.brillouin, pnum) for tri in ltris]
    sf, sfidx = const_energy_line_in_patches(ltris, ladjs, lpats, 0., ts.dispersion[1])
    plt = draw_lines(ts.brillouin.edges)
    draw_lines!(plt, sf; lcidx=sfidx)
    sf2, sf2idx = const_energy_line_in_patches(ltris, ladjs, lpats, 2.0, ts.dispersion[1])
    draw_lines!(plt, sf2; lcidx=sf2idx)
    @savepng plt "surface3"
end

