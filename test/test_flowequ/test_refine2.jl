"""
测试refine Surface
"""

using Test
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation
using MARY_fRG.Drawers
using MARY_fRG.Fermi.Surface


#@testset "靠近费米面" begin
#    qs = common_square_lattice(0.20)
#    Γ4 = TFGamma4(qs, 8.0, 16, 20)
#    surf = const_energy_triangle(Γ4.ltris, 0.5, Γ4.model.dispersion[1])
#    plt = draw_polygon(surf; pcidx=ones(Int64, length(surf)))
#    @savepng plt "surf_refine2"
#    surf2 = refine_to_surface(surf, 0.5, Γ4.model.dispersion[1])
#    plt = draw_polygon(surf2; pcidx=ones(Int64, length(surf2)))
#    @savepng plt "surf_refine22"
#    surf3 = refine_to_surface(surf2, 0.5, Γ4.model.dispersion[1])
#    plt = draw_polygon(surf3; pcidx=ones(Int64, length(surf3)))
#    @savepng plt "surf_refine23"
#end
#
#
#
#@testset "细化费米面" begin
#    qs = common_square_lattice(0.20)
#    Γ4 = TFGamma4(qs, 8.0, 16, 20)
#    surf = const_energy_triangle(Γ4.ltris, 0.5, Γ4.model.dispersion[1])
#    plt = draw_points([tri.center for tri in surf])
#    @savepng plt "surf_refine24"
#    surf2 = refine_list_triangle(surf, 1)
#    plt = draw_points([tri.center for tri in surf2])
#    @savepng plt "surf_refine25"
#    surf3 = refine_list_triangle(surf2, 3)
#    plt = draw_points([tri.center for tri in surf3])
#    @savepng plt "surf_refine26"
#end
#
#
#@testset "addition" begin
#    qs = common_square_lattice(0.20)
#    Γ4 = TFGamma4(qs, 8.0, 16, 100)
#    posi, nega, away = TFGamma4_addition_ltris_mt(Γ4, 2.0)
#    plt = draw_points([tri.center for tri in away])
#    draw_points!(plt, [tri.center for tri in posi])
#    draw_points!(plt, [tri.center for tri in nega])
#    @savepng plt "surf_refine27"
#    posi, nega, away = TFGamma4_addition_ltris_mt(Γ4, 4.0)
#    plt = draw_points([tri.center for tri in away])
#    draw_points!(plt, [tri.center for tri in posi])
#    draw_points!(plt, [tri.center for tri in nega])
#    @savepng plt "surf_refine28"
#    posi, nega, away = TFGamma4_addition_ltris_mt(Γ4, 8.0)
#    plt = draw_points([tri.center for tri in away])
#    draw_points!(plt, [tri.center for tri in posi])
#    draw_points!(plt, [tri.center for tri in nega])
#    @savepng plt "surf_refine29"
#    pt1 = posi[1].center
#    pt2 = posi[2].center
#    println(Γ4.model.dispersion[1](pt1.x, pt1.y))
#    println(Γ4.model.dispersion[1](pt2.x, pt2.y))
#end


qs = common_square_lattice(0.20)
Γ4 = TFGamma4(qs, 8.0, 16, 100)
nΓ4 = TFGamma4_addition_ltris_mt(Γ4, 4.0)
plt = draw_points([tri.center for tri in nΓ4.ltris])
@savepng plt "surf_refine210"
