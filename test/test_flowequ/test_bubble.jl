using Test
using MARY_fRG.Basics
using MARY_fRG.Fermi, MARY_fRG.Triangulated
using MARY_fRG.Fermi.Surface, MARY_fRG.Fermi.Patch
using MARY_fRG.FlowEquation


@testset "ec测试结果" begin
    qs = common_square_lattice(0.)
    pnum = 44
    ltris, ladjs = split_square(qs.brillouin, 40)
    lpats = group_ltris_into_patches_mt(ltris, qs.brillouin, pnum)
    sf = const_energy_line(ltris, ladjs, 0., qs.dispersion[1])
    res = pi_αβ_plus_ec(
        sf, sf, 4pi^2, 0.2, Point2D(0., 0.2), qs.dispersion[1], qs.kadd
    )
    @test isapprox(res, 0.116678333187, atol=1e-12)
    res = pi_αβ_minus_ec(
        sf, sf, 4pi^2, 0.2, Point2D(-0.15, 0.), qs.dispersion[1], qs.kadd
    )
    @test isapprox(res, 0.107264645635, atol=1e-12)
    #
    pnum = 8
    lpats = group_ltris_into_patches_mt(ltris, qs.brillouin, pnum)
    patcs = patches_under_vonhove(qs.brillouin, qs.dispersion[1], pnum)
    mpatcs = [patcs]
    #Γ4 = Gamma4(
    #    qs, 8.0, zeros(1, 1, 1, 1, pnum, pnum, pnum),
    #    zeros(Int64, 1, 1, 1, 1, pnum, pnum, pnum),
    #    pnum, mpatcs,
    #    ltris, lpats, ladjs, nothing
    #)
    Γ4 = ECGamma4(qs, 8.0, pnum, 40)
    @time all_bubble_ec_mt(Γ4, 1.0)
end


@testset "tf测试结果" begin
    qs = common_square_lattice(0.)
    pnum = 44
    ltris, ladjs = split_square(qs.brillouin, 40)
    lpats = group_ltris_into_patches_mt(ltris, qs.brillouin, pnum)
    res = pi_αβ_plus_tf(
        ltris, 4pi^2, area(ltris[1]), 0.2, Point2D(0., 0.2),
        qs.dispersion[1], qs.dispersion[1], qs.kadd
    )
    @test isapprox(res, 0.0503171, atol=1e-6)
    res = pi_αβ_minus_tf(
        ltris, 4pi^2, area(ltris[1]), 0.2, Point2D(-0.15, 0.),
        qs.dispersion[1], qs.dispersion[1], qs.kadd
    )
    @test isapprox(res, 0.2063251, atol=1e-6)
    #
    pnum = 8
    lpats = group_ltris_into_patches_mt(ltris, qs.brillouin, pnum)
    patcs = patches_under_vonhove(qs.brillouin, qs.dispersion[1], pnum)
    mpatcs = [patcs]
    #Γ4 = Gamma4(
    #    qs, 8.0, zeros(1, 1, 1, 1, pnum, pnum, pnum),
    #    zeros(Int64, 1, 1, 1, 1, pnum, pnum, pnum),
    #    pnum, mpatcs,
    #    ltris, lpats, ladjs, nothing
    #)
    Γ4 = TFGamma4(qs, 8.0, pnum, 40)
    @time all_bubble_tf_mt(Γ4, 1.0)
    #println(Γ4.ltris_pat)
    @time all_bubble_tf_mt(Γ4, 1.0)
    #println(Γ4.ltris_pat)
end
