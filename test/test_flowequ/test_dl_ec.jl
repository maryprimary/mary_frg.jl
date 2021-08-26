"""
测试能量截断的求导
"""

using Test
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation


@testset "能量截断" begin
    model = common_square_lattice(0.)
    Γ4 = ECGamma4(model, 4.0, 8, 40)
    Γ4.V .+= 1.0
    println(Γ4.patches)
    bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex =
    all_bubble_ec_mt(Γ4, 1.0)
    dl = dl_ec_mt(Γ4, bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 2], -0.1289112, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 3], -0.0304875, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 4], 0.0032081, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 5], 0.2625107, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 2], 0.0336956, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 3], 0.0304875, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 4], -0.098423, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 5], 0., atol=1e-6)
end

