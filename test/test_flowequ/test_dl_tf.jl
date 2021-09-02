"""
测试温度流求导
"""

using Test
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation


@testset "温度流" begin
    model = common_square_lattice(0.)
    Γ4 = TFGamma4(model, 4.0, 8, 40)
    Γ4.V .+= 1.0
    bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, 1.0)
    @time dl = dl_tf_mt(Γ4, bubb_pp, bubb_fs, bubb_ex)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 2], -0.0127333, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 3], -0.0038054, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 4], 0.0143136, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 2, 2, 5], 0.0441989, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 2], 0.0181191, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 3], 0.0038054, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 4], -0.0089279, atol=1e-6)
    @test isapprox(dl[1, 1, 1, 1, 4, 2, 5], 0., atol=1e-6)
end

