"""
测试单带的正方格子
"""

using Plots
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation

"""
从头开始运行代码
"""
function run()
    model = common_square_lattice(0.20)
    Γ4 = ECGamma4(
        model, 4.0, 16, 50
    )
    Γ4.V .+= 1.0
    lval = 0.
    lstep = 0.01
    for idx in 1:1:701
        bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex = all_bubble_ec_mt(Γ4, lval)
        dl = dl_ec_mt(Γ4, bubb_pp, bubb_fs, bubb_nfs, bubb_ex, bubb_nex)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma4"*string(idx))
        end
        lval += lstep
    end
end


run()
