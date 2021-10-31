"""
测试单带的正方格子
"""

using Plots
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation


#阻止画图
ENV["GKSwstype"] = "100"


"""
从头开始运行代码
"""
function run_ec()
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


#run_ec()


function run_tf()
    model = common_square_lattice(0.00)
    Γ4 = TFGamma4(
        model, 8.0, 16, 100
    )
    Γ4.V .+= 1.0
    lval = 0.
    lstep = 0.01
    for idx in 1:1:701
        bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, lval)
        dl = dl_tf_mt(Γ4, bubb_pp, bubb_fs, bubb_ex)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
            png(plt, "Gamma4"*string(idx))
        end
        lval += lstep
    end
end


run_tf()
