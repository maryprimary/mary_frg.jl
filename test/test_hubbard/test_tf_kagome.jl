"""
测试上VHS的kagome
"""

using Plots
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation

#阻止画图
ENV["GKSwstype"] = "100"

function run_tf()
    model = upperband_kagome_lattice(0.)
    Γ4 = TFGamma4CMPLX(
        model, 6.0, 24, 223
    )
    println(length(Γ4.ltris))
    Γ4.V[1, 1, 1, 1, :, :, :] .+= get_kagome_U_mt(Γ4, 6.0)
    println(
        all(isfinite.(Γ4.V[1, 1, 1, 1, :, :, :]))
    )
    #plt = heatmap(Γ4.V[1, 1, 1, 1, :, :, 1])
    #png(plt, "Gamma41")
    #return
    lval = 0.
    lstep = 0.01
    for idx in 1:1:3501
        bubb_pp, bubb_fs, bubb_ex = all_bubble_tf_mt(Γ4, lval)
        dl = dl_tf_mt(Γ4, bubb_pp, bubb_fs, bubb_ex)
        Γ4.V .+= dl .* lstep
        if idx % 50 == 1
            plt = heatmap(real(Γ4.V[1, 1, 1, 1, :, :, 1]))
            png(plt, "Gamma4"*string(idx))
            plt = heatmap(imag(Γ4.V[1, 1, 1, 1, :, :, 1]))
            png(plt, "iGamma4"*string(idx))
        end
        lval += lstep
    end
end


run_tf()
