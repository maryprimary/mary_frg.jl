"""
测试ult时的计算
"""

using Plots
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation
using Profile


"""
获取极低温的时候的bubble
"""
function get_ult_bubb(Γ4; maxnum=10000)
    Γ4_ult = fliter_away_surface(Γ4, 2.0)
    Γ4_ult = refine_to_surface(Γ4_ult)
    while length(Γ4_ult.ltris) < maxnum
        Γ4_ult = refine_to_surface(Γ4_ult)
    end
    #plt = draw_points([tri.center for tri in Γ4_ult.ltris])
    #@savepng plt "ult_surface"
    maxi, mini = engpeak_to_surface(Γ4_ult)
    bubb_pp_u, bubb_fs_u, bubb_ex_u = @profile all_bubble_tf_ult_mt(Γ4_ult, maxi)
    return bubb_pp_u, bubb_fs_u, bubb_ex_u
end


"""
运行
"""
function run()
    model = upperband_kagome_lattice(0.)
    Γ4 = TFGamma4(
        model, 6.0, 24, 200
    )
    bpu, bfu, beu = get_ult_bubb(Γ4)
    #
    println(all(isfinite.(bpu.V)))
    println(all(isfinite.(bfu.V)))
    println(all(isfinite.(beu.V)))
    #
    umat = bpu.V[1, 1, 1, 1, :, :, :]
    maxv, maxi = findmin(isfinite.(umat))
    println(maxi, maxv)
    println(umat[maxi])
end


run()

run()


outf = open("kag_profile", "w")
Profile.print(outf, format = :flat)
close(outf)

