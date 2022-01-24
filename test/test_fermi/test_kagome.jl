"""
测试kagome的相关功能
"""

using MARY_fRG.Basics
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Patch
using MARY_fRG.Fermi.Surface
using MARY_fRG.Drawers
using Plots
using MARY_fRG.FlowEquation



"""
运行
"""
function run()
    nu1, nu2, nu3 = get_kagome_ν(0.1, 0.3)
    println(sum(nu1 .* nu1))
    println(sum(nu1 .* nu2))
    println(sum(nu1 .* nu3))
    println(nu1[1]*nu1[1] + nu2[1]*nu2[1] + nu3[1]*nu3[1])
    println(nu1[1]*nu1[2] + nu2[1]*nu2[2] + nu3[1]*nu3[2])
    #费米面的路径
    mpt = Point2D(0, 2*π / sqrt(3))
    mpt2 = Point2D(-π, π / sqrt(3))
    mpt3 = Point2D(-π, -π / sqrt(3))
    mpt4 = Point2D(0, -2*π / sqrt(3))
    mpt5 = Point2D(π, -π / sqrt(3))
    mpt6 = Point2D(π, π / sqrt(3))
    #mpt = Point2D(4π/3-0.8, 0)
    #mpt2 = Point2D(2π/3-0.4, 2π/√3-0.4*√3)
    #mpt3 = Point2D(-2π/3+0.4, 2π/√3-0.4*√3)
    #mpt4 = Point2D(-4π/3+0.8, 0)
    #mpt5 = Point2D(-2π/3+0.4, -2π/√3+0.4*√3)
    #mpt6 = Point2D(2π/3-0.4, -2π/√3+0.4*√3)
    xvals = [0.0]
    pps = [mpt]
    for idx in 0:1:18
        omega = (idx+1) / 19
        push!(xvals, omega)
        push!(pps, middle_point(mpt, mpt2; sc1=1-omega, sc2=omega))
    end
    #println(pps)
    #nu1, nu2, nu3 = get_kagome_ν(mpt2.x, mpt2.y)
    #println(nu1)
    #println(pps[11])
    #nu1, nu2, nu3 = get_kagome_ν(pps[11].x, pps[11].y)
    #println(nu2[1])
    #return
    push!(xvals, 1.0)
    push!(pps, mpt2)
    for idx in 0:1:18
        omega = (idx+1) / 19
        push!(xvals, 1+omega)
        push!(pps, middle_point(mpt2, mpt3; sc1=1-omega, sc2=omega))
    end
    push!(xvals, 2.0)
    push!(pps, mpt3)
    for idx in 0:1:18
        omega = (idx+1) / 19
        push!(xvals, 2+omega)
        push!(pps, middle_point(mpt3, mpt4, sc1=1-omega, sc2=omega))
    end
    push!(xvals, 3.0)
    push!(pps, mpt4)
    for idx in 0:1:18
        omega = (idx+1) / 19
        push!(xvals, 3+omega)
        push!(pps, middle_point(mpt4, mpt5, sc1=1-omega, sc2=omega))
    end
    push!(xvals, 4.0)
    push!(pps, mpt5)
    for idx in 0:1:18
        omega = (idx+1) / 19
        push!(xvals, 4+omega)
        push!(pps, middle_point(mpt5, mpt6, sc1=1-omega, sc2=omega))
    end
    push!(xvals, 5.0)
    push!(pps, mpt6)
    for idx in 0:1:18
        omega = (idx+1) / 19
        push!(xvals, 5+omega)
        push!(pps, middle_point(mpt6, mpt, sc1=1-omega, sc2=omega))
    end
    push!(xvals, 6.0)
    push!(pps, mpt)
    println(pps)
    sitea = []
    siteb = []
    sitec = []
    for kkt in pps
        nu1, nu2, nu3 = get_kagome_ν(kkt.x, kkt.y)
        push!(sitea, nu2[1])
        push!(siteb, nu2[2])
        push!(sitec, nu2[3])
        #push!(sitea, abs(nu2[1]))
        #push!(siteb, abs(nu2[2]))
        #push!(sitec, abs(nu2[3]))
    end
    println(siteb)
    plt = plot(xvals, sitea)
    plot!(plt, xvals, siteb)
    plot!(plt, xvals, sitec)
    png(plt, "kag_band_contrib1")
    #绘制一下d能带的贡献
    sitea = []
    siteb = []
    sitec = []
    #nu1, nu2, nu3 = get_nu_numpy(pps[10].coord[0], pps[10].coord[1])
    #print(nu2[1])
    #print(check_patches_converge(pps[10]))
    for kkt in pps
        nu1, nu2, nu3 = get_kagome_ν(kkt.x, kkt.y)
        #push!(sitea, nu3[1])
        #push!(siteb, nu3[2])
        #push!(sitec, nu3[3])
        push!(sitea, abs(nu3[1]))
        push!(siteb, abs(nu3[2]))
        push!(sitec, abs(nu3[3]))
    end
    plt = plot(xvals, sitea)
    plot!(plt, xvals, siteb)
    plot!(plt, xvals, sitec)
    png(plt, "kag_band_contrib2")
    #print(pps[10])
    #print(siteb[10])
    #nu1, nu2, nu3 = get_nu_numpy(pps[10].coord[0], pps[10].coord[1])
    #print(nu2[1])
    #print(check_patches_converge(pps[10]))
end


run()


function run2()
    mat = zeros(999, 999)
    for idx in CartesianIndices(mat)
        idx1, idx2 = Tuple(idx)
        kx = 4*π*(idx1 - 500) / 500 / 3
        ky = 4*π*(idx2 - 500) / 500 / 3
        x::BigFloat = BigFloat(kx) / BigFloat(4)
        y::BigFloat = sqrt(BigFloat(3)) * BigFloat(ky) / BigFloat(4)
        sqr = 2 * cos(2*x - 2*y)
        sqr += 2 * cos(2*x + 2*y)
        sqr += 2 * cos(4*x) + 3
        mat[idx1, idx2] = 1 - sqrt(sqr)
    end
    heatmap(mat)
    png("band2")
    sitea = []
    siteb = []
    sitec = []
    for aidx in 1:1:1000
        ang = 2*π*(aidx - 1) / 1000
        kx = 1.75*cos(ang)
        ky = 1.75*sin(ang)
        nu1, nu2, nu3 = get_kagome_ν(kx, ky)
        push!(sitea, abs(nu2[1]))
        push!(siteb, abs(nu2[2]))
        push!(sitec, abs(nu2[3]))
    end
    plt = plot(1:1:1000, sitea)
    plot!(plt, 1:1:1000, siteb)
    plot!(plt, 1:1:1000, sitec)
    png(plt, "kag_band_contrib3")
end

#run2()

"""
相互作用
"""
function interact()
    ks = upperband_kagome_lattice(0.)
    pats = patches_under_vonhove(ks.brillouin, ks.dispersion[1], 24)
    #pats = [
    #    Point2D(0.1, -0.1),
    #    Point2D(1., 1.),
    #    Point2D(-1., -1.)
    #]
    plt = draw_lines(ks.brillouin.edges)
    draw_points!(plt, pats)
    @savepng plt "kag_surf1"
    #pats = get_von_hove_patches(24)
    #vmat = band_v(6.0, pats)
    #draw_heatmap(vmat[:, :, 0])
    #println(vmat[1, 2, 1], vmat[2, 1, 2])
    #println(vmat[2, 0, 0], vmat[0, 2, 2])
    #println(vmat[2, 0, 1], vmat[0, 2, 1])
    println(pats)
    umat = get_kagome_U_mt(6.0, pats)
    println(all(isfinite.(umat)))
    maxv, maxi = findmin(isfinite.(umat))
    println(maxi, maxv)
    println(umat[maxi])
    #
    #k1v = pats[maxi[1]]
    #k2v = pats[maxi[2]]
    #k3v = pats[maxi[3]]
    #k4v = Fermi.hexagon_kadd(k1v, k2v)
    #k4v = Fermi.hexagon_kadd(k4v, -k3v)
    #_, k1nu, _ = get_kagome_ν(k1v.x, k1v.y)
    #_, k2nu, _ = get_kagome_ν(k2v.x, k2v.y)
    #_, k3nu, _ = get_kagome_ν(k3v.x, k3v.y)
    #_, k4nu, _ = get_kagome_ν(k4v.x, k4v.y)
    #println(k1v)
    #println(k2v)
    #println(k3v)
    #println(k4v)
    #println(k1nu)
    #println(k2nu)
    #println(k3nu)
    #println(k4nu)
    #
    plt = heatmap(umat[:, :, 1])
    png(plt, "kag_interact")
    #amat = umat + vmat
    #draw_heatmap(amat[:, :, 0])
    #smat = get_sublattice_u(umat, pats)
    #print("site")
    #draw_heatmap(smat[0, 0, 0, 0, :, :, 0])
    #draw_heatmap(smat[1, 1, 1, 1, :, :, 2])
    #draw_heatmap(smat[0, 1, 1, 0, :, :, 0])
    #draw_heatmap(smat[2, 2, 2, 2, :, :, 5])
end

#interact()


#get_kagome_ν(0., 0.)

"""
计算ult
"""
function get_ult_bubb(Γ4; maxnum=1000000)
    Γ4_ult = fliter_away_surface(Γ4, 2.0)
    Γ4_ult = refine_to_surface(Γ4_ult)
    while length(Γ4_ult.ltris) < maxnum
        Γ4_ult = refine_to_surface(Γ4_ult)
    end
    #
    maxi, mini = engpeak_to_surface(Γ4_ult)
    bubb_pp_u, bubb_fs_u, bubb_ex_u = all_bubble_tf_ult_mt(Γ4_ult, maxi)
    return bubb_pp_u, bubb_fs_u, bubb_ex_u
end

function kprofile()
    model = upperband_kagome_lattice(0.)
    Γ4 = TFGamma4(
        model, 6.0, 24, 200
    )
    @time bpu, bfu, beu = get_ult_bubb(Γ4)
    @time bpu, bfu, beu = get_ult_bubb(Γ4)
end

#kprofile()
