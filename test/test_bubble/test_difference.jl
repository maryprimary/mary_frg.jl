"""测试不同的粒度，bubble的区别"""


using Plots
using MARY_fRG.Basics
using MARY_fRG.Drawers
using MARY_fRG.Fermi
using MARY_fRG.Fermi.Patch
using MARY_fRG.FlowEquation


"""筛选tris"""
function fliter(Γ4::Gamma4{T, P}, lval) where {T, P}
    lamb = Γ4.λ_0 * exp(-lval)
    resltris::Vector{P} = []
    reslpats::Vector{Int64} = []
    cert = 1#length(Γ4.ltris) > 2000000 ? 5 : 15
    mintri = missing
    minpat = missing
    mineng = 10.
    for (tri, pat) in zip(Γ4.ltris, Γ4.lpats)
        if pat != 1
            continue
        end
        engs = [Γ4.model.dispersion[bidx](
            tri.center.x, tri.center.y) / lamb for bidx in 1:1:Γ4.model.bandnum
        ]
        #如果还有可能有贡献
        if minimum(abs.(engs)) < cert #&& minimum(abs.(engs2)) > cert # && pat == 1# && minimum(abs.(engs)) > 0.1
            push!(resltris, tri)
            push!(reslpats, pat)
        end
        #if minimum(abs.(engs)) < mineng
        #    mintri = tri
        #    minpat = pat
        #    mineng = minimum(abs.(engs))
        #end
    end
    #
    #resltris = [mintri]#resltris[1:1]
    #reslpats = [minpat]#reslpats[1:1]
    #
    ltris_pat = Vector{typeof(resltris)}(undef, Γ4.patchnum)
    for idx in 1:1:Γ4.patchnum
        ltris_pat[idx] = []
    end
    #
    for (tri, pat) in zip(resltris, reslpats)
        push!(ltris_pat[pat], tri)
    end
    #
    return Gamma4(
        Γ4.model,
        Γ4.λ_0,
        Γ4.V,
        Γ4.k4tab,
        Γ4.patchnum,
        Γ4.patches,
        resltris, reslpats, nothing, ltris_pat
    )
end


function refine(Γ4::Gamma4{T, P}) where {T, P}
    resltris::Vector{P} = []
    reslpats::Vector{Int64} = []
    resltris = refine_list_triangle(Γ4.ltris, 5)
    reslpats = group_ltris_into_patches_mt(resltris, Γ4.model.brillouin, Γ4.patchnum)
    #
    #resltris = [mintri]#resltris[1:1]
    #reslpats = [minpat]#reslpats[1:1]
    #
    ltris_pat = Vector{typeof(resltris)}(undef, Γ4.patchnum)
    for idx in 1:1:Γ4.patchnum
        ltris_pat[idx] = []
    end
    #
    for (tri, pat) in zip(resltris, reslpats)
        push!(ltris_pat[pat], tri)
    end
    #
    return Gamma4(
        Γ4.model,
        Γ4.λ_0,
        Γ4.V,
        Γ4.k4tab,
        Γ4.patchnum,
        Γ4.patches,
        resltris, reslpats, nothing, ltris_pat
    )
end


function run_rftf()
    model = common_square_lattice(0.20)
    disp = model.dispersion[1]
    lval = 15.
    #500
    Γ4_5 = TFGamma4(
        model, 8.0, 16, 500
    )
    Γ4_5.V .+= 1.0
    Γ4_5_1 = fliter(Γ4_5, 11.0)
    Γ4_5_2 = fliter(Γ4_5, 12.0)
    println(length(Γ4_5_1.ltris))
    println(length(Γ4_5_2.ltris))
    #fliter2
    bubbval_fs = Array{Float64, 2}(undef, Γ4_5_2.patchnum, Γ4_5_2.patchnum)
    lamb = Γ4_5_2.λ_0 * exp(-lval)
    println("lamb ", lamb)
    Threads.@threads for idxs in CartesianIndices(bubbval_fs)
        i2, i3 = Tuple(idxs)
        k2, k3 = Γ4_5_2.patches[1][i2], Γ4_5_2.patches[1][i3]
        q_fs = Γ4_5_2.model.kadd(k3, -k2)
        bubbres = pi_αβ_plus_tf(
            Γ4_5_2.ltris, area(Γ4_5_2.model.brillouin),
            lamb, q_fs, Γ4_5_2.model.dispersion[1],
            Γ4_5_2.model.dispersion[1], Γ4_5_2.model.kadd
        )
        bubbval_fs[i2, i3] = bubbres
    end
    plt = heatmap(1e5*bubbval_fs)
    @savepng plt "bubb_5_2_2"
    #fliter1
    bubbval_fs = Array{Float64, 2}(undef, Γ4_5_1.patchnum, Γ4_5_1.patchnum)
    lamb = Γ4_5_1.λ_0 * exp(-lval)
    Threads.@threads for idxs in CartesianIndices(bubbval_fs)
        i2, i3 = Tuple(idxs)
        k2, k3 = Γ4_5_1.patches[1][i2], Γ4_5_1.patches[1][i3]
        q_fs = Γ4_5_1.model.kadd(k3, -k2)
        bubbres = pi_αβ_plus_tf(
            Γ4_5_1.ltris, area(Γ4_5_1.model.brillouin),
            lamb, q_fs, Γ4_5_1.model.dispersion[1],
            Γ4_5_1.model.dispersion[1], Γ4_5_1.model.kadd
        )
        bubbval_fs[i2, i3] = bubbres
    end
    plt = heatmap(1e5*bubbval_fs)
    @savepng plt "bubb_5_2_1"
    #fliter3
    ltris3::Vector{Basics.AbstractTriangle{:RT}} = []
    ltris4::Vector{Basics.AbstractTriangle{:RT}} = []
    for tri in Γ4_5_1.ltris
        if !any(map((x)->x==tri, Γ4_5_2.ltris))
            push!(ltris4, tri)
        end
    end
    #push!(ltris3, Γ4_5_2.ltris[1])
    #println(Γ4_5_2.ltris[1])
    #println(disp(Γ4_5_2.ltris[1].center.x, Γ4_5_2.ltris[1].center.y))
    #push!(ltris3, Γ4_5_2.ltris[2])
    push!(ltris3, Γ4_5_2.ltris[3])
    #push!(ltris3, ltris4[1])
    push!(ltris3, ltris4[2])
    println(disp(ltris4[2].center.x, ltris4[2].center.y))
    #push!(ltris3, ltris4[3])
    println(length(ltris3))
    bubbval_fs = Array{Float64, 2}(undef, Γ4_5_1.patchnum, Γ4_5_1.patchnum)
    lamb = Γ4_5_1.λ_0 * exp(-lval)
    Threads.@threads for idxs in CartesianIndices(bubbval_fs)
        i2, i3 = Tuple(idxs)
        k2, k3 = Γ4_5_1.patches[1][i2], Γ4_5_1.patches[1][i3]
        q_fs = Γ4_5_1.model.kadd(k3, -k2)
        bubbres = pi_αβ_plus_tf(
            ltris3, area(Γ4_5_1.model.brillouin),
            lamb, q_fs, Γ4_5_1.model.dispersion[1],
            Γ4_5_1.model.dispersion[1], Γ4_5_1.model.kadd
        )
        bubbval_fs[i2, i3] = bubbres
    end
    plt = heatmap(1e5*bubbval_fs)
    @savepng plt "bubb_5_2_3"
    return
    #plt = draw_points([tri.center for tri in Γ4_5.ltris])
    #abses = [abs(disp(tri.center.x, tri.center.y)) for tri in Γ4_5.ltris]
    #println(minimum(abses))
    #abscount = 0
    #for ae in abses
    #    if ae < 0.01
    #        abscount += 1
    #    end
    #end
    #println(abscount)
    #@savepng plt "diffe1"
    println(Γ4_5.ltris[1])
    bubb_pp_5, bubb_fs_5, bubb_ex_5 = all_bubble_tf_mt(Γ4_5, lval)
    println("500时的贡献 ", sum(abs.(bubb_fs_5.V)))
    println("500时最大", findmax(abs.(bubb_fs_5.V)))
    k1, k2 = Γ4_5.patches[1][12], Γ4_5.patches[1][4]
    println(k1)
    println(k2)
    kprim = Γ4_5.model.kadd(k1, k2)
    println(kprim)
    kprim = Γ4_5.model.kadd(kprim, -Γ4_5.ltris[1].center)
    println(kprim)
    println("kprim disp ", Γ4_5.model.dispersion[1](kprim.x, kprim.y))
    
    return
    ###
    #100
    Γ4_1 = TFGamma4(
        model, 8.0, 16, 800
    )
    Γ4_1.V .+= 1.0
    #Γ4_1 = fliter_away_surface(Γ4_1, 10.0)
    #plt = draw_points([tri.center for tri in Γ4_1.ltris])
    #@savepng plt "diffe3"
    #Γ4_1 = TFGamma4_addition_ltris_mt(Γ4_1, 10.0)
    Γ4_1 = fliter(Γ4_1, 12.0)
    #println(Γ4_1.ltris)
    #plt = draw_polygon(Γ4_1.ltris)
    #plot!(plt, aspect_ratio=:equal)
    #@savepng plt "diffe4"
    Γ4_1 = refine(Γ4_1)
    #println(Γ4_1.ltris)
    plt = draw_polygon(Γ4_1.ltris)
    #plot!(plt, aspect_ratio=:equal)
    #@savepng plt "diffe5"
    Γ4_1 = fliter(Γ4_1, 12.0)
    abses = [abs(disp(tri.center.x, tri.center.y)) for tri in Γ4_1.ltris]
    println(minimum(abses))
    abscount = 0
    for ae in abses
        if ae < 0.01
            abscount += 1
        end
    end
    println(abscount)
    #println(minimum([abs(disp(tri.center.x, tri.center.y)) for tri in Γ4_1.ltris]))
    plt = draw_points([tri.center for tri in Γ4_1.ltris])
    @savepng plt "diffe2"
    println(Γ4_1.ltris[1])
    k1, k2 = Γ4_1.patches[1][11], Γ4_1.patches[1][3]
    println(k1)
    println(k2)
    kprim = Γ4_1.model.kadd(k1, k2)
    println(kprim)
    kprim = Γ4_1.model.kadd(kprim, -Γ4_1.ltris[1].center)
    println(kprim)
    println("kprim disp ", Γ4_1.model.dispersion[1](kprim.x, kprim.y))
    bubb_pp_1, bubb_fs_1, bubb_ex_1 = all_bubble_tf_mt(Γ4_1, lval)
    println("200时的贡献 ", sum(abs.(bubb_pp_1.V)))
    println("200时最大", findmax(abs.(bubb_pp_1.V)))
    #500和100的差
    diffe = bubb_pp_5.V - bubb_pp_1.V
    diffe = abs.(diffe)
    println(sum(diffe))
    println(sum(abs.(bubb_pp_1.V)))
    #refine
    #Γ4_r = TFGamma4_addition_ltris_mt(Γ4_1, lval)
    #bubb_pp_r, bubb_fs_r, bubb_ex_r = all_bubble_tf_mt(Γ4_r, lval)
    ##500和100的差
    #diffe = bubb_pp_5.V - bubb_pp_r.V
    #diffe = abs.(diffe)
    #println(sum(diffe))
end


#run_rftf()



function run_rf_ult()
    model = common_square_lattice(0.20)
    disp = model.dispersion[1]
    lval = 17.
    #500
    Γ4_5 = TFGamma4(
        model, 8.0, 16, 500
    )
    Γ4_5.V .+= 1.0
    #Γ4_5 = fliter(Γ4_5, 12.0)
    #plt = draw_points([tri.center for tri in Γ4_5.ltris])
    #abses = [abs(disp(tri.center.x, tri.center.y)) for tri in Γ4_5.ltris]
    #println(minimum(abses))
    #abscount = 0
    #for ae in abses
    #    if ae < 0.01
    #        abscount += 1
    #    end
    #end
    #println(abscount)
    #@savepng plt "diffe1"
    println(Γ4_5.ltris[1])
    bubb_pp_5, bubb_fs_5, bubb_ex_5 = all_bubble_tf_mt(Γ4_5, lval)
    println("500时的贡献 ", sum(abs.(bubb_pp_5.V)))
    println("500时最大", findmax(abs.(bubb_pp_5.V)))
    ###
    #100
    Γ4_1 = TFGamma4(
        model, 8.0, 16, 800
    )
    Γ4_1.V .+= 1.0
    bubb_pp_1, bubb_fs_1, bubb_ex_1 = all_bubble_tf_mt(Γ4_1, lval)
    println("200时的贡献 ", sum(abs.(bubb_pp_1.V)))
    println("200时最大", findmax(abs.(bubb_pp_1.V)))
    ###
    #ult
    Γ4_1 = fliter_away_surface(Γ4_1, 6.0)
    Γ4_1 = refine_to_surface(Γ4_1)
    while length(Γ4_1.ltris) < 1000000
        Γ4_1 = refine_to_surface(Γ4_1)
    end
    plt = draw_points([tri.center for tri in Γ4_1.ltris])
    @savepng plt "diffe6"
    println("tri大小", length(Γ4_1.ltris))
    maxi, mini = engpeak_to_surface(Γ4_1)
    println(maxi, " ", mini)
    bubb_pp_u, bubb_fs_u, bubb_ex_u = all_bubble_tf_ult_mt(Γ4_1, maxi)
    println("ult时的贡献 ", sum(abs.(bubb_pp_u.V)))
    println("ult时最大", findmax(abs.(bubb_pp_u.V)))
    #500和100的差
    diffe = bubb_pp_5.V - bubb_pp_1.V
    diffe = abs.(diffe)
    println(sum(diffe))
    println(sum(abs.(bubb_pp_1.V)))
    #500个ult的差
    diffe = bubb_pp_5.V - bubb_pp_u.V
    diffe = abs.(diffe)
    println(sum(diffe))
    println(sum(abs.(bubb_pp_u.V)))
    #ult和500的画图
    plt = heatmap(1e5*bubb_pp_u.V[1, 1, 1, 1, 1, :, :])
    png(plt, "bubb_u1")
    plt = heatmap(1e5*bubb_pp_u.V[1, 1, 1, 1, 2, :, :])
    png(plt, "bubb_u2")
    plt = heatmap(1e5*bubb_pp_u.V[1, 1, 1, 1, 3, :, :])
    png(plt, "bubb_u3")
    plt = heatmap(1e5*bubb_pp_u.V[1, 1, 1, 1, 4, :, :])
    png(plt, "bubb_u4")
    plt = heatmap(1e5*bubb_pp_5.V[1, 1, 1, 1, 1, :, :])
    png(plt, "bubb_51")
    plt = heatmap(1e5*bubb_pp_5.V[1, 1, 1, 1, 2, :, :])
    png(plt, "bubb_52")
    plt = heatmap(1e5*bubb_pp_5.V[1, 1, 1, 1, 3, :, :])
    png(plt, "bubb_53")
    plt = heatmap(1e5*bubb_pp_5.V[1, 1, 1, 1, 4, :, :])
    png(plt, "bubb_54")
    println(Γ4_1.patches[1][1], Γ4_1.patches[1][4], Γ4_1.patches[1][9])
end


run_rf_ult()
