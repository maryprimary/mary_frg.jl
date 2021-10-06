"""测试不同的粒度，bubble的区别"""


using Plots
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
    lval = 10.
    #500
    Γ4_5 = TFGamma4(
        model, 8.0, 16, 500
    )
    Γ4_5.V .+= 1.0
    Γ4_5 = fliter(Γ4_5, 12.0)
    plt = draw_points([tri.center for tri in Γ4_5.ltris])
    abses = [abs(disp(tri.center.x, tri.center.y)) for tri in Γ4_5.ltris]
    println(minimum(abses))
    abscount = 0
    for ae in abses
        if ae < 0.01
            abscount += 1
        end
    end
    println(abscount)
    @savepng plt "diffe1"
    println(Γ4_5.ltris[1])
    bubb_pp_5, bubb_fs_5, bubb_ex_5 = all_bubble_tf_mt(Γ4_5, lval)
    println("500时的贡献 ", sum(abs.(bubb_pp_5.V)))
    println("500时最大", findmax(abs.(bubb_pp_5.V)))
    k1, k2 = Γ4_5.patches[1][12], Γ4_5.patches[1][4]
    println(k1)
    println(k2)
    kprim = Γ4_5.model.kadd(k1, k2)
    println(kprim)
    kprim = Γ4_5.model.kadd(kprim, -Γ4_5.ltris[1].center)
    println(kprim)
    println("kprim disp ", Γ4_5.model.dispersion[1](kprim.x, kprim.y))
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


run_rftf()
