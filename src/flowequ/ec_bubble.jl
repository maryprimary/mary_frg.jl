"""
温度流的Bubble积分
"""

"""
获取所有的Bubble(pp, fs, nfs, ex, nex)
qpp = k1+k2; qfs = k3-k2; qex = k1-k3
"""
function all_bubble_ec_mt(Γ4::Gamma4, lval)
    #
    lamb = Γ4.λ_0 * exp(-lval)
    brlu_area = area(Γ4.model.brillouin)
    #首先得到需要被积分的等能线
    posimat = Matrix{Vector{Segment}}(undef, Γ4.model.bandnum, Γ4.patchnum)
    negamat = Matrix{Vector{Segment}}(undef, Γ4.model.bandnum, Γ4.patchnum)
    #bidx-能带的指标，nidx-patch的指标
    for idxs in CartesianIndices(posimat)
        (bidx, nidx) = Tuple(idxs)
        posimat[bidx, nidx] = []
        negamat[bidx, nidx] = []
    end
    #
    for bidx in 1:1:Γ4.model.bandnum
        posi, peidx = const_energy_line_in_patches(
            Γ4.ltris, Γ4.ladjs, Γ4.lpats, lamb, Γ4.model.dispersion[bidx]
        )
        #将正能量的线加到对应的位置上
        for (edg, pidx) in zip(posi, peidx)
            push!(posimat[bidx, pidx], edg)
        end
        #
        nega, neidx = const_energy_line_in_patches(
            Γ4.ltris, Γ4.ladjs, Γ4.lpats, -lamb, Γ4.model.dispersion[bidx]
        )
        #将负能量的线加到对应的位置上
        for (edg, pidx) in zip(nega, neidx)
            push!(negamat[bidx, pidx], edg)
        end
    end
    #计算pp散射的bubble
    bubbval = Array{Float64, 7}(
        undef,
        #alpha, beta, b1, b2
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k1(b1), k2(b2)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    Threads.@threads for idxs in CartesianIndices(bubbval)
        alpha, beta, b1, b2, i_n, i1, i2 = Tuple(idxs)
        k1, k2 = Γ4.patches[b1][i1], Γ4.patches[b2][i2]
        qpp = Γ4.model.kadd(k1, k2)
        bubbres = pi_αβ_minus_ec(
            posimat[alpha, i_n], negamat[alpha, i_n],
            brlu_area, lamb,
            qpp, Γ4.model.dispersion[beta]
        )
        bubbval[alpha, beta, b1, b2, i_n, i1, i2] = bubbres
    end
    bubb_qpp = Bubble{:ec, :minus}(lval, bubbval)
    #计算fs和-fs
    bubbval_fs = Array{Float64, 7}(
        undef,
        #alpha, beta, b2, b3
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k2(b2), k3(b3)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    bubbval_nfs = Array{Float64, 7}(
        undef,
        #alpha, beta, b2, b3
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k2(b2), k3(b3)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    Threads.@threads for idxs in CartesianIndices(bubbval_fs)
        alpha, beta, b2, b3, i_n, i2, i3 = Tuple(idxs)
        k2, k3 = Γ4.patches[b2][i2], Γ4.patches[b3][i3]
        qfs = Γ4.model.kadd(k3, -k2)
        nqfs = -qfs
        bubbres_fs = pi_αβ_plus_ec(
            posimat[alpha, i_n], negamat[alpha, i_n],
            brlu_area, lamb,
            qfs, Γ4.model.dispersion[beta]
        )
        bubbres_nfs = pi_αβ_plus_ec(
            posimat[alpha, i_n], negamat[alpha, i_n],
            brlu_area, lamb,
            nqfs, Γ4.model.dispersion[beta]
        )
        bubbval_fs[alpha, beta, b2, b3, i_n, i2, i3] = bubbres_fs
        bubbval_nfs[alpha, beta, b2, b3, i_n, i2, i3] = bubbres_nfs
    end
    bubb_qfs = Bubble{:ec, :plus}(lval, bubbval_fs)
    bubb_nqfs = Bubble{:ec, :plus}(lval, bubbval_nfs)
    #计算ex和-ex
    bubbval_ex = Array{Float64, 7}(
        undef,
        #alpha, beta, b1, b3
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k1(b1), k3(b3)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    bubbval_nex = Array{Float64, 7}(
        undef,
        #alpha, beta, b1, b3
        Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum, Γ4.model.bandnum,
        #nidx(alpha), k1(b1), k3(b3)
        Γ4.patchnum, Γ4.patchnum, Γ4.patchnum
    )
    Threads.@threads for idxs in CartesianIndices(bubbval_ex)
        alpha, beta, b1, b3, i_n, i1, i3 = Tuple(idxs)
        k1, k3 = Γ4.patches[b1][i1], Γ4.patches[b3][i3]
        qex = Γ4.model.kadd(k1, -k3)
        nqex = -qex
        bubbres_ex = pi_αβ_plus_ec(
            posimat[alpha, i_n], negamat[alpha, i_n],
            brlu_area, lamb,
            qex, Γ4.model.dispersion[beta]
        )
        bubbres_nex = pi_αβ_plus_ec(
            posimat[alpha, i_n], negamat[alpha, i_n],
            brlu_area, lamb,
            nqex, Γ4.model.dispersion[beta]
        )
        bubbval_ex[alpha, beta, b1, b3, i_n, i1, i3] = bubbres_ex
        bubbval_nex[alpha, beta, b1, b3, i_n, i1, i3] = bubbres_nex
    end
    bubb_qex = Bubble{:ec, :plus}(lval, bubbval_ex)
    bubb_nqex = Bubble{:ec, :plus}(lval, bubbval_nex)
    return bubb_qpp, bubb_qfs, bubb_nqfs, bubb_qex, bubb_nqex
end



"""
使用能量cutoff作为flow parameter的bubble\n
posi是能量为+LAMBDA的边，nega是能量为-LAMBDA的边, lamb是LAMBDA\n
disp是色散关系，注意这个需要是beta带的，qval是需要平移的大小，应该用一个Point来包装，\n
kshf是动量相加的函数, 这个函数应该能处理好到第一布里渊区的映射\n
area是第一布里渊区的面积\n
```(10.112)本身已经处理好了动量守恒，k, k-q是需要满足动量守恒的关系的，而处理好```
```k-q到第一布里渊区的映射就处理好了Umklapp```
"""
function pi_αβ_plus_ec(
    posi::Vector{Basics.Segment}, nega::Vector{Basics.Segment},
    area::Float64, lamb::Float64, 
    qval::Point2D, 
    disp::Function
)
    """
    10.112中的 PI^+(n, q) = +LAMBDA (2pi)^-2 beta^-1 Int_{k in k_n} G'(k)G(k - Q)
    其中有一个beta是频率积分带来的，2pi^2是动量积分带来的
    G(k)=CITA(LAMBDA < abs(disp(k))) / i*omega - disp(k)
    G'(k)=-DELTA(abs(disp(k))-LAMBDA) / i*omege - disp(k)
    在零温的情况下10.112中的频率部分可以积分出来，此后的k都是不包含频率的
    = +LAMBDA (2pi)^-2 Int_{k in k_n} CITA() -DELTA()
    { beta^-1 sum_{omega} [(i*omega-disp(k))(i*omega-disp(k - q))]^-1 }
    花括号中的内容求和完之后等于 - CITA(-disp(k)disp(k-q)) / (abs(disp(k)) + abs(disp(k-p)))
    积分会变成
    = +LAMBDA (2pi)^-2 Int_{k in k_n} DELTA(abs(disp(k))-LAMBDA) CITA(LAMBDA<abs(disp(k-q)))
    CITA(-disp(k)disp(k-q)) / (abs(disp(k)) + abs(disp(k-p)))
    因为采用的能量cutoff中有一个 DELTA(abs(disp(k))-LAMBDA)，disp(k)等于正的或者负的LAMBDA
    而CITA(-disp(k)disp(k-q))限制了disp(k)和disp(k-q)符号相反
    所以上式变成
    (第一项disp(k)=LAMBDA>0，于是disp(k-q)<0，而且abs(disp(k))=-disp(k)>LAMBDA)
    (第二项类似，分子中的abs(disp(k))都可以直接换成LAMBDA，abs(disp(k-q))也都知道符号)
    = +LAMBDA (2pi)^-2 Int_{k in kn} {
        DELTA(disp(k)-LAMBDA)CITA(-disp(k-q)-LAMBDA) / (LAMBDA - disp(k - q))
        DELTA(disp(k)+LAMBDA)CITA(disp(k-q)-LAMBDA) / (LAMBDA + disp(k - q)) }
    还可以从积分里面把DELTA给积分掉，这样对于二维平面的积分也会变成对
    disp(k) = LAMBDA 或者 -LAMBDA的线的积分
    = +LAMBDA (2pi)^-2 *
    [Int_{disp(k) = +LAMBDA} CITA(-disp(k-q)-LAMBDA) / (LAMBDA - disp(k - q))]
    +[Int_{disp(k) = -LAMBDA} CITA(disp(k-q)-LAMBDA) / (LAMBDA + disp(k - q))  ]
    """
    nega_q = -qval
    #积分正LAMBDA的线
    intposi = 0.
    for edg in posi
        kval = middle_point(edg.pt1, edg.pt2)
        kprim = kval + nega_q
        #kprim = kadd(kval, nega_q)
        #CITA
        disp_kprim = disp(kprim.x, kprim.y)
        if -disp_kprim < lamb
            continue
        end
        #线积分，计算线元的长度
        intposi += edg.length / (lamb - disp_kprim)
    end
    #积分负LAMBDA的线
    intnega = 0.
    for edg in nega
        kval = middle_point(edg.pt1, edg.pt2)
        kprim = kval + nega_q
        #kprim = kadd(kval, nega_q)
        #CITA
        disp_kprim = disp(kprim.x, kprim.y)
        if disp_kprim < lamb
            continue
        end
        intnega += edg.length / (lamb + disp_kprim)
    end
    #乘上系数
    result = lamb * (intposi + intnega) / area
    return result
end



"""
使用能量cutoff作为flow parameter的bubble\n
posi是能量为+LAMBDA的边，nega是能量为-LAMBDA的边, lamb是LAMBDA\n
disp是色散关系，注意这个需要是beta带的，qval是需要平移的大小，应该用一个Point来包装，\n
kshf是动量相加的函数, 这个函数应该能处理好到第一布里渊区的映射\n
area是第一布里渊区的面积\n
```(10.112)本身已经处理好了动量守恒，k, k-q是需要满足动量守恒的关系的，而处理好```
```k-q到第一布里渊区的映射就处理好了Umklapp```
因为算能量的时候，超出第一布里渊区也没有关系，所以不再需要kadd
"""
function pi_αβ_minus_ec(
    posi::Vector{Basics.Segment}, nega::Vector{Basics.Segment},
    area::Float64, lamb::Float64, 
    qval::Point2D, 
    disp::Function
)
    """
    10.112中的 PI^-(n, q) = -LAMBDA (2pi)^-2 beta^-1 Int_{k in k_n} G'(k)G(- k + Q)
    = -LAMBDA (2pi)^-2 Int_{k in k_n} CITA() -DELTA()
    { beta^-1 sum_{omega} [(i*omega-disp(k))(-i*omega-disp(-k + q))]^-1 }
    在零温下这个频率积分等于，注意-k那里把频率也给反过来了
    +CITA(+disp(k)disp(-k+q)) / (abs(disp(k)) + abs(disp(-k+q)))
    原式就等于
    = LAMBDA (2pi)^-2 Int_{k in k_n} {
        DELTA(abs(disp(k))-LAMBDA) CITA(abs(disp(-k+q)-LAMBDA))
        CITA(disp(k)disp(-k+q)) / (abs(disp(k)) + abs(disp(-k+q))) }
    第二个CITA限制了disp(k)和disp(-k+q)同号，积分积掉DELTA，分类讨论正负
    = LAMBDA (2pi)^-2 {
        Int_{disp(k) = +LAMBDA} CITA(disp(-k+q) - LAMBDA) / (LAMBDA + disp(-k+q)) +
        Int_{disp(k) = -LAMBDA} CITA(-disp(-k+q) -LAMBDA) / (LAMBDA - disp(-k+q))
    }
    """
    #积分正LAMBDA的线
    intposi = 0.
    for edg in posi
        kval = middle_point(edg.pt1, edg.pt2)
        kprim = qval - kval
        #kprim = kadd(-kval, qval)
        #CITA
        disp_kprim = disp(kprim.x, kprim.y)
        if disp_kprim < lamb
            continue
        end
        #线积分，计算线元的长度
        intposi += edg.length / (lamb + disp_kprim)
    end
    #积分负LAMBDA的线
    intnega = 0.
    for edg in nega
        kval = middle_point(edg.pt1, edg.pt2)
        kprim = qval - kval
        #kprim = kadd(-kval, qval)
        #CITA
        disp_kprim = disp(kprim.x, kprim.y)
        if -disp_kprim < lamb
            continue
        end
        intnega += edg.length / (lamb - disp_kprim)
    end
    #乘上系数
    result = lamb * (intposi + intnega) / area
    return result
end




