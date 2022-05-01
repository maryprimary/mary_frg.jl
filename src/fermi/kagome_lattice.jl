"""
kagome lattice的相关功能
"""

using ..Fermi
#using ..Fermi.Patch

using LinearAlgebra
#using GenericLinearAlgebra

export upperband_kagome_lattice
export wuxx_av3sb5_kagome_lattice
#export upperband_kagome_lattice2
export get_kagome_ν, get_kagome_U_mt
export get_wuxx_U_mt

export dispersion


"""
上面那个能带
"""
macro upperband(kx, ky)
    return esc(quote
        xval = $kx / 4
        yval = √3 * ($ky) / 4
        sqr = 2 * cos(2*xval - 2*yval)
        sqr += 2 * cos(2*xval + 2*yval)
        sqr += 2 * cos(4*xval) + 3
        if sqr < 0
            sqr = 0
        end
        eng = - 1 + sqrt(sqr)
    end)
end


"""
下面那个能带
"""
macro lowerband(kx, ky)
    return esc(quote
        xval = $kx / 4
        yval = √3 * ($ky) / 4
        sqr = 2 * cos(2*xval - 2*yval)
        sqr += 2 * cos(2*xval + 2*yval)
        sqr += 2 * cos(4*xval) + 3
        if sqr < 0
            sqr = 0
        end
        eng = - 1 - sqrt(sqr)
    end) 
end


"""
上半个能带的kagome lattice
"""
function upperband_kagome_lattice(μ)
    #disp(x, y) = (@upperband x y) + μ
    return TriangularSystem{:KAGOME}(
        [μ]
    )
end


"""
单带kagome的色散
"""
function dispersion(model::TriangularSystem{:KAGOME}, ::Int64, x, y)
    return (@upperband x y) + model.μs[1]
end


"""
两个能带的A V_3 Sb_5
"""
function wuxx_av3sb5_kagome_lattice()
    #disp1(x, y) = 0.5*(@upperband x y) - 0.055
    #disp2(x, y) = (@lowerband x y) + 2.182
    return TriangularSystem{:WUXX}(
        [-0.055, 2.182]
    )
end


"""
wuxxmodel的色散
"""
function dispersion(model::TriangularSystem{:WUXX}, bidx::Int64, x, y)
    if bidx == 1
        return 0.5*(@upperband x y) + model.μs[1]
    else
        return (@lowerband x y) + model.μs[2]
    end
end



"""
获取nu变换的矩阵
子格子位置是
 B
A C
"""
function get_kagome_ν(kx, ky)
    if kx == 0. && ky == 0.
        kx = 1e-8
        ky = -1e-8
    end
    #x::BigFloat = BigFloat(kx) / BigFloat(4)
    #y::BigFloat = sqrt(BigFloat(3)) * BigFloat(ky) / BigFloat(4)
    #mat = zeros(3, 3)
    #mat[1, 2] = 2cos(x+y)
    #mat[1, 3] = 2cos(2x)
    #mat[2, 3] = 2cos(-x+y)
    #mat[2, 1] = 2cos(x+y)
    #mat[3, 1] = 2cos(2x)
    #mat[3, 2] = 2cos(-x+y)
    mat = zeros(3, 3)
    mat[1, 2] = 2cos(kx/4+√3*ky/4)#1 + exp(-im*(kx/2+√3*ky/2))
    mat[1, 3] = 2cos(kx/2)#1 + exp(-im*(kx))
    mat[2, 3] = 2cos(kx/4-√3*ky/4)#1 + exp(-im*(kx/2 - √3*ky/2))
    mat[2, 1] = 2cos(kx/4+√3*ky/4)#1 + exp(im*(kx/2+√3*ky/2))
    mat[3, 1] = 2cos(kx/2)#1 + exp(im*(kx))
    mat[3, 2] = 2cos(kx/4-√3*ky/4)#1 + exp(im*(kx/2 - √3*ky/2))
    evals, evecs = eigen(mat)
    #println(transpose(evecs)*mat*evecs)
    #println(evals)
    nu1 = evecs[:, 1]
    if abs(minimum(nu1)) > maximum(nu1) + eps()
        nu1 = -nu1
    end
    nu2 = evecs[:, 2]
    if abs(minimum(nu2)) > maximum(nu2) + eps()
        nu2 = -nu2
    end
    nu3 = evecs[:, 3]
    if abs(minimum(nu3)) > maximum(nu3) + eps()
        nu3 = -nu3
    end
    return nu1, nu2, nu3
    #sqr::BigFloat = 2 * cos(2*x - 2*y)
    #sqr::BigFloat += 2 * cos(2*x + 2*y)
    #sqr::BigFloat += 2 * cos(4*x) + 3
    #if sqr < 0
    #    sqr = 0
    #end
    ##第一组
    #nu1::Vector{BigFloat} = zeros(BigFloat, 3)
    #nu1[1] = -cos(2*x)/2 + cos(2*y)/2
    #nu1[2] = -cos(x-y)/2 + cos(3*x+y)/2
    #nu1[3] = 1 - cos(x+y)^2
    ##归一化
    #n1 = nu1[1]^2 + nu1[2]^2 + nu1[3]^2
    #n1 = 1.0 / sqrt(n1)
    ##println(nu1)
    #nu1 = n1 * nu1
    ##println(nu1)
    ##println("=======")
    ##
    #nu2::Vector{BigFloat} = zeros(BigFloat, 3)
    ##第二组本征态
    #lamb2::BigFloat = 1 - sqrt(sqr)
    #nu2[1] = 0.5*(lamb2^2) - 2*(cos(x-y)^2)
    #nu2[2] = 2*cos(2*x)*cos(x-y) + cos(x+y)*lamb2
    #nu2[3] = cos(2*x)*lamb2 + 2*cos(x-y)*cos(x+y)
    ##归一化
    ##println(nu2)
    #n2 = nu2[1]^2 + nu2[2]^2 + nu2[3]^2
    ##n2 = 1.0 / sqrt(n2)
    ##nu2 = n2 * nu2
    #if abs(minimum(nu2)) > maximum(nu2)
    #    nu2 = -nu2
    #end
    ##println(nu2)
    ##println("**********")
    ##
    #nu3::Vector{BigFloat} = zeros(BigFloat, 3)
    ##第三组本征态
    #lamb3::BigFloat = 1 + sqrt(sqr)
    #nu3[1] = 0.5*(lamb3^2) - 2*(cos(x-y)^2)
    #nu3[2] = 2*cos(2*x)*cos(x-y) + cos(x+y)*lamb3
    #nu3[3] = cos(2*x)*lamb3 + 2*cos(x-y)*cos(x+y)
    ##归一化
    ##println(nu3)
    #n3 = nu3[1]^2 + nu3[2]^2 + nu3[3]^2
    #n3 = 1.0 / sqrt(n3)
    #nu3 = n3 * nu3
    ##println(nu3)
    ##println("~~~~~~~~")
    #return nu1, nu2, nu3
end



"""
获取相互作用在能带下的表示
"""
function get_kagome_U_mt(Γ4, uval)
    pinfos = Γ4.patches[1]
    npat = length(pinfos)
    ret = zeros(npat, npat, npat)
    @Threads.threads for idxs in CartesianIndices(ret)
        k1i, k2i, k3i = Tuple(idxs)
        k1v = pinfos[k1i]
        k2v = pinfos[k2i]
        k3v = pinfos[k3i]
        #k4v = Fermi.hexagon_kadd2(k1v, k2v)
        #k4v = Fermi.hexagon_kadd2(k4v, -k3v)
        k4i = Γ4.k4tab[1, 1, 1, 1, k1i, k2i, k3i]
        k4v = pinfos[k4i]
        #p这个带只需要第二个
        _, k1nu, _ = get_kagome_ν(k1v.x, k1v.y)
        _, k2nu, _ = get_kagome_ν(k2v.x, k2v.y)
        _, k3nu, _ = get_kagome_ν(k3v.x, k3v.y)
        _, k4nu, _ = get_kagome_ν(k4v.x, k4v.y)
        #计算数值
        for sidx in 1:1:3
            ret[k1i, k2i, k3i] += uval * k1nu[sidx] *
            k2nu[sidx] * k3nu[sidx] * k4nu[sidx]
        end
    end
    return ret
end


"""
获取A V_3 Sb_5相互作用在能带下的表示
"""
function get_wuxx_U_mt(Γ4, u1val, u2val, upval; usesymm=true)
    pinfos1 = Γ4.patches[1]
    pinfos2 = Γ4.patches[2]
    npat = length(pinfos1)
    if length(pinfos2) != npat
        throw(error("pat数量对不上"))
    end
    #
    ret = zeros(2, 2, 2, 2, npat, npat, npat)
    place_holder = Array{Int8, 3}(undef, npat, npat, npat)
    #六度对称
    symmsector = Int64(npat // 6)
    if usesymm
        place_holder = Array{Int8, 3}(undef, symmsector, npat, npat)
    end
    #能带内的相互作用
    #第一条能带是靠上的能带
    @Threads.threads for idxs in CartesianIndices(place_holder)
        k1i, k2i, k3i = Tuple(idxs)
        k1v = pinfos1[k1i]
        k2v = pinfos1[k2i]
        k3v = pinfos1[k3i]
        #k4v = Fermi.hexagon_kadd2(k1v, k2v)
        #k4v = Fermi.hexagon_kadd2(k4v, -k3v)
        k4i = Γ4.k4tab[1, 1, 1, 1, k1i, k2i, k3i]#Patch.find_patch_index_hexa(k4v, EqHexagon(Point2D(0., 0.), 4pi/3.0), npat)
        k4v = pinfos1[k4i]
        #p这个带只需要第二个
        _, k1nu, _ = get_kagome_ν(k1v.x, k1v.y)
        _, k2nu, _ = get_kagome_ν(k2v.x, k2v.y)
        _, k3nu, _ = get_kagome_ν(k3v.x, k3v.y)
        _, k4nu, _ = get_kagome_ν(k4v.x, k4v.y)
        #计算数值
        for sidx in 1:1:3
            ret[1, 1, 1, 1, k1i, k2i, k3i] += u1val * k1nu[sidx] *
            k2nu[sidx] * k3nu[sidx] * k4nu[sidx]
        end
    end
    #第二条是靠下的能带
    @Threads.threads for idxs in CartesianIndices(place_holder)
        k1i, k2i, k3i = Tuple(idxs)
        k1v = pinfos2[k1i]
        k2v = pinfos2[k2i]
        k3v = pinfos2[k3i]
        #k4v = Fermi.hexagon_kadd2(k1v, k2v)
        #k4v = Fermi.hexagon_kadd2(k4v, -k3v)
        #k4i = Patch.find_patch_index_hexa(k4v, EqHexagon(Point2D(0., 0.), 4pi/3.0), npat)
        k4i = Γ4.k4tab[2, 2, 2, 2, k1i, k2i, k3i]
        k4v = pinfos2[k4i]
        #d这个带只需要第三个
        _, _, k1nu = get_kagome_ν(k1v.x, k1v.y)
        _, _, k2nu = get_kagome_ν(k2v.x, k2v.y)
        _, _, k3nu = get_kagome_ν(k3v.x, k3v.y)
        _, _, k4nu = get_kagome_ν(k4v.x, k4v.y)
        #计算数值
        for sidx in 1:1:3
            ret[2, 2, 2, 2, k1i, k2i, k3i] += u2val * k1nu[sidx] *
            k2nu[sidx] * k3nu[sidx] * k4nu[sidx]
        end
    end
    #能带之间的相互作用
    @Threads.threads for idxs in CartesianIndices(place_holder)
        k1i, k2i, k3i = Tuple(idxs)
        #第一个能带的patches
        k1v1 = pinfos1[k1i]
        k2v1 = pinfos1[k2i]
        k3v1 = pinfos1[k3i]
        #第一个能带的系数
        _, k11nu, _ = get_kagome_ν(k1v1.x, k1v1.y)
        _, k21nu, _ = get_kagome_ν(k2v1.x, k2v1.y)
        _, k31nu, _ = get_kagome_ν(k3v1.x, k3v1.y)
        #第二个能带的patches
        k1v2 = pinfos2[k1i]
        k2v2 = pinfos2[k2i]
        k3v2 = pinfos2[k3i]
        #第二个能带的系数
        _, _, k12nu = get_kagome_ν(k1v2.x, k1v2.y)
        _, _, k22nu = get_kagome_ν(k2v2.x, k2v2.y)
        _, _, k32nu = get_kagome_ν(k3v2.x, k3v2.y)
        #两个k4和系数
        #k4v1 = Fermi.hexagon_kadd2(k1v1, k2v2)
        #k4v1 = Fermi.hexagon_kadd2(k4v1, -k3v2)
        #k4i1 = Patch.find_patch_index_hexa(k4v1, EqHexagon(Point2D(0., 0.), 4pi/3.0), npat)
        k4i1 = Γ4.k4tab[1, 2, 2, 1, k1i, k2i, k3i]
        k4v1 = pinfos1[k4i1]
        _, k41nu, _ = get_kagome_ν(k4v1.x, k4v1.y)
        #k4v2 = Fermi.hexagon_kadd2(k1v2, k2v1)
        #k4v2 = Fermi.hexagon_kadd2(k4v2, -k3v1)
        #k4i2 = Patch.find_patch_index_hexa(k4v2, EqHexagon(Point2D(0., 0.), 4pi/3.0), npat)
        k4i2 = Γ4.k4tab[2, 1, 1, 2, k1i, k2i, k3i]
        k4v2 = pinfos2[k4i2]
        _, _, k42nu = get_kagome_ν(k4v2.x, k4v2.y)
        #计算数值
        #对于同能带U来说，自旋的求和up，dn和dn，up重复两次，再除2，正好一倍
        #对于不同能带U来说，自旋求和本来就有，能带指标1221,2112的重复刚好和除2抵消。
        for sidx in 1:1:3
            ret[1, 2, 2, 1, k1i, k2i, k3i] += upval * k11nu[sidx] *
            k22nu[sidx] * k32nu[sidx] * k41nu[sidx]
            ret[2, 1, 1, 2, k1i, k2i, k3i] += upval * k12nu[sidx] *
            k21nu[sidx] * k31nu[sidx] * k42nu[sidx]
        end
    end
    #不利用对称性则直接返回
    if !usesymm
        return ret
    end
    #利用对称性
    symm_holder = Array{Int8, 3}(undef, 5, npat, npat)
    @Threads.threads for idxs in CartesianIndices(symm_holder)
        sec, k2i, k3i = Tuple(idxs)
        offset = sec*symmsector
        nk2i = k2i - offset
        nk2i = nk2i < 1 ? nk2i + npat : nk2i
        nk3i = k3i - offset
        nk3i = nk3i < 1 ? nk3i + npat : nk3i
        #
        ret[1, 1, 1, 1, 1+offset:symmsector+offset, k2i, k3i] =
        ret[1, 1, 1, 1, 1:symmsector, nk2i, nk3i]
        #
        ret[2, 2, 2, 2, 1+offset:symmsector+offset, k2i, k3i] =
        ret[2, 2, 2, 2, 1:symmsector, nk2i, nk3i]
        #
        ret[1, 2, 2, 1, 1+offset:symmsector+offset, k2i, k3i] =
        ret[1, 2, 2, 1, 1:symmsector, nk2i, nk3i]
        #
        ret[2, 1, 1, 2, 1+offset:symmsector+offset, k2i, k3i] =
        ret[2, 1, 1, 2, 1:symmsector, nk2i, nk3i]
        #
    end
    retnew = zeros(2, 2, 2, 2, npat, npat, npat)
    #反转对称性
    @Threads.threads for idxs in CartesianIndices(retnew)
        b1, b2, b3, b4, i1, i2, i3 = Tuple(idxs)
        i4 = Γ4.symmtab[b1, b2, b3, b4, i1, i2, i3]
        if i4 == -1
            retnew[b1, b2, b3, b4, i1, i2, i3] = ret[b1, b2, b3, b4, i1, i2, i3]
        else
            retnew[b1, b2, b3, b4, i1, i2, i3] =  0.25 * (
                ret[b1, b2, b3, b4, i1, i2, i3] +
                ret[b4, b3, b2, b1, i4, i3, i2] +
                ret[b2, b1, b4, b3, i2, i1, i4] +
                ret[b3, b4, b1, b2, i3, i4, i1]
            )
        end
    end
    return retnew
end


"""
wuxx格子上的近邻V
1/2 (n_{1i} n_{1j} + n_{1j}n_{1i} + n_{2i} n_{2j} + n_{2j}n_{2i}
n_{1i} n_{2j} + n_{2j}n_{1i} + n_{1i} n_{2j} + n_{2j}n_{1i}
)
"""
function get_wuxx_V_mt()
    
end

