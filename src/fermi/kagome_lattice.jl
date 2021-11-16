"""
kagome lattice的相关功能
"""

module Kagome

using ..Fermi

export upperband_kagome_lattice
export wuxx_av3sb5_kagome_lattice
#export upperband_kagome_lattice2
export get_kagome_ν, get_kagome_U_mt
export get_wuxx_U_mt


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
    disp(x, y) = (@upperband x y) + μ
    return TriangularSystem{:KAGOME}(
        [disp]
    )
end


"""
两个能带的A V_3 Sb_5
"""
function wuxx_av3sb5_kagome_lattice()
    disp1(x, y) = 0.5*(@upperband x y) - 0.055
    disp2(x, y) = (@lowerband x y) + 2.182
    return TriangularSystem{:WUXX}(
        [disp1, disp2]
    )
end



"""
获取nu变换的矩阵
"""
function get_kagome_ν(kx, ky)
    if kx == 0. && ky == 0.
        kx = 1e-10
        ky = -1e-10
    end
    x::BigFloat = kx / BigFloat(4)
    y::BigFloat = BigFloat(√3) * ky / BigFloat(4)
    sqr::BigFloat = 2 * cos(2*x - 2*y)
    sqr::BigFloat += 2 * cos(2*x + 2*y)
    sqr::BigFloat += 2 * cos(4*x) + 3
    if sqr < 0
        sqr = 0
    end
    #第一组
    nu1::Vector{BigFloat} = zeros(BigFloat, 3)
    nu1[1] = -cos(2*x)/2 + cos(2*y)/2
    nu1[2] = -cos(x-y)/2 + cos(3*x+y)/2
    nu1[3] = 1 - cos(x+y)^2
    #归一化
    n1 = nu1[1]^2 + nu1[2]^2 + nu1[3]^2
    n1 = 1.0 / sqrt(n1)
    #println(nu1)
    nu1 = n1 * nu1
    #println(nu1)
    #println("=======")
    #
    nu2::Vector{BigFloat} = zeros(BigFloat, 3)
    #第二组本征态
    lamb2::BigFloat = 1 - sqrt(sqr)
    nu2[1] = 0.5*(lamb2^2) - 2*(cos(x-y)^2)
    nu2[2] = 2*cos(2*x)*cos(x-y) + cos(x+y)*lamb2
    nu2[3] = cos(2*x)*lamb2 + 2*cos(x-y)*cos(x+y)
    #归一化
    #println(nu2)
    n2 = nu2[1]^2 + nu2[2]^2 + nu2[3]^2
    n2 = 1.0 / sqrt(n2)
    nu2 = n2 * nu2
    #println(nu2)
    #println("**********")
    #
    nu3::Vector{BigFloat} = zeros(BigFloat, 3)
    #第三组本征态
    lamb3::BigFloat = 1 + sqrt(sqr)
    nu3[1] = 0.5*(lamb3^2) - 2*(cos(x-y)^2)
    nu3[2] = 2*cos(2*x)*cos(x-y) + cos(x+y)*lamb3
    nu3[3] = cos(2*x)*lamb3 + 2*cos(x-y)*cos(x+y)
    #归一化
    #println(nu3)
    n3 = nu3[1]^2 + nu3[2]^2 + nu3[3]^2
    n3 = 1.0 / sqrt(n3)
    nu3 = n3 * nu3
    #println(nu3)
    #println("~~~~~~~~")
    return nu1, nu2, nu3
end



"""
获取相互作用在能带下的表示
"""
function get_kagome_U_mt(uval, pinfos)
    npat = length(pinfos)
    ret = zeros(npat, npat, npat)
    @Threads.threads for idxs in CartesianIndices(ret)
        k1i, k2i, k3i = Tuple(idxs)
        k1v = pinfos[k1i]
        k2v = pinfos[k2i]
        k3v = pinfos[k3i]
        k4v = Fermi.hexagon_kadd2(k1v, k2v)
        k4v = Fermi.hexagon_kadd2(k4v, -k3v)
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
function get_wuxx_U_mt(u1val, u2val, upval, pinfos1, pinfos2)
    npat = length(pinfos1)
    if length(pinfos2) != npat
        throw(error("pat数量对不上"))
    end
    #
    ret = zeros(2, 2, 2, 2, npat, npat, npat)
    place_holder = Array{Int8, 3}(undef, npat, npat, npat)
    #能带内的相互作用
    #第一条能带是靠上的能带
    @Threads.threads for idxs in CartesianIndices(place_holder)
        k1i, k2i, k3i = Tuple(idxs)
        k1v = pinfos1[k1i]
        k2v = pinfos1[k2i]
        k3v = pinfos1[k3i]
        k4v = Fermi.hexagon_kadd2(k1v, k2v)
        k4v = Fermi.hexagon_kadd2(k4v, -k3v)
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
        k4v = Fermi.hexagon_kadd2(k1v, k2v)
        k4v = Fermi.hexagon_kadd2(k4v, -k3v)
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
        k4v1 = Fermi.hexagon_kadd2(k1v1, k2v1)
        k4v1 = Fermi.hexagon_kadd2(k4v1, -k3v1)
        #第一个能带的系数
        _, k11nu, _ = get_kagome_ν(k1v1.x, k1v1.y)
        _, k21nu, _ = get_kagome_ν(k2v1.x, k2v1.y)
        _, k31nu, _ = get_kagome_ν(k3v1.x, k3v1.y)
        _, k41nu, _ = get_kagome_ν(k4v1.x, k4v1.y)
        #第二个能带的patches
        k1v2 = pinfos2[k1i]
        k2v2 = pinfos2[k2i]
        k3v2 = pinfos2[k3i]
        k4v2 = Fermi.hexagon_kadd2(k1v2, k2v2)
        k4v2 = Fermi.hexagon_kadd2(k4v2, -k3v2)
        #第二个能带的系数
        _, _, k12nu = get_kagome_ν(k1v2.x, k1v2.y)
        _, _, k22nu = get_kagome_ν(k2v2.x, k2v2.y)
        _, _, k32nu = get_kagome_ν(k3v2.x, k3v2.y)
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
    return ret
end



end # end module
