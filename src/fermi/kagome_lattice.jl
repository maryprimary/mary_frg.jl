"""
kagome lattice的相关功能
"""

module Kagome

using ..Fermi

export upperband_kagome_lattice
#export upperband_kagome_lattice2
export get_kagome_ν, get_kagome_U_mt


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
            if isapprox(sqr, 0., atol=1e-12)
                sqr = 0
            else
                throw(error("能带有错误"))
            end
        end
        eng = - 1 + sqrt(sqr)
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



end # end module
