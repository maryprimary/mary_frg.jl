"""
正方晶系相关
"""

export common_square_lattice
export dispersion


"""
普通的正方格子
"""
function common_square_lattice(μ)
    #disp(x, y) = -2*(cos(x) + cos(y)) + μ
    return QuadrateSystem{:COM}(
        [μ]
    )
end


"""
普通正方格子的色散
"""
function dispersion(model::QuadrateSystem{:COM}, ::Int64, x, y)
    return -2*(cos(x) + cos(y)) + model.μs[1]
end

