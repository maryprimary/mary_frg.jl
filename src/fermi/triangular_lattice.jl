"""
三角晶系相关
"""

export common_triangle_lattice
export dispersion


"""
普通的三角格子
"""
function common_triangle_lattice(μ)
    #disp(x, y) = -2cos(x) - 4cos(sqrt(3)*y/2.0)*cos(x/2.0) + μ
    return TriangularSystem{:COM}(
        [μ]
    )
end


"""
普通三角格子的色散
"""
function dispersion(model::TriangularSystem{:COM}, ::Int64, x, y)
    return -2cos(x) - 4cos(sqrt(3)*y/2.0)*cos(x/2.0) + model.μs[1]
end

