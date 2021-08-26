"""
三角晶系相关
"""


"""
普通的三角格子
"""
function common_triangle_lattice(μ)
    disp(x, y) = -2cos(x) - 4cos(sqrt(3)*y/2.0)*cos(x/2.0) + μ
    return TriangularSystem{:COM}(
        [disp]
    )
end

