"""
正方晶系相关
"""


"""
普通的正方格子
"""
function common_square_lattice(μ)
    disp(x, y) = -2*(cos(x) + cos(y)) + μ
    return QuadrateSystem{:COM}(
        [disp]
    )
end
