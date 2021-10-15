"""
二维模型的抽象类型
"""


module Fermi
    
abstract type Abstract2DModel{T} end


using ..Basics

export QuadrateSystem, TriangularSystem

"""
正方晶系
"""
struct QuadrateSystem{T} <: Abstract2DModel{T}
    brillouin :: Basics.AbstractRectangle{:AXISSQUARE}
    dispersion :: Vector{Function}
    bandnum :: Int64
    kadd :: Function
    QuadrateSystem{T}(disps) where T = begin
        new{T}(
            Square(Point2D(0., 0.), 2pi),
            disps,
            length(disps),
            square_kadd
        )
    end
end


"""
正方晶系的动量平移
"""
function square_kadd(pt1::Point2D, pt2::Point2D)
    destx = pt1.x + pt2.x
    desty = pt1.y + pt2.y
    #如果大于pi，就像0靠拢
    while abs(destx) > pi
        destx = destx - sign(destx) * 2pi
    end
    while abs(desty) > pi
        desty = desty - sign(desty) * 2pi
    end
    return Point2D(destx, desty)
end


"""
三角晶系，这里三角晶系的的初基矢量的选择
a1=(-1/2, √3/2); a2=(1/2, √3/2);
b1=π(-2, 2/√3); b2=π(2, 2/√3);
"""
struct TriangularSystem{T} <: Abstract2DModel{T}
    brillouin :: Basics.AbstractHexagon{:EQ}
    dispersion :: Vector{Function}
    bandnum :: Int64
    kadd :: Function
    TriangularSystem{T}(disps) where T = begin
        new{T}(
            EqHexagon(Point2D(0., 0.), 4pi/3.0),
            disps,
            length(disps),
            hexagon_kadd
        )
    end
end


"""
三角晶系的动量平移
"""
function hexagon_kadd(pt1::Point2D, pt2::Point2D)
    #找到离目标最近的一个第一布里渊区的中心
    #如果这个点距离中心是最近的，那么他就在第一布里渊区里面了，如果
    #离其他的点近，就平移一次，再检查一遍
    b1x = -2pi
    b1y = 2pi/sqrt(3)
    brlu_cnts = Vector{Point2D}(undef, 7)
    brlu_cnts[1] = Point2D(0., 0.)
    brlu_cnts[2] = Point2D(b1x, b1y)
    brlu_cnts[3] = Point2D(0., 2*b1y)
    brlu_cnts[4] = Point2D(-b1x, b1y)
    brlu_cnts[5] = Point2D(-b1x, -b1y)
    brlu_cnts[6] = Point2D(0., -2*b1y)
    brlu_cnts[7] = Point2D(b1x, -b1y)
    #
    dest = pt1 + pt2
    dis2cents = [dest - cnt for cnt in brlu_cnts]
    dis2cents = [pt.x^2 + pt.y^2 for pt in dis2cents]
    #
    (minv, mina) = findmin(dis2cents)
    while mina != 1
        dest = dest - brlu_cnts[mina]
        dis2cents = [dest - cnt for cnt in brlu_cnts]
        dis2cents = [pt.x^2 + pt.y^2 for pt in dis2cents]
        (minv, mina) = findmin(dis2cents)
    end
    return dest
end



#
#一些内建的模型
#不能在include上面加"""document"""

include("quadrate_lattice.jl")
export common_square_lattice


include("triangular_lattice.jl")
export common_triangle_lattice


include("kagome_lattice.jl")
export upperband_kagome_lattice


include("surface.jl")


include("patch.jl")

end


