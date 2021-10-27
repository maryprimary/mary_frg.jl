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
            hexagon_kadd2
        )
    end
end


"""
三角晶系的动量平移
"""
function hexagon_kadd(pt1::Point2D, pt2::Point2D)
    @warn "hexagon_kadd很慢，用hexagon_kadd2代替"
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


"""
三角晶系的动量平移
b1=π(-2, 2/√3); b2=π(2, 2/√3);
"""
function hexagon_kadd2(pt1::Point2D, pt2::Point2D)
    dest = pt1 + pt2
    destx, desty = dest.x, dest.y
    #用一个长方型的框包起来
    while abs(desty) > π*2/√3
        desty = desty - sign(desty) * π*4/√3
    end
    while abs(destx) > 2*π
        destx = destx - sign(destx) * π*4
    end
    #如果在竖线上，肯定已经在第一布里渊区
    if isapprox(destx, 0.)
        return Point2D(destx, desty)
    end
    #如果在横线上
    if isapprox(desty, 0.)
        if destx > pi*4/3
            return Point2D(destx-π*2, π*2/√3)
        end
        if desty < -pi*4/3
            return Point2D(destx+π*2, π*2/√3)
        end
        return Point2D(destx, desty)
    end
    #如果在四个角落上
    lslope = -√3 * sign(destx) * sign(desty)
    y0 = sign(desty)*π*4/√3
    #println(lslope)
    #println(y0)
    #println(lslope*destx + y0)
    if desty*sign(desty) > (lslope*destx + y0)*sign(desty)
        destx = destx - sign(destx)*π*2
        desty = desty - sign(desty)*π*2/√3
    end
    return Point2D(destx, desty)
end



##"""不比kadd2要快
##三角晶系的动量平移
##b1=π(-2, 2/√3); b2=π(2, 2/√3);
##"""
##function hexagon_kadd3(pt1::Point2D, pt2::Point2D)
##    dest = pt1 + pt2
##    destx, desty = dest.x, dest.y
##    #用一个长方型的框包起来
##    while abs(desty) > π*2/√3
##        desty = desty - sign(desty) * π*4/√3
##    end
##    while abs(destx) > 2*π
##        destx = destx - sign(destx) * π*4
##    end
##    #如果在竖线上，肯定已经在第一布里渊区
##    if isapprox(destx, 0.)
##        return Point2D(destx, desty)
##    end
##    #
##    if sign(destx) * sign(desty) > 0
##        inp = 0.5*(√3)*destx + 0.5*desty
##        if inp > π*2/√3
##            return Point2D(destx - π*2, desty - π*2/√3)
##        elseif inp < -π*2/√3
##            return Point2D(destx + π*2, desty + π*2/√3)
##        end
##        return Point2D(destx, desty)
##    else
##        inp = -0.5*(√3)*destx + 0.5*desty
##        if inp > π*2/√3
##            return Point2D(destx + π*2, desty - π*2/√3)
##        elseif inp < -π*2/√3
##            return Point2D(destx - π*2, desty + π*2/√3)
##        end
##        return Point2D(destx, desty)
##    end
##end


#
#一些内建的模型
#不能在include上面加"""document"""

include("quadrate_lattice.jl")
export common_square_lattice


include("triangular_lattice.jl")
export common_triangle_lattice


include("kagome_lattice.jl")


include("surface.jl")


include("patch.jl")

end


