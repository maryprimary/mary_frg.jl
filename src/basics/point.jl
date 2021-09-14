"""
基本的点
"""

using Printf

struct Point2D
    x :: Float64
    y :: Float64
end

#function Base.show(io::IO, pt::Point2D)
#    @printf(io, "2D Point (%.6f, %.6f)\n", pt.x, pt.y)
#end

function shift_point(pt1::Point2D, pt2::Point2D)
    return Point2D(pt1.x + pt2.x, pt1.y + pt2.y)
end

Base.:+(pt1::Point2D, pt2::Point2D) = shift_point(pt1, pt2)

Base.:-(pt1::Point2D) = Point2D(-pt1.x, -pt1.y)

Base.:-(pt1::Point2D, pt2::Point2D) = shift_point(pt1, -pt2)

Base.:*(scal::N, pt1::Point2D) where N <: Real = Point2D(scal*pt1.x, scal*pt1.y)

"""
两个点的中点
"""
function middle_point(pt1::Point2D, pt2::Point2D; sc1=0.5, sc2=0.5) :: Point2D
    scs = sc1 + sc2
    xpt = pt1.x * sc1 + pt2.x * sc2
    xpt /= scs
    ypt = pt1.y * sc1 + pt2.y * sc2
    ypt /= scs
    return Point2D(xpt, ypt)
end


"""获取坐标系中点的角度，从x轴逆时针开始算"""
function absolute_angle(pt::Point2D) 
    ang = atan(pt.y, pt.x)
    #负半轴加2pi
    if pt.y < 0
        ang = ang + 2pi
    end
    return ang
end
