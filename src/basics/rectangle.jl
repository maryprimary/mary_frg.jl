"""
长方形有关的功能
"""


struct AbstractRectangle{T} <: AbstractPolygon{T}
    vertex :: Array{Point2D, 1}
    edges :: Array{Segment, 1}
    center :: Point2D
    width :: Float64
    height :: Float64
end


"""
构造一个平行坐标轴的长方形
"""
function Rectangle(center::Point2D, width, height)
    verts = [
        Point2D(center.x - width/2, center.y + height/2),
        Point2D(center.x - width/2, center.y - height/2),
        Point2D(center.x + width/2, center.y - height/2),
        Point2D(center.x + width/2, center.y + height/2)
    ]
    edges = [
        Segment(verts[1], verts[2]),
        Segment(verts[2], verts[3]),
        Segment(verts[3], verts[4]),
        Segment(verts[4], verts[1])
    ]
    return AbstractRectangle{:AXISRECT}(verts, edges, center, width, height)
end


function Square(center::Point2D, width)
    verts = [
        Point2D(center.x - width/2, center.y + width/2),
        Point2D(center.x - width/2, center.y - width/2),
        Point2D(center.x + width/2, center.y - width/2),
        Point2D(center.x + width/2, center.y + width/2)
    ]
    edges = [
        Segment(verts[1], verts[2]),
        Segment(verts[2], verts[3]),
        Segment(verts[3], verts[4]),
        Segment(verts[4], verts[1])
    ]
    return AbstractRectangle{:AXISSQUARE}(verts, edges, center, width, width)
end


"""
长方型的面积
"""
function area(rect::AbstractRectangle{N}) where N
    return rect.height * rect.width
end


#display函数在中间加::MIME"text/plain"
function Base.show(io::IO, tri::AbstractRectangle{N}) where N
    name = Dict(
        :AXISRECT => "Rectangle",
        :AXISSQUARE => "Square"
    )
    @printf(io, "%s\n center: %s\n width: %.6f\n height: %.6f",
    name[N], tri.center, tri.width, tri.height)
end

