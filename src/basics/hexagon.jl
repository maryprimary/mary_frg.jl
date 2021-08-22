"""
六边型相关的功能
"""


struct AbstractHexagon{T} <: AbstractPolygon{T}
    vertex :: Array{Point2D, 1}
    edges :: Array{Segment, 1}
    center :: Point2D
end


"""
创建一个等边的六角型，平行与坐标轴
"""
function EqHexagon(center::Point2D, edgelen)
    cnx, cny = center.x, center.y
    xlim = edgelen / 2.0
    ylim = sqrt(3) * xlim
    verts = [
        Point2D(cnx+2xlim, cny),
        Point2D(cnx+xlim, cny+ylim),
        Point2D(cnx-xlim, cny+ylim),
        Point2D(cnx-2xlim, cny),
        Point2D(cnx-xlim, cny-ylim),
        Point2D(cnx+xlim, cny-ylim)
    ]
    edges = [Segment(verts[idx], verts[idx+1]) for idx in 1:1:5]
    push!(edges, Segment(verts[6], verts[1]))
    return AbstractHexagon{:EQ}(verts, edges, center)
end


"""
六边型的面积
"""
function area(hexa::AbstractHexagon{N}) where N
    if N == :EQ
        return (hexa.edges[1].length)^2 * 3sqrt(3) / 2
    end
    throw(error("not implement"))
end


#display函数在中间加::MIME"text/plain"
function Base.show(io::IO, hexa::AbstractHexagon{N}) where N
    name = Dict(
        :EQ => "EqHexagon",
    )
    @printf(io, "%s\n center: %s\n edgelen: %.6f",
    name[N], hexa.center, hexa.edges[1].length)
end
