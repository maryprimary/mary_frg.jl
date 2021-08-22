"""
三角形有关的功能
"""

using Printf


struct AbstractTriangle{T} <: AbstractPolygon{T}
    vertex :: Array{Point2D, 1}
    edges :: Array{Segment, 1}
    center :: Point2D
end


"""
构造一个普通的三角形
"""
function Triangle(pt1, pt2, pt3)
    verts = [pt1, pt2, pt3]
    edges = Array{Segment, 1}(undef, 3)
    edges[1] = Segment(pt1, pt2)
    edges[2] = Segment(pt2, pt3)
    edges[3] = Segment(pt3, pt1)
    center = Point2D(
        sum([pt1.x, pt2.x, pt3.x]) / 3.0,
        sum([ver.y for ver in verts]) / 3.0
    )
    return AbstractTriangle{:COM}(verts, edges, center)
end

"""
构造一个普通的三角形
"""
function Triangle(verts::Array{Point2D, 1})
    @assert length(verts) == 3
    return Triangle(verts[1], verts[2], verts[3])
end



"""
构造一个直角三角形，第一个角必须是直角
"""
function RtTriangle(pt1, pt2, pt3)
    verts = [pt1, pt2, pt3]
    edges = Array{Segment, 1}(undef, 3)
    edges[1] = Segment(pt1, pt2)
    edges[2] = Segment(pt2, pt3)
    edges[3] = Segment(pt3, pt1)
    center = middle_point(
        pt1,
        middle_point(pt2, pt3),
        sc1=1., sc2=1.414
    )
    return AbstractTriangle{:RT}(verts, edges, center)
end

"""
构造一个直角三角形，会检查第一个角是否是直角
"""
function RtTriangle(verts::Array{Point2D, 1})
    @assert length(verts) == 3
    vec1 = verts[2] - verts[1]
    vec2 = verts[3] - verts[1]
    @assert abs(vec1.x * vec2.x + vec1.y * vec2.y) < 1e-6
    return RtTriangle(verts[1], verts[2], verts[3])
end


"""
构造一个等边三角形
"""
function EqTriangle(pt1, pt2, pt3)
    verts = [pt1, pt2, pt3]
    edges = Array{Segment, 1}(undef, 3)
    edges[1] = Segment(pt1, pt2)
    edges[2] = Segment(pt2, pt3)
    edges[3] = Segment(pt3, pt1)
    center = Point2D(
        sum([pt1.x, pt2.x, pt3.x]) / 3.0,
        sum([ver.y for ver in verts]) / 3.0
    )
    return AbstractTriangle{:EQ}(verts, edges, center)
end


"""
构造一个等边三角形，会检查边长
"""
function EqTriangle(verts::Array{Point2D, 1})
    @assert length(verts) == 3
    ret = EqTriangle(verts[1], verts[2], verts[3])
    @assert abs(ret.edges[1].length - ret.edges[2].length) < 1e-6
    @assert abs(ret.edges[3].length - ret.edges[2].length) < 1e-6
    return ret
end

"""
返回一个三角形的面积
"""
function area(tri::AbstractTriangle{N}) where N
    vec1 = tri.vertex[2] - tri.vertex[1]
    vec2 = tri.vertex[3] - tri.vertex[1]
    return 0.5abs(vec1.x * vec2.y - vec1.y * vec2.x)
end

#display函数在中间加::MIME"text/plain"
function Base.show(io::IO, tri::AbstractTriangle{N}) where N
    name = Dict(
        :COM => "Triangle",
        :RT => "RTTriangle",
        :EQ => "EqTriangle"
    )
    @printf(io, "%s\n center: %s\n vertex: %s",
    name[N], tri.center, tri.vertex)
end
