"""
将三角切分的结果再进行切分
"""
module Refine
    
using ..Basics

export split_triangle, refine_triangles, find_adjs_by_adjoint


"""
将三角形重新切分成四个，可以用在直角三角形上
依旧保证出来的四个是直角三角形，且第一个角是直角
"""
@generated function split_triangle(rtri::Basics.AbstractTriangle{T}) where T
    type2func = Dict(
        :COM => Triangle,
        :RT => RtTriangle,
        :EQ => EqTriangle
    )
    func = type2func[T]
    if T == :RT
        return quote
            #新的四个直角三角形的四个直角
            rtvexs = Vector{Point2D}(undef, 4)
            rtvexs[1] = rtri.vertex[1]
            rtvexs[2] = middle_point(rtri.vertex[1], rtri.vertex[2])
            rtvexs[3] = middle_point(rtri.vertex[2], rtri.vertex[3])
            rtvexs[4] = middle_point(rtri.vertex[3], rtri.vertex[1])
            #
            newtris = Vector{Basics.AbstractTriangle{T}}(undef, 4)
            newtris[1] = RtTriangle(rtvexs[2], rtvexs[3], rtvexs[1])
            newtris[2] = RtTriangle(rtvexs[2], rtri.vertex[2], rtvexs[3])
            newtris[3] = RtTriangle(rtvexs[4], rtvexs[3], rtri.vertex[3])
            newtris[4] = RtTriangle(rtvexs[4], rtvexs[3], rtvexs[1])
            return newtris
        end
    end
    return quote
        #新的四个直角三角形的四个直角
        rtvexs = Vector{Point2D}(undef, 4)
        rtvexs[1] = rtri.vertex[1]
        rtvexs[2] = middle_point(rtri.vertex[1], rtri.vertex[2])
        rtvexs[3] = middle_point(rtri.vertex[2], rtri.vertex[3])
        rtvexs[4] = middle_point(rtri.vertex[3], rtri.vertex[1])
        #
        newtris = Vector{Basics.AbstractTriangle{T}}(undef, 4)
        newtris[1] = ($func)(rtvexs[1], rtvexs[2], rtvexs[4])
        newtris[2] = ($func)(rtvexs[2], rtri.vertex[2], rtvexs[3])
        newtris[3] = ($func)(rtvexs[3], rtvexs[4], rtvexs[2])
        newtris[4] = ($func)(rtvexs[4], rtvexs[3], rtri.vertex[3])
        return newtris
    end # end quote
end


"""
将一组的三角形细化
"""
function refine_triangles(ltris::Vector{T}) where T <: Basics.AbstractTriangle
    newtris = Matrix{T}(undef, length(ltris), 4)
    for (idx, tri) in enumerate(ltris)
        newtris[idx, :] = split_triangle(tri)
    end
    newtris = reshape(newtris, length(ltris)*4)
    return newtris
end


"""
判断一条边是不是某个三角形的
"""
macro isedge(edge, tri)
    return esc(quote
        topt1 = [vex - ($edge).pt1 for vex in ($tri).vertex]
        topt1 = [sqrt(pt.x^2 + pt.y^2) for pt in topt1]
        haspt1 = isapprox(minimum(topt1), 0., atol=1e-6)
        topt2 = [vex - ($edge).pt2 for vex in ($tri).vertex]
        topt2 = [sqrt(pt.x^2 + pt.y^2) for pt in topt2]
        haspt2 = isapprox(minimum(topt2), 0., atol=1e-6)
        haspt1 && haspt2
    end)
end


"""
通过判断边是否重合来判断两个三角形是否挨在一起
"""
function find_adjs_by_adjoint(ltris::Vector{T}) where T <: Basics.AbstractTriangle
    ladjs = Vector{Tuple{
        Union{Missing, T}, Union{Missing, T}, Union{Missing, T}
    }}(undef, length(ltris))
    for idx in 1:1:length(ltris)
        ladjs[idx] = (missing, missing, missing)
    end
    @Threads.threads for idx in 1:1:length(ltris)
        tri = ltris[idx]
        fadjs = Vector{Union{Missing, T}}(undef, 3)
        fadjs[:] = [missing, missing, missing]
        for (adjidx, adjtri) in enumerate(ltris)
            if idx == adjidx
                continue
            end
            if @isedge tri.edges[1] adjtri
                fadjs[1] = adjtri
                continue
            end
            if @isedge tri.edges[2] adjtri
                fadjs[2] = adjtri
                continue
            end
            if @isedge tri.edges[3] adjtri
                fadjs[3] = adjtri
                continue
            end
        end
        ladjs[idx] = Tuple(fadjs)
    end
    return ladjs
end



end # end module

