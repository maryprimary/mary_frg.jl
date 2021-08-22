"""
绘制相关的功能
"""
module Drawers
    
using Plots
using ..Basics

export @savepng
export draw_lines, draw_lines!
export draw_polygon, draw_polygon!
export draw_points, draw_points!
export draw_figure

"""
绘制散点图
"""
function draw_points(pts::Array{Point2D, 1}; kw...)
    fig = plot()
    draw_points!(fig, pts; kw...)
    return fig
end

"""
绘制散点图
"""
function draw_points!(plt, pts::Array{Point2D, 1}; kw...)
    scatter!(plt, [pt.x for pt in pts], [pt.y for pt in pts];
    legend=false, kw...)
    return plt
end


"""绘制线"""
function draw_lines(lns::Array{Segment, 1}; lcidx=nothing, kw...)
    fig = plot()
    draw_lines!(fig, lns; lcidx=lcidx, kw...)
    return fig
end


"""绘制线"""
function draw_lines!(plt, lns::Array{Segment, 1}; lcidx=nothing, kw...)
    plot!(plt, legend=false, kw...)
    colors = [:red, :orange, :yellow, :green, :grays, :blue, :purple]
    custom_color = !isnothing(lcidx)
    for (idx, seg) in enumerate(lns)
        cidx = custom_color ? lcidx[idx] : idx
        plot!(plt, [seg.pt1.x, seg.pt2.x], [seg.pt1.y, seg.pt2.y],
        color=colors[(cidx-1)%7+1])
    end
end


"""绘制多边形"""
function draw_polygon(pls::Array{<:AbstractPolygon, 1}; pcidx=nothing, kw...)
    fig = plot()
    draw_polygon!(fig, pls; pcidx=pcidx, kw...)
    return fig
end


"""绘制多边形"""
function draw_polygon!(plt, pls::Array{<:AbstractPolygon, 1}; pcidx=nothing, kw...)
    plot!(plt, legend=false, kw...)
    colors = [:red, :orange, :yellow, :green, :grays, :blue, :purple]
    custom_color = !isnothing(pcidx)
    for (idx, ply) in enumerate(pls)
        cidx = custom_color ? pcidx[idx] : idx
        ptxs = [pt.x for pt in ply.vertex]
        push!(ptxs, ply.vertex[1].x)
        ptys = [pt.y for pt in ply.vertex]
        push!(ptys, ply.vertex[1].y)
        plot!(plt,
            ptxs, ptys,
            lw=0,
            fill=(0, colors[(cidx-1)%7+1]),
            fillalpha=0.5
        )
    end
end


"""
绘制各个成分都有的图形
"""
function draw_figure(
    pts::Union{Nothing, Array{Point2D, 1}},
    lns::Union{Nothing, Array{Segment, 1}},
    pls::Union{Nothing, Array{<:AbstractPolygon, 1}};
    ptext=nothing, lcidx=nothing, pcidx=nothing
    )
    plt = plot()
    isnothing(pts) ? nothing : draw_points!(plt, pts, text=ptext)
    isnothing(lns) ? nothing : draw_lines!(plt, lns, lcidx=lcidx)
    isnothing(pls) ? nothing : draw_polygon!(plt, pls, pcidx=pcidx)
    return plt
end


"""
保存图片
"""
macro savepng(plt, name)
    esc(quote
        $png($plt, $name)
    end)
end


end
