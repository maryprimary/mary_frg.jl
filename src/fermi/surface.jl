"""
和费米面相关的功能
"""
module Surface
    

using ..Fermi
using ..Basics


export @onsurface
export const_energy_line, const_energy_line_in_patches
export const_energy_triangle


"""
获取等能面
"""
function const_energy_line(
    ltris::Vector{T},
    ladjs::Vector{Tuple{Union{Missing, T}, Union{Missing, T}, Union{Missing, T}}},
    eng::Float64,
    disp::Function) where T <: Basics.AbstractTriangle 
    #计算每个位置上面的能量
    eng_dict = Dict()
    for tri in ltris
        #能量这里减去需要寻找的能量大小
        eng_dict[tri] = disp(tri.center.x, tri.center.y) - eng
    end
    #
    edges::Vector{Basics.Segment} = []
    for (idx, tri) in enumerate(ltris)
        #寻找从负数到正数的边界
        #跳过所有的正数
        if eng_dict[tri] > 0.
            continue
        end
        #查看相邻的几个有没有大于零的
        for (eidx, adj) in enumerate(ladjs[idx])
            #已经到头的时候是不继续添加的
            if ismissing(adj)
                continue
            end
            if eng_dict[adj] > 0
                push!(edges, tri.edges[eidx])
            end
        end
    end
    return edges
end



"""
获取费米面，但是要分好属于哪一组
"""
function const_energy_line_in_patches(
    ltris::Vector{T},
    ladjs::Vector{Tuple{Union{Missing, T}, Union{Missing, T}, Union{Missing, T}}},
    lpats::Vector{Int64},
    eng::Float64,
    disp::Function) where T <: Basics.AbstractTriangle 
    eng_dict = Dict()
    for tri in ltris
        #能量这里减去需要寻找的能量大小
        eng_dict[tri] = disp(tri.center.x, tri.center.y) - eng
    end
    #
    edges::Vector{Basics.Segment} = []
    epidx::Vector{Int64} = []
    for (idx, tri) in enumerate(ltris)
        #寻找从负数到正数的边界
        if eng_dict[tri] > 0.
            continue
        end
        #查看相邻的几个有没有大于零的
        for (eidx, adj) in enumerate(ladjs[idx])
            #已经到头的时候是不继续添加的
            if ismissing(adj)
                continue
            end
            if eng_dict[adj] > 0
                push!(edges, tri.edges[eidx])
                #如果大于0和小于0的不在一个patch，会产生一些混淆。但是没明显的影响。
                push!(epidx, lpats[idx])
            end
        end
    end
    return edges, epidx
end


"""
判断小三角是否穿过费米面
"""
macro onsurface(tri, disp, eng)
    return esc(quote
        ver1 = ($tri).vertex[1]
        sgn1 = sign($disp(ver1.x, ver1.y)-$eng)
        ver2 = ($tri).vertex[2]
        sgn2 = sign($disp(ver2.x, ver2.y)-$eng)
        ver3 = ($tri).vertex[3]
        sgn3 = sign($disp(ver3.x, ver3.y)-$eng)
        sgn1 != sgn2 || sgn1 != sgn3
    end)
end



"""
穿过费米面的小三角
"""
function const_energy_triangle(
    ltris::Vector{T},
    eng::Float64,
    disp::Function) where T <: Basics.AbstractTriangle
    #
    #如果一个小三角型的三个顶点的能量符号不同，那就是穿过了费米面
    #
    edges::Vector{T} = []
    for tri in ltris
        if (@onsurface tri disp eng)
            push!(edges, tri)
        end
    end
    return edges
end



end # end module



