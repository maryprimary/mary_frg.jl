"""
和费米面相关的功能
"""
module Surface
    

using ..Fermi
using ..Basics


export const_energy_line, const_energy_line_in_patches

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


end # end module



