"""
将一个图形进行三角切分
"""
module Triangulated

using Printf
using ..Basics

export split_square, split_hexagon
export save_triangulated, load_triangulated




"""
将一个正方形切分成4个直角三角形
"""
function split_single_square(squ::Basics.AbstractRectangle{:AXISSQUARE})
    rtver = squ.center
    top = RtTriangle(rtver, squ.vertex[4], squ.vertex[1])
    left = RtTriangle(rtver, squ.vertex[1], squ.vertex[2])
    btm = RtTriangle(rtver, squ.vertex[2], squ.vertex[3])
    right = RtTriangle(rtver, squ.vertex[3], squ.vertex[4])
    return [top, left, btm, right]
end


const RtTriOrMissing = Union{Basics.AbstractTriangle{:RT}, Missing}
const EqTriOrMissing = Union{Basics.AbstractTriangle{:EQ}, Missing}

"""
将一个大正方形切分成很多的直角三角形
"""
function split_square(squ::Basics.AbstractRectangle{:AXISSQUARE}, nps::Int64)
    lsqus = Matrix{Basics.AbstractRectangle{:AXISSQUARE}}(undef, (nps, nps))
    width = squ.width / nps
    statx, staty = squ.vertex[2].x, squ.vertex[2].y
    #从左下角开始，先水平的顺序进行标号
    # [idx2]
    # nps
    # ...
    # 1   x  x
    # 0   x  x
    #     0  1  ...  nps [idx1]
    for idx12 in CartesianIndices(lsqus)
        idx1, idx2 = Tuple(idx12)
        xct = (idx1 + 0.5) * width + statx
        yct = (idx2 + 0.5) * width + staty
        lsqus[idx1, idx2] = Square(Point2D(xct, yct), width)
    end
    #将每个小正方切分
    ltris = Array{Basics.AbstractTriangle{:RT}, 3}(undef, (nps, nps, 4))
    for idx12 in CartesianIndices(lsqus)
        idx1, idx2 = Tuple(idx12)
        ltris[idx1, idx2, :] = split_single_square(lsqus[idx1, idx2])
    end
    #找到每个小三角形挨着的小三角形
    #这个时候的顺序应该是和小三角形的边的顺序一致的
    ladjs = Array{Tuple{RtTriOrMissing, RtTriOrMissing, RtTriOrMissing}, 3}(
        undef, (nps, nps, 4)
    )
    for idx12 in CartesianIndices(lsqus)
        idx1, idx2 = Tuple(idx12)
        #上面的那个小三角，第一条边是右侧，第二条是上边
        adjs2 = idx2 < nps ? ltris[idx1, idx2+1, 3] : missing
        adjs = (ltris[idx1, idx2, 4], adjs2, ltris[idx1, idx2, 2])
        ladjs[idx1, idx2, 1] = adjs
        #左面的小三角，第一条边是上边，第二条是左面
        adjs2 = idx1 > 1 ? ltris[idx1-1, idx2, 4] : missing
        adjs = (ltris[idx1, idx2, 1], adjs2, ltris[idx1, idx2, 3])
        ladjs[idx1, idx2, 2] = adjs
        #下面的小三角，第一条边是左面，第二条是下面
        adjs2 = idx2 > 1 ? ltris[idx1, idx2-1, 1] : missing
        adjs = (ltris[idx1, idx2, 2], adjs2, ltris[idx1, idx2, 4])
        ladjs[idx1, idx2, 3] = adjs
        #右面的小三角，第一条边是下面，第二条是右边
        adjs2 = idx1 < nps ? ltris[idx1+1, idx2, 2] : missing
        adjs = (ltris[idx1, idx2, 3], adjs2, ltris[idx1, idx2, 1])
        ladjs[idx1, idx2, 4] = adjs
    end
    #将小三角和它对应的近邻都reshape
    ltris = reshape(ltris, nps*nps*4)
    ladjs = reshape(ladjs, nps*nps*4)
    return ltris, ladjs
end


"""
切分一个水平的等腰梯形，按照底边的差进行切分，所以角度必须是60度
的等腰梯形
"""
function split_trapezoid(lseg::Segment, sseg::Segment)
    dlen = lseg.length - sseg.length
    part = sseg.length / dlen
    @assert isapprox(part, round(part); atol=1e-6)
    part = Int64(round(part))
    #两个边上平均分得的所有点
    lpt = Vector{Point2D}(undef, part+2)
    spt = Vector{Point2D}(undef, part+1)
    #
    spt[1] = sseg.pt1
    lpt[1] = lseg.pt1
    #
    for idx in 1:1:part
        coef = idx / part
        coefb = 1 - coef
        spt[idx+1] = middle_point(sseg.pt1, sseg.pt2, sc1=coefb, sc2=coef)
        coef = idx / (part+1)
        coefb = 1 - coef
        lpt[idx+1] = middle_point(lseg.pt1, lseg.pt2, sc1=coefb, sc2=coef)
    end
    lpt[end] = lseg.pt2
    #
    #从长的开始计，第二个点到倒数第二个点
    eqtris = Vector{Basics.AbstractTriangle{:EQ}}(undef, 1+2*part)
    eqtris[1] = EqTriangle(lpt[1], spt[1], lpt[2])
    for idx in 1:1:part
        eqtris[2*idx] = EqTriangle(lpt[idx+1], spt[idx+1], spt[idx])
        eqtris[2*idx+1] = EqTriangle(lpt[idx+1], spt[idx+1], lpt[idx+2])
    end
    return eqtris
end


"""
切分一个水平的六边型
"""
function split_hexagon(hexa::Basics.AbstractHexagon{:EQ}, nps::Int64)
    #六边形分成两个等腰梯形
    topedge = hexa.edges[2]
    midedge = Segment(hexa.vertex[1], hexa.vertex[4])
    btmedge = hexa.edges[5]
    #
    #所有的小三角形
    eqtris = Vector{Basics.AbstractTriangle{:EQ}}(undef, 6*nps^2)
    #
    #先处理上面半个梯形
    #注意这里由于Hexagon是逆时针对顶点进行的排序
    #所以梯形和三角形增长的方向是从右往左
    llines = Vector{Segment}(undef, nps+1)
    llines[1] = topedge
    for idx in 1:1:nps
        coef = idx / nps
        coefb = 1 - coef
        llines[1+idx] = Segment(
            middle_point(topedge.pt1, midedge.pt1, sc1=coefb, sc2=coef),
            middle_point(topedge.pt2, midedge.pt2, sc1=coefb, sc2=coef)
        )
    end
    #llines中一共有nps+1个线
    #注意由于梯形的点是逆时针的，所以这边三角形的增长也是从右向左的
    sptr = 1::Int64
    for idx in 1:1:nps
        ets = split_trapezoid(llines[idx+1], llines[idx])
        @assert length(ets) == 2*nps + 1 + 2*(idx-1)
        eqtris[sptr:sptr+length(ets)-1] = ets
        sptr += length(ets)
    end
    @assert sptr == (length(eqtris) // 2)+1
    #再处理下面半个梯形的
    #这里需要注意下半个梯形的点的顺序问题
    llines = Vector{Segment}(undef, nps+1)
    llines[1] = midedge
    for idx in 1:1:nps
        coef = idx / nps
        coefb = 1 - coef
        llines[1+idx] = Segment(
            middle_point(midedge.pt1, btmedge.pt2, sc1=coefb, sc2=coef),
            middle_point(midedge.pt2, btmedge.pt1, sc1=coefb, sc2=coef)
        )
    end
    #llines中共有nps+1条线
    for idx in 1:1:nps
        ets = split_trapezoid(llines[idx], llines[idx+1])
        @assert length(ets) == 4*nps - 1 - 2*(idx-1)
        eqtris[sptr:sptr+length(ets)-1] = ets
        sptr += length(ets)
    end
    @assert sptr == length(eqtris)+1
    #开始处理相邻关系
    ladjs = Vector{Tuple{EqTriOrMissing, EqTriOrMissing, EqTriOrMissing}}(
        undef, length(eqtris)
    )
    # Warning: 这个ladjs中的相邻的三角形的顺序一定和Eqtriangle的边的顺序是一样的
    #先是nps个上半层，其中三角形的数量是2*nps+1+2*idx
    #idx从0到nps-1
    #先处理最上面的一层的相邻
    #最上面一层就是从[1,2*nps+1]这个范围内的所有
    for idx in 1:1:(2*nps+1)
        #每个三角都有三条边，于是有三个相邻
        #这些是底边在下面的相邻的有下一个，和下一行正对的
        #上半部分中，底边在下面的三角形的点是按照右-左-下的顺序来的
        #底边在上面的三角形是左-上-右
        if idx % 2 == 1
            right = idx == 1 ? missing : eqtris[idx-1]
            left = idx == 2*nps+1 ? missing : eqtris[idx+1]
            ladjs[idx] = (right, left, eqtris[idx+2*nps+2])
        else
            ladjs[idx] = (eqtris[idx+1], missing, eqtris[idx-1])
        end
    end
    #之后，还有nps-1个上半部分的梯形
    sptr = 2*nps + 1
    for npsidx in 2:1:nps
        trinum = 2*nps+1+2*(npsidx-1)
        #如果到了最后一个，下面的一行的三角形的数量和这一行是一样多的
        #计算的时候不需要多添加一个1
        offset = npsidx == nps ? 0 : 1
        for idx in 1:1:trinum
            tot_idx = sptr + idx
            #上半部分中，底边在下面的三角形的点是按照右-左-下的顺序来的
            #底边在上面的三角形是左-上-右
            if idx % 2 == 1
                right = idx == 1 ? missing : eqtris[tot_idx-1]
                left = idx == trinum ? missing : eqtris[tot_idx+1]
                ladjs[tot_idx] = (right, left, eqtris[tot_idx+trinum+offset])
            else
                ladjs[tot_idx] =
                (eqtris[tot_idx+1], eqtris[tot_idx-trinum+1], eqtris[tot_idx-1])
            end
        end
        sptr += trinum
    end
    #下半部分从后往前开始处理
    #最下面的一个梯形
    for idx in 1:1:2*nps+1
        #每个三角都有三条边，于是有三个相邻
        #这些是底边在上面的，相邻的有下一个，上面对的，还有上一个
        #下半部分中，底边在上面的三角形是按照右-左-上的顺序来的
        #底边在下面的三角形是按照左-下-右的顺序来的
        tot_idx = length(eqtris) - (2*nps+1) + idx
        if idx % 2 == 1
            rightone = idx == 1 ? missing : eqtris[tot_idx-1]
            leftone = idx == 2*nps+1 ? missing : eqtris[tot_idx+1]
            ladjs[tot_idx] = (rightone, leftone, eqtris[tot_idx-2-2*nps])
        else
            ladjs[tot_idx] = (eqtris[tot_idx+1], missing, eqtris[tot_idx-1])
        end
    end
    sptr = length(eqtris) - (2*nps+1)
    for npsidx in 2:1:nps
        offset = npsidx != nps ? 1 : 0
        trinum = 2*nps+1+2*(npsidx-1)
        for idx in 1:1:trinum
            tot_idx = sptr - trinum + idx
            #下半部分中，底边在上面的三角形是按照右-左-上的顺序来的
            if idx % 2 == 1
                rightone = idx == 1 ? missing : eqtris[tot_idx-1]
                leftone = idx == trinum ? missing : eqtris[tot_idx+1]
                ladjs[tot_idx] = (rightone, leftone, eqtris[tot_idx-offset-trinum])
            #下半部分中，底边在下面的三角形是按照左-下-右的顺序来的
            else
                ladjs[tot_idx] =
                (eqtris[tot_idx+1], eqtris[tot_idx+trinum-1], eqtris[tot_idx-1])
            end
        end
        sptr -= trinum
    end
    return eqtris, ladjs
end



"""
平行于坐标轴的图形的字符串
"""
macro axis_polygon_string(poly)
    type2name = Dict(
        Basics.AbstractRectangle{:AXISSQUARE} => "Square",
        Basics.AbstractHexagon{:EQ} => "EqHexagon"
    )
    return quote
        #宏里面不要加return，代换回去以后直接就把函数返回了
        par1 = @sprintf "%s\n" $(esc(type2name))[typeof($(esc(poly)))]
        par2 = @sprintf "%.12f,%.12f\n" $(esc(poly)).center.x $(esc(poly)).center.y
        par3 = @sprintf "%.12f\n" $(esc(poly)).edges[1].length
        par1 * par2 * par3
    end
end 


"""
三角型的字符串
"""
macro triangle_string(tri)
    type2name = Dict(
        Basics.AbstractTriangle{:COM} => "Triangle",
        Basics.AbstractTriangle{:EQ} => "EqTriangle",
        Basics.AbstractTriangle{:RT} => "RtTriangle"
    )
    return esc(quote
        par1 = @sprintf "%s:" $(type2name)[typeof($tri)]
        pt1 = @sprintf "%.12f,%.12f:" ($tri).vertex[1].x $(tri).vertex[1].y
        pt2 = @sprintf "%.12f,%.12f:" ($tri).vertex[2].x ($tri).vertex[2].y
        pt3 = @sprintf "%.12f,%.12f\n" ($tri).vertex[3].x ($tri).vertex[3].y
        par1*pt1*pt2*pt3
    end)
end

"""
保存结果到文件, 会分别保存到三个文件里面
"""
function save_triangulated(
    poly::AbstractPolygon,
    ltris::Vector{T},
    ladjs::Vector{Tuple{Union{Missing, T}, Union{Missing, T}, Union{Missing, T}}},
    prefix
) where T <: Basics.AbstractTriangle
    @assert isa(poly, Union{
        Basics.AbstractRectangle{:AXISSQUARE},
        Basics.AbstractHexagon{:EQ}}
        )
    #保存描述图形的文件
    plystr = @axis_polygon_string poly
    outf = open(prefix*".ply", "w")
    write(outf, plystr)
    close(outf)
    #保存描述三角形的文件
    outf = open(prefix*".tri", "w")
    tri2idx = Dict{
        Union{Missing, AbstractPolygon},
        Int64
    }(missing => -1)
    for (idx, tri) in enumerate(ltris)
        tristr = @triangle_string tri
        tri2idx[tri] = idx
        write(outf, tristr)
    end
    close(outf)
    #保存描述相邻三角型的文件
    outf2 = open(prefix*".adj", "w")
    for idx in 1:1:length(ladjs)
        adjs = ladjs[idx]
        aidx1, aidx2, aidx3 = tri2idx[adjs[1]], tri2idx[adjs[2]], tri2idx[adjs[3]]
        write(outf2, @sprintf "%d:%d:%d\n" aidx1 aidx2 aidx3)
    end
    close(outf2)
end



"""
读取配置
"""
function load_triangulated(prefix)
    #从名称到类型
    name2type = Dict(
        "Triangle" => Basics.AbstractTriangle{:COM},
        "EqTriangle" => Basics.AbstractTriangle{:EQ},
        "RtTriangle" => Basics.AbstractTriangle{:RT}
    )
    #读取多边形
    inf = open(prefix*".ply", "r")
    plystr = readlines(inf)
    cntstr = split(plystr[2], ',')
    poly = eval(Symbol(plystr[1]))(
        Point2D(parse(Float64, cntstr[1]), parse(Float64, cntstr[2])),
        parse(Float64, plystr[3])
    )
    close(inf)
    #读取切分的三角形
    inf = open(prefix*".tri", "r")
    ltristr = readlines(inf)
    typestr = split(ltristr[1], ':')[1]
    ltris = Vector{name2type[typestr]}(undef, length(ltristr))
    for (idx, lstr) in enumerate(ltristr)
        context = split(lstr, ':')
        pt1 = split(context[2], ',')
        pt2 = split(context[3], ',')
        pt3 = split(context[4], ',')
        ltris[idx] = eval(Symbol(context[1]))(
            Point2D(parse(Float64, pt1[1]), parse(Float64, pt1[2])),
            Point2D(parse(Float64, pt2[1]), parse(Float64, pt2[2])),
            Point2D(parse(Float64, pt3[1]), parse(Float64, pt3[2]))
        )
    end
    close(inf)
    #读取近邻的配置
    inf = open(prefix*".adj", "r")
    ladjstr = readlines(inf)
    typeormis = Union{Missing, name2type[typestr]}
    ladjs = Vector{Tuple{typeormis, typeormis, typeormis}}(
        undef, length(ladjstr)
    )
    for (idx, lstr) in enumerate(ladjstr)
        context = split(lstr, ':')
        adj1 = context[1] == "-1" ? missing : ltris[parse(Int64, context[1])]
        adj2 = context[2] == "-1" ? missing : ltris[parse(Int64, context[2])]
        adj3 = context[3] == "-1" ? missing : ltris[parse(Int64, context[3])]
        ladjs[idx] = (adj1, adj2, adj3)
    end
    close(inf)
    return poly, ltris, ladjs
end


end # end module

