"""
和patch相关的功能
"""
module Patch

using ..Fermi
using ..Basics


export find_patch_index_hexa, find_patch_index_squa
export group_ltris_into_patches_mt
export patches_under_vonhove
export get_k4_table


"""
在六角型的区域里找到属于哪个patch
"""
function find_patch_index_hexa(pt::Point2D, hexa::Basics.AbstractHexagon{:EQ}, pnum::Int64)
    if abs(pt.x) < 1e-8 && abs(pt.y) < 1e-8
        return missing
    end
    #把六边型分成六瓣
    pre_pnum = Int64(pnum // 6)
    #找到属于哪个以后，先转到水平的左侧那个上面
    angle = absolute_angle(pt)
    #println(angle)
    partidx = (3 * angle / pi) + 0.5
    partidx = Int64(floor(partidx)) % 6
    #println(partidx)
    rotvec1 = [cos((partidx-3)*pi/3), sin((partidx-3)*pi/3)]
    rotvec2 = [-sin((partidx-3)*pi/3), cos((partidx-3)*pi/3)]
    #trieng = -2cos(pt.x) - 4cos(sqrt(3)*pt.y/2.0)*cos(pt.x/2.0)
    #将点旋转
    rpt = Point2D(
        pt.x * rotvec1[1] + pt.y * rotvec1[2],
        pt.x * rotvec2[1] + pt.y * rotvec2[2]
    )
    #计算属于哪个patch
    trieng = -2cos(rpt.x) - 4cos(sqrt(3)*rpt.y/2.0)*cos(rpt.x/2.0) - 2.
    #利用三角格子的von Hove filling来判断在里面还是外面
    #如果是在外面
    if trieng > 0
        #EqHexagon的右上角就是第四个角
        #为了和内侧有相同的角度，这里水平伸长三倍，然后放到另一侧计算角度
        spt = rpt - hexa.vertex[4]
        spt = Point2D(-3spt.x, -spt.y)
        relang = (7pi / 6) - absolute_angle(spt)
        relidx = relang * pre_pnum / (pi/3)
        #由于精度原因，relang有可能小于0的极小数，floor就变成-1了
        relidx = Int64(floor(relidx)) + 1
    else
        relang = absolute_angle(rpt) - 5pi / 6
        relidx = relang * pre_pnum / (pi/3)
        relidx = Int64(floor(relidx)) + 1
        #if relidx > 1
        #    println(relang)
        #    println(partidx)
        #end
    end
    #由于一个patch边界上的可能会判断为下一个patch（floor以后刚好是pre_pnum)
    #所以是有可能大于pnum数量的
    #0的话代表relidx是-1,是上一个patch，转回去就是最后一个
    pidx = (relidx + pre_pnum * partidx) % pnum
    return pidx == 0 ? pnum : pidx
end



"""
在正方形的区域中找到属于哪个patch
"""
function find_patch_index_squa(pt::Point2D,
    squa::Basics.AbstractRectangle{:AXISSQUARE}, pnum::Int64)
    if abs(pt.x) < 1e-8 && abs(pt.y) < 1e-8
        return missing
    end
    #把正方形分成4份
    pre_pnum = Int64(pnum // 4)
    #先转到左上角的那个里面
    #右上角的话，在x轴有个角度的突变
    angle = absolute_angle(pt)
    partidx = 2 * angle / pi
    partidx = Int64(floor(partidx)) % 4
    #
    rotvec1 = [cos((partidx-1)*pi/2), sin((partidx-1)*pi/2)]
    rotvec2 = [-sin((partidx-1)*pi/2), cos((partidx-1)*pi/2)]
    #将点旋转
    rpt = Point2D(
        pt.x * rotvec1[1] + pt.y * rotvec1[2],
        pt.x * rotvec2[1] + pt.y * rotvec2[2]
    )
    #根据正方格子的von Hove filling判断里外
    squeng = -2*(cos(rpt.x) + cos(rpt.y))
    #如果是在里面
    if squeng > 0
        #Square的左上角就是第一个角
        relang = (2pi) - absolute_angle(rpt - squa.vertex[1])
        relidx = relang * pre_pnum / (pi / 2)
        #由于精度原因，relang有可能小于0的极小数，floor就变成-1了
        relidx = Int64(floor(relidx)) + 1
    else
        relang = absolute_angle(rpt) - (pi / 2)
        relidx = relang * pre_pnum / (pi / 2)
        relidx = Int64(floor(relidx)) + 1
    end
    pidx = (relidx + pre_pnum * partidx) % pnum
    return pidx == 0 ? pnum : pidx
end


"""
把所有的小三角型分好组
"""
function group_ltris_into_patches_mt(
    ltris::Vector{T}, brlu::P, pnum
) where {T <: Basics.AbstractTriangle, P <: AbstractPolygon}
    find_algo = Dict(
        Basics.AbstractRectangle{:AXISSQUARE} => find_patch_index_squa,
        Basics.AbstractHexagon{:EQ} => find_patch_index_hexa
    )[P]
    #
    lpats = Vector{Int64}(undef, length(ltris))
    Threads.@threads for idx in 1:1:length(ltris) # (idx, tri) in enumerate(ltris)
        tri = ltris[idx]
        #lpats[idx] = find_algo(tri.center, brlu, pnum)
        tpat = find_algo(tri.center, brlu, pnum)
        if ismissing(tpat)
            tpat = 0
        end
        lpats[idx] = tpat
    end
    return lpats
end


"""
二分法求0点
"""
function bisect_fermi_suface(angle, model, bandidx, rad1, rad2)
    cosa = cos(angle)
    sina = sin(angle)
    lv = dispersion(model, bandidx, rad1 * cosa, rad1 * sina)
    rv = dispersion(model, bandidx, rad2 * cosa, rad2 * sina)
    @assert lv * rv <= 0
    if lv < 0
        left, right = rad1, rad2
    else
        left, right = rad2, rad1
    end
    middle = (left + right) / 2
    mv = dispersion(model, bandidx, middle * cos(angle), middle * sin(angle))
    count = 0
    while !isapprox(mv, 0., atol=1e-12)
        if mv < 0
            left = middle
        else
            right = middle
        end
        middle = (left + right) / 2
        mv = dispersion(model, bandidx, middle * cosa, middle * sina)
        count += 1
        if count > 50
            throw(error("can not find root"))
        end
    end
    return middle
end


"""
获得向里面的费米面上的patches，所以不需要关心具体的布里渊区形状
需要自己确认费米面在von Hove fill 以下
"""
function patches_under_vonhove(
    model::T, bandidx, pnum
    ) where T <: Fermi.Abstract2DModel
    dang = 2pi / pnum
    radius = isa(model.brillouin, Basics.AbstractHexagon{:EQ}) ? 2pi*sqrt(3)/3 : pi
    sang = isa(model.brillouin, Basics.AbstractHexagon{:EQ}) ? -pi / 6 : 0
    patches = Vector{Point2D}(undef, pnum)
    for idx = 1:1:pnum
        #稍微偏离一些防止正好卡在patch边界的k4
        ang = (idx - π/6) * dang + sang
        root = bisect_fermi_suface(ang, model, bandidx, 0., radius)
        patches[idx] = Point2D(root*cos(ang), root*sin(ang))
    end
    return patches
end

    

#"""
#根据动量守恒得到k4所属的patches
#"""
#function get_k4_table(sys::T, patcs::Vector{Point2D}
#    )where T <: Union{QuadrateSystem, TriangularSystem}
#    println(T)
#    find_algo = T <: QuadrateSystem ? find_patch_index_squa :
#    find_patch_index_hexa
#    pnum = length(patcs)
#    k4tab = Array{Int64, 3}(undef, pnum, pnum, pnum)
#    Threads.@threads for idxs in CartesianIndices(k4tab)
#        i1, i2, i3 = Tuple(idxs)
#        k4p = sys.kadd(patcs[i1], patcs[i2])
#        k4p = sys.kadd(k4p, -patcs[i3])
#        k4tab[i1, i2, i3] = find_algo(k4p, sys.brillouin, pnum)
#    end
#    return k4tab
#end


end # end module

