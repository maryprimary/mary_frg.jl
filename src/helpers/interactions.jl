"""
一些常见的相互作用
"""

module Interactions

using ..Basics


export triangular_system_heisenberg


"""
三角格子上交换相互作用的Γ4
"""
function triangular_system_heisenberg(patcs::Vector{Point2D}, kadd, J)
    pnum = length(patcs)
    result = zeros(pnum, pnum, pnum)
    co1 = 0.5
    co2 = sqrt(3) * 0.5
    for idxs in CartesianIndices(result)
        i1, i2, i3 = Tuple(idxs)
        k1, k2, k3 = patcs[i1], patcs[i2], patcs[i3]
        q1 = kadd(k3, -k1)
        q2 = kadd(k3, -k2)
        phase = cos(-co1*q1.x + co2*q1.y)
        phase += cos(co1*q1.x + co2*q1.y)
        phase += cos(q1.x)
        phase += 0.5*cos(-co1*q2.x + co2*q2.y)
        phase += 0.5*cos(co1*q2.x + co2*q2.y)
        phase += 0.5*cos(q2.x)
        result[i1, i2, i3] = -J * phase
    end
    return result
end

    
end
