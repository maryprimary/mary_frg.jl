"""
Bubble积分
qpp = k1+k2; qfs = k3-k2; qex = k1-k3
"""
struct Bubble{T, P}
    lval :: Float64
    #这里面的顺序是α, β, b_1, b_2, k_n, k_1, k_2
    #k1，k2, b1, b2会用来决定qpp,qfs等数值
    V :: Array{Float64, 7}
end




