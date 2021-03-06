"""
绘制bubble贡献的图
"""


using Plots

"""
bubb的贡献
"""
function bubb(x, y)
    eps_k = x
    neps_kp = y
    lamb = 1e-2
    if abs(eps_k - neps_kp) < 1.e-10
        #如果两个数值比较接近, Pi^{-}和Pi^{+}的公式完全一样，就是第二个能量要加个负号
        # lim (eps_k -> -eps_kp) Pi^{-} =
        # 1/T (e^{eps/T} (-eps/T*e^{eps/T} + eps/T + e^{eps/T} + 1)) 
        #/ (e^{eps/T} + 1)^3
        bval = eps_k / lamb
        if bval > 25
            #@info "数值不稳定"
            return 0.
        end
        expb = exp(bval)
        num = expb * (-bval * expb + bval + expb + 1)
        den = (1+expb)^3
        d_val = num / den / lamb
    else
        if (eps_k / lamb) > 25
            #@info "数值不稳定"
            num_left = 0.
        else
            #e^{epsilon_k / T}
            exp_k_t = exp(eps_k / lamb)
            num_left = eps_k / lamb * exp_k_t / ((1 + exp_k_t)^2)
        end
        if (neps_kp / lamb) > 25
            #@info "数值不稳定"
            num_righ = 0.
        else
            #e^{-epsilon_{-k+q} / T}
            exp_nkp_t = exp(neps_kp / lamb)
            num_righ = neps_kp / lamb * exp_nkp_t / 
            ((1 + exp_nkp_t)^2)
        end
        d_val = (num_left - num_righ) / (eps_k - neps_kp)
    end # end if 小区域的贡献计算完
    return d_val
end


surface(-4:0.1:4, -4:0.1:4, bubb)

savefig("bubb_contrib2.png")



