#=
测试对称性
=#

using MARY_fRG.Basics
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation


function run()
    k1c = Point2D(1., 1.)
    k1tri = RtTriangle(k1c, k1c + Point2D(0.01, 0.), k1c+Point2D(1., 1.01))
    qfs = Point2D(-1.55, 2.401)
    model = wuxx_av3sb5_kagome_lattice()
    println(FlowEquation.pi_αβ_plus_tf(model, [k1tri], 1., 1., qfs, 1, 1))
    k2c = k1c - qfs
    k2tri = RtTriangle(k2c, k2c + Point2D(0.01, 0.), k2c+Point2D(1., 1.01))
    qfs = -qfs
    println(FlowEquation.pi_αβ_plus_tf(model, [k2tri], 1., 1., qfs, 1, 1))
end


run()

