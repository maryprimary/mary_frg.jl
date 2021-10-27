"""
测试包围费米面
"""


using Test
using MARY_fRG.Fermi
using MARY_fRG.FlowEquation
using MARY_fRG.Drawers
using MARY_fRG.Fermi.Surface


function squ()
    qs = common_square_lattice(0.20)
    Γ4 = TFGamma4(qs, 8.0, 16, 20)
    plt = draw_polygon(Γ4.ltris)
    @savepng plt "refine3_1"
    #
    Γ4 = refine_to_surface(Γ4)
    plt = draw_polygon(Γ4.ltris)
    @savepng plt "refine3_2"
    #
    Γ4 = refine_to_surface(Γ4)
    plt = draw_polygon(Γ4.ltris)
    @savepng plt "refine3_3"
end

squ()
