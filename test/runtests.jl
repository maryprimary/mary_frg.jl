using Test

using MARY_fRG


include("test_basics/test_points.jl")

include("test_basics/test_lines.jl")

include("test_basics/test_triangles.jl")

include("test_basics/test_rectangle.jl")

include("test_basics/test_hexagon.jl")

#绘图的测试
#include("test_drawers/test_polygon.jl")

include("test_triangulated/test_square.jl")

include("test_triangulated/test_hexagon.jl")

#refine
#include("test_triangulated/test_refine.jl")
