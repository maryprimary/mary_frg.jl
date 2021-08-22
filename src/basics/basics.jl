"""
包含最基础功能的模块
"""
module Basics
    include("point.jl")
    export Point2D
    export middle_point, absolute_angle

    include("line.jl")
    export Segment

    include("polygon.jl")
    export AbstractPolygon

    include("triangle.jl")
    export Triangle, area
    export RtTriangle, EqTriangle

    include("rectangle.jl")
    export Rectangle, Square, area

    include("hexagon.jl")
    export EqHexagon, area


end