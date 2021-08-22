"""
和线有关的功能
"""

struct Segment
    pt1 :: Point2D
    pt2 :: Point2D
    length :: Float64
    Segment(pt1::Point2D, pt2::Point2D) = begin
        len = sqrt((pt1.x - pt2.x)^2 + (pt1.y - pt2.y)^2)
        new(pt1, pt2, len)
    end
end

