using Test
using MARY_fRG.Basics


@testset "长方形" begin
    rect1 = Rectangle(Point2D(0., 0.), 1, 2)
    println(rect1)
    @test area(rect1) == 2
end


@testset "正方形" begin
    squ1 = Square(Point2D(0., 0.), 2)
    println(squ1)
    @test area(squ1) == 4
end
