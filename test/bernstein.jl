module TestBernstein

using Test
using QuadLearnData

@testset "Single Bernstein polynomial evaluation" begin
    I = Interval(1.0, 3.0)

    # test basis function at x=0.0
    x = I.a
    @test bernstein(I,2,0,x) ≈ 1.0
    @test bernstein(I,2,1,x) ≈ 0.0
    @test bernstein(I,2,2,x) ≈ 0.0

    # test basis function at x=1.0
    x = I.b
    @test bernstein(I,2,0,x) ≈ 0.0
    @test bernstein(I,2,1,x) ≈ 0.0
    @test bernstein(I,2,2,x) ≈ 1.0

    # test basis function at x=2.0
    x = I.a + (I.b-I.a)/2
    @test bernstein(I,2,0,x) ≈ 0.25
    @test bernstein(I,2,1,x) ≈ 0.5
    @test bernstein(I,2,2,x) ≈ 0.25
end

@testset "Basis evaluation single point" begin
    # test partition of unity different orders
    # for single point evaluation
    I = Interval(1.0, 3.0)
    X = LinRange(I.a, I.b, 4)
    for p in 2:4
        for x in X
            B = bernstein(I, 2, x)
            @test sum(B) .≈ 1
            @test all(B .>= 0.0)
        end
    end
end

@testset "Basis evaluation multiple points" begin
    # test partition of unity different orders
    # for multiple point evaluation
    I = Interval(1.0, 3.0)
    x = LinRange(I.a, I.b, 10)
    for p in 1:4
        B = bernstein(I, p, x)
        @test all(sum(B, dims=1) .≈ 1)
        @test all(B .>= 0.0)
    end
end

end