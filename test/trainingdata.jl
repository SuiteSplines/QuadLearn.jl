module TestTrainingData

using Test, Revise

using SortedSequences, CartesianProducts, QuadLearnData, StaticArrays  
using StaticArrays, LinearAlgebra

@testset "Check bernstein integrals uncut element" begin
    element = Interval(-1.0,1.0) ⨱ Interval(-1.0,1.0)
    M = integrate_bernstein_basis(4, element)
    @test sum(M) ≈ 4.0
end

@testset "Check bernstein integrals cut element" begin
    element = Interval(-1.0,1.0) ⨱ Interval(-1.0,1.0)
    M = integrate_bernstein_basis(4, implicit_circle_def(3.0), element)
    @test sum(M) ≈ 4.0

    r = 0.75
    M = integrate_bernstein_basis(4, implicit_circle_def(r), element; n=64)
    @test isapprox(sum(M) - pi * r^2, 0.0, atol=1e-13)
end

end # TestTrainingData