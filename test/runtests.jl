using Test

tests = [
    "bernstein",
    "trainingdata",
    "models",
    "parametrization",
    "patchaccessor",
]
isdir("$(@__DIR__)/../models/data") && push!(tests, "quadrature")

@testset "QuadLearnData Unit-testing" begin

    # run all unit-tests
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end

end # @testset
