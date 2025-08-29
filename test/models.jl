module TestModels

using Test

using LinearAlgebra, QuadLearnData
using SortedSequences: IncreasingRange
using Algoim, AbstractMappings, IgaBase, ImmersedSplines

@testset "Check interpolation of data" begin
    
    # train a small model using spline interpolation
    model, traindata,  testdata = train_model(;
        order=2,
        interp_p=4,
        R=IncreasingRange(2.0, 3.0, 10),
        θ̂=refined_increasing_vector(n=(20,), p_l=(1.0,), p_r=(4.0,), breakpoints=(0.0, 1.0)),
        ẑ=refined_increasing_vector(n=(20,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
        case=6
    )

    # evaluate the model on the training input data
    x = traindata.input
    y = model(x)

    # the interpolation data is exactly interpolated
    # so error should be to machine precision
    @test all(abs.(y - traindata.y) .< 10e-15)
end

end # TestModels