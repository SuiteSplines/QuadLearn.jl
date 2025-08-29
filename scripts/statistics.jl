using Revise, Test

using QuadLearnData, UnivariateSplines, CartesianProducts, KroneckerProducts
using LinearAlgebra, Plots, StatsBase, StatsPlots
using Statistics
using DataFrames

using FileIO, JLD2

# radius
R = 2.0

err_mean, err_std, err_max = Float64[], Float64[], Float64[]

for case = 1:15
    # this characterizes the filename of the model
    filename_specifier = string("_R", Int(R), "_q", q, "_p", p, "_case", case)

    # load data
    testdata = FileIO.load(string("models/data/testdata", filename_specifier, ".jld2"), "testdata")
    model = FileIO.load(string("models/data/modeldata", filename_specifier, ".jld2"), "model")

    # compute the testerror
    x = testdata.input
    y = model(x)
    testerror = y - testdata.y

    push!(err_mean, mean(testerror))
    push!(err_std, std(testerror))
    push!(err_max, maximum(abs.(testerror)))
end



# statistics
DataFrame(
    "mean(testerror)" => err_mean,
    "std(testerror)" => err_std,
    "max(abs(testerror))" => err_max,
) 