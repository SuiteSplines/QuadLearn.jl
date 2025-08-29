using Revise, Test

using QuadLearnData, UnivariateSplines, CartesianProducts, KroneckerProducts
using LinearAlgebra #, Plots, StatsBase, StatsPlots
using GLMakie
using Statistics
using DataFrames
using Dates
using Printf

using FileIO, JLD2

q = 2    # polynomial degree of the target space for quadrature
p = q*q  # polynomial degree of the spline interpolations space
R = 5.0  # radius interval [R, R + 1.0], note: Float64(R) === Int(R)!
case = 6 # selected cut case 

# this will characterize the filename of the saved files
filename_specifier = string("_R", Int(R), "_q", q, "_p", p, "_case", case)

# load case data
traindata = FileIO.load(string("models/data/traindata", filename_specifier, ".jld2"), "traindata")
testdata = FileIO.load(string("models/data/testdata", filename_specifier, ".jld2"), "testdata")
model = FileIO.load(string("models/data/modeldata", filename_specifier, ".jld2"), "model")
    
# Evaluate the model on the training input data and compute the difference. 
# This should be within machine precision (least squares opti)
x_train = traindata.input
y_train = model(x_train)
@test isapprox(norm(y_train - traindata.y), 0.0; atol=1e-12)

# Evaluate the model on the testdata input data and compute the difference. 
# Hopefully this is in the order of 10-5 - 10-8 (at "midpoints" of parameter space)
x_test = testdata.input
y_test = model(x_test)
#@test isapprox(norm(y_test - testdata.y), 0.0; atol=1e-5)

# compute the testerror and reshape it into a 3D array such that we can make some
# slices and make nice plots
testerror = y_test - testdata.y
trainerror = y_train - traindata.y

# get sampling points for test and train data
R_test, R_train = testdata.input.data[1], traindata.input.data[1]
θ̂_test, θ̂_train = testdata.input.data[2], traindata.input.data[2]
ẑ_test, ẑ_train = testdata.input.data[3], traindata.input.data[3]
k = 3

# collect k'th moment data in test
T_test = reshape(testdata.y[k,:], length.(testdata.input.data))
t_test = reshape(y_test[k,:], length.(testdata.input.data))
e_test = reshape(testerror[k,:], length.(testdata.input.data))
f_T_test(k) = reshape(testdata.y[k,:], length.(testdata.input.data))
f_t_test(k) = reshape(y_test[k,:], length.(testdata.input.data))
f_e_test(k) = reshape(testerror[k,:], length.(testdata.input.data))

# collect k'th moment data in training
T_train = reshape(traindata.y[k,:], length.(traindata.input.data))
t_train = reshape(y_train[k,:], length.(traindata.input.data))
e_train = reshape(trainerror[k,:], length.(traindata.input.data))
f_T_train(k) = reshape(traindata.y[k,:], length.(traindata.input.data))
f_t_train(k) = reshape(y_train[k,:], length.(traindata.input.data))
f_e_train(k) = reshape(trainerror[k,:], length.(traindata.input.data))

# Take a slice in (ẑ,θ̂) plane and plot results at fixed radius R

scientific_latex(x) = begin
    base, exp = split(@sprintf("%.1e", x), "e")
    exp_val = parse(Int, exp)  # this removes leading zeros automatically
    #latex_str = L"" * base * L"\cdot 10^{" * string(exp_val) * L"}"
    latex_str = L"%$(base) \times 10^{ %$(string(exp_val)) }"
end

# response surface wireframe
GLMakie.activate!()
with_theme(theme_latexfonts()) do
    f = Figure(fontsize=35, figure_padding=100)
    resize!(f, 2400, 1800)
    #display(f)
    linds = LinearIndices((q+1,q+1))
    cinds = CartesianIndices((q+1,q+1))
    gridcolor = RGBf(0.7,0.7,0.7)
    for k = 1:(q+1)^2
        ax = Axis3(f[cinds[k][1], cinds[k][2]],
            zticks=WilkinsonTicks(4, k_min=3),
            ztickformat = v -> scientific_latex.(v),
            zlabeloffset = 115,
            xlabelrotation = 0,
            ylabelrotation = 0,
            zlabelrotation = 0,
            xlabeloffset = 60,
            ylabeloffset = 60,
            azimuth = π/8,
            elevation = π/8,
            title = L"\alpha = %$k",
            xgridcolor = gridcolor,
            ygridcolor = gridcolor,
        )
        ax.xlabel = L"\hat{\theta}"
        ax.ylabel = L"\hat{z}"
        ax.zlabel = ""
        #Makie.wireframe!(ax, θ̂_train, ẑ_train, f_T_train(linds[k])[1,:,:], color=:black)
        Makie.wireframe!(ax, θ̂_test, ẑ_test, f_e_test(linds[k])[1,:,:], color=:black)
    end
    for k = 1:q
        colgap!(f.layout, k, Relative(0.025))
        rowgap!(f.layout, k, Relative(0.025))
    end
    save("test.pdf", f; update=true, backend=CairoMakie)
end


# collect labels
test_labels = vcat(map(k -> repeat(["M$k"], size(testerror, 2)), 1:(q+1)^2)...)
train_labels = vcat(map(k -> repeat(["M$k"], size(trainerror, 2)), 1:(q+1)^2)...)
abs.(testerror')[:]
#
#    ylims, ytics = (1e-16, 1e-5), map(k -> 10.0^(-k), 6:16)
#    h3_test = Plots.violin(test_labels, abs.(testerror'), leg = false, yscale=:log10, ylims = ylims, yticks = ytics)


# test error violin plot
GLMakie.activate!()
with_theme(theme_latexfonts()) do
    f = Figure(fontsize=18)
    display(f)
    resize!(f, 350,350)
    gridcolor = RGBf(0.7,0.7,0.7)
    ax = Axis(f[1,1],
        yscale=log10,
        xticks = 1:(q+1)^2,
        yticks = LogTicks(WilkinsonTicks(11,k_min=11,k_max=11)),
        title = L"\text{Case %$case}",
        xgridcolor = gridcolor,
        ygridcolor = gridcolor,
    )
    ylims!(10^-19, 10^-6)
    ax.xlabel = L"\alpha"

    colors = resample_cmap(:tab10, 9)
    test_colors = vcat([fill(colors[k], size(testerror,2)) for k in 1:(q+1)^2]...)
    test_labels = vcat([fill(k, size(testerror,2)) for k in 1:(q+1)^2]...)
    violin!(ax, test_labels, abs.(testerror'[:]); scale=:width, datalimits=(0,Inf), color = test_colors, strokewidth=1.5)
    save("test.pdf", f; update=true, backend=CairoMakie)
end
