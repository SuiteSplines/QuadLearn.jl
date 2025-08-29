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

case_data = []
for case = 1:15
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

    # collect data
    push!(case_data, (
        testerror = testerror,
    ))
end

scientific_latex(x) = begin
    base, exp = split(@sprintf("%.1e", x), "e")
    exp_val = parse(Int, exp)  # this removes leading zeros automatically
    #latex_str = L"" * base * L"\cdot 10^{" * string(exp_val) * L"}"
    latex_str = L"%$(base) \times 10^{ %$(string(exp_val)) }"
end

# test error violin plot
GLMakie.activate!()
with_theme(theme_latexfonts()) do
    f = Figure(fontsize=25, figure_padding=20)
    resize!(f, 1800, 2000)
    linds = LinearIndices((5,3))
    cinds = CartesianIndices((5,3))
    gridcolor = RGBf(0.7,0.7,0.7)
    case = 1
    for case = 1:15
        ax = Axis(f[cinds[case][2], cinds[case][1]],
            yscale=log10,
            xticks = 1:(q+1)^2,
            yticks = LogTicks(WilkinsonTicks(11,k_min=11,k_max=11)),
            title = L"\text{Case %$case}",
            xgridcolor = gridcolor,
            ygridcolor = gridcolor,
        )
        ylims!(10^-19, 10^-5)
        ax.xlabel = L"\alpha"

        colors = resample_cmap(:tab10, 9)
        test_colors = vcat([fill(colors[k], size(case_data[case].testerror,2)) for k in 1:(q+1)^2]...)
        test_labels = vcat([fill(k, size(case_data[case].testerror,2)) for k in 1:(q+1)^2]...)
        violin!(ax, test_labels, abs.((case_data[case].testerror')[:]); scale=:width, datalimits=(0,Inf), color = test_colors, strokewidth=1.5)
    end
    display(f)
    #save("test.png", f; update=true)
    save("test.pdf", f; update=true, backend=CairoMakie)
end





using StatsPlots, DataFrames
df = DataFrame(A = 1:10, B =rand(10))
plotd = @df df StatsPlots.plot(:A, :B);
savefig(plotd,"test.pdf")

# plot the error for each of the moments as a violin plot

begin
    Plots.default(
        fontfamily = "Computer Modern",
        guidefont = (26,),
        titlefont = (26,),
        tickfont = (24,),
    )
    layout = Plots.@layout [grid(5, 3)]
    f = Plots.plot(layout = layout, size = (2600, 3300))
    ylims, ytics = (1e-19, 1e-5), map(k -> 10.0^(-k), 6:2:19)
    for case = 1:15
        test_labels = map(k -> repeat(["$k"], size(case_data[case].testerror, 2)), 1:(q+1)^2)
        h3_test = Plots.violin!(f, test_labels, abs.(case_data[case].testerror'),
            leg = false,
            yscale=:log10,
            ylims = ylims,
            yticks = ytics,
            linewidth = 2,
            framestyle = :box,
            #xguidefontvalign = :top,
            xlabel = (case == 14) ? L"\alpha" : "",
            ylabel = (case == 7) ? "Pointwise error distribution" : "",
            yguidefonthalign = :left,
            left_margin = (case == 7) ? 3.0Plots.cm : 2Plots.cm,
            bottom_margin = (case == 14) ? 1.5Plots.cm : 0.4Plots.cm,
            title = "Case $case",
            xtickfontvalign = :bottom,
            ytickfonthalign = :center,
            gridcolor = RGB(0.1,0.1,0.1),
            gridlinewidth = 2,
            subplot=case,
        )
    end
    savefig(h3_test,"test.pdf")
end

