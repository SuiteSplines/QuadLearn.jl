using Revise, Test

using QuadLearnData, UnivariateSplines, CartesianProducts, KroneckerProducts
using LinearAlgebra, Plots, StatsBase, StatsPlots
using Statistics
using DataFrames
using Dates

using FileIO, JLD2


# choose wheter to load an existing model or train a new model
train = true

q = 3     # polynomial degree of the target space for quadrature
p = q*q     # polynomial degree of the spline interpolations space
R = 5.0   # radius interval [R, R + 1.0], note: Float64(R) === Int(R)!
#case = 1 # selected cut case 

for case = 1:15
    # this will characterize the filename of the saved files
    filename_specifier = string("_R", Int(R), "_q", q, "_p", p, "_case", case)

    if train == true
        @info "Training $filename_specifier @ $(Dates.format(now(), "HH:MM"))"

        if case == 1
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(30,), p_l=(4.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(30,), p_l=(4.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=1
            )
        
        elseif case == 2
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(6.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(6.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=2
            )
        
        elseif case == 3
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(4.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=3
            )
        
        elseif case == 4
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(4.0,), breakpoints=(0.0, 1.0)),
                case=4
            )
        
        elseif case == 5
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(1.5,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=5
            )
        
        elseif case == 6
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(4.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=6
            )
        
        elseif case == 7
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(1.5,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=7
            )
        
        elseif case == 8
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(2.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(2.0,), breakpoints=(0.0, 1.0)),
                case=8
            )
        
        elseif case == 9
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(4.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(4.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=9
            )
        elseif case == 10
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 50),
                θ̂=refined_increasing_vector(n=(20,), p_l=(8.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(20,), p_l=(8.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=10
            )
        
        elseif case == 11
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(15,), p_l=(2.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(30,), p_l=(6.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=11
            )
        
        elseif case == 12
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(4.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=12
            )
        
        elseif case == 13
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(4.0,), breakpoints=(0.0, 1.0)),
                case=13
            )
        elseif case == 14
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                R=IncreasingRange(R, R+1.0, 15),
                θ̂=refined_increasing_vector(n=(25,), p_l=(2.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                ẑ=refined_increasing_vector(n=(25,), p_l=(1.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                case=14
            )
        
        elseif case == 15
        
            model, traindata,  testdata = train_model(;
                order=q,
                interp_p=p,
                #R=IncreasingRange(R, R+1.0, 10),
                #θ̂=refined_increasing_vector(n=(20,), p_l=(2.0,), p_r=(1.0,), breakpoints=(0.0, 1.0)),
                #ẑ=refined_increasing_vector(n=(20,), p_l=(1.0,), p_r=(1.5,), breakpoints=(0.0, 1.0)),
                R=IncreasingRange(R, R+1.0, 30),
                θ̂=IncreasingRange(0.0, 1.0, 30),
                ẑ=IncreasingRange(0.0, 1.0, 30),
                case=15
            )

        else
            error("invalid case id")
        end
        # save model
        FileIO.save(string("models/data/traindata", filename_specifier, ".jld2"), "traindata", traindata)
        FileIO.save(string("models/data/testdata", filename_specifier, ".jld2"), "testdata", testdata)
        FileIO.save(string("models/data/modeldata", filename_specifier, ".jld2"), "model", model)
    else
        # load data
        traindata = FileIO.load(string("models/data/traindata", filename_specifier, ".jld2"), "traindata")
        testdata = FileIO.load(string("models/data/testdata", filename_specifier, ".jld2"), "testdata")
        model = FileIO.load(string("models/data/modeldata", filename_specifier, ".jld2"), "model")
    end

    @info "... postprocessing $filename_specifier"
    
    # Evaluate the model on the training input data and compute the difference. 
    # This should be within machine precision (least squares opti)
    x_train = traindata.input
    y_train = model(x_train)
    @test isapprox(norm(y_train - traindata.y), 0.0; atol=1e-12)

    # Evaluate the model on the testdata input data and compute the difference. 
    # Hopefully this is in the order of 10-5 - 10-8 (at "midpoints" of parameter space)
    x_test = testdata.input
    y_test = model(x_test)
    @test isapprox(norm(y_test - testdata.y), 0.0; atol=1e-2)

    # compute the testerror and reshape it into a 3D array such that we can make some
    # slices and make nice plots
    testerror = y_test - testdata.y
    trainerror = y_train - traindata.y
    R_test, R_train = testdata.input.data[1], traindata.input.data[1]
    θ̂_test, θ̂_train = testdata.input.data[2], traindata.input.data[2]
    ẑ_test, ẑ_train = testdata.input.data[3], traindata.input.data[3]
    k = 3

    T_test = reshape(testdata.y[k,:], length.(testdata.input.data))
    t_test = reshape(y_test[k,:], length.(testdata.input.data))
    e_test = reshape(testerror[k,:], length.(testdata.input.data))

    T_train = reshape(traindata.y[k,:], length.(traindata.input.data))
    t_train = reshape(y_train[k,:], length.(traindata.input.data))
    e_train = reshape(trainerror[k,:], length.(traindata.input.data))

    # Take a slice in 'ẑ' and 'θ̂' and plot results at fixed "R"
    # We first show the response surface
    h1_test = Plots.wireframe(ẑ_test, θ̂_test, t_test[1,:,:], xlabel="ẑ", ylabel="θ̂")
    h1_train = Plots.wireframe(ẑ_train, θ̂_train, t_train[1,:,:], xlabel="ẑ", ylabel="θ̂")
    savefig(h1_test, string("models/results/test_landscape", filename_specifier, ".png"))
    savefig(h1_train, string("models/results/train_landscape", filename_specifier, ".png"))

    # Now we show the error (check error over parameters, adapt spaces (probably))
    h2_test = Plots.wireframe(ẑ_test, θ̂_test, e_test[1,:,:], xlabel="ẑ", ylabel="θ̂")
    h2_train = Plots.wireframe(ẑ_train, θ̂_train, e_train[1,:,:], xlabel="ẑ", ylabel="θ̂")
    savefig(h2_test, string("models/results/test_landscape_error", filename_specifier, ".png"))
    savefig(h2_train, string("models/results/train_landscape_error", filename_specifier, ".png"))

    # plot the error for each of the moments as a violin plot
    test_labels = map(k -> repeat(["M$k"], size(testerror, 2)), 1:(q+1)^2)
    train_labels = map(k -> repeat(["M$k"], size(trainerror, 2)), 1:(q+1)^2)

    ylims, ytics = (1e-16, 1e-5), map(k -> 10.0^(-k), 6:16)
    h3_test = Plots.violin(test_labels, abs.(testerror'), leg = false, yscale=:log10, ylims = ylims, yticks = ytics)

    ylims, ytics = (1e-25, 1e-13), map(k -> 10.0^(-k), 14:26)
    h3_train = Plots.violin(train_labels, abs.(trainerror'), leg = false, yscale=:log10, ylims = ylims, yticks = ytics)
    savefig(h3_test,string("models/results/test_model_accuracy", filename_specifier, ".png"))
    savefig(h3_train,string("models/results/train_model_accuracy", filename_specifier, ".png"))
end

