export QuadLearnModel, train_model, init_cached_model_loader

using ChebyshevApprox

struct QuadLearnModel
    space::TensorProduct{3,SplineSpace{Float64}} # B-spline space of the model
    data::Matrix{Float64} # B-spline coefficients of the model
end

function train_model(; order, interp_p=4, R, θ̂, ẑ, case)

    # construct spline space
    space_R = SplineSpace(interp_p, R, Int[interp_p+1, ones(length(R)-2)..., interp_p+1])
    space_θ̂ = SplineSpace(interp_p, θ̂, Int[interp_p+1, ones(length(θ̂)-2)..., interp_p+1])
    space_ẑ = SplineSpace(interp_p, ẑ, Int[interp_p+1, ones(length(ẑ)-2)..., interp_p+1])
    space = space_R ⨷ space_θ̂ ⨷ space_ẑ

    # compute greville interpolation grid
    x = CartesianProduct(s -> grevillepoints(s), space)

    # generate training data
    traindata, testdata = getdata(order=order, R = x.data[1], θ̂ = x.data[2], ẑ = x.data[3], case=case)
    
    # compute collocation matrices
    B_R = Matrix(ders_bspline_interpolation_matrix(space_R, x.data[1], 1)[1])
    B_θ̂ = Matrix(ders_bspline_interpolation_matrix(space_θ̂, x.data[2], 1)[1])
    B_ẑ = Matrix(ders_bspline_interpolation_matrix(space_ẑ, x.data[3], 1)[1])

    # compute interpolation matrix
    B = KroneckerProduct(B_ẑ, B_θ̂, B_R)

    ## compute spline coefficients that interpolate the trainingset
    # SVD
    #D = svd(B)
    #data = D.V * ((traindata.y * D.U)' ./ D.S)

    # direct
    data = B \ traindata.y'

    return QuadLearnModel(space, data), traindata, testdata
end

function (model::QuadLearnModel)(R::Real, θ̂::Real, ẑ::Real)

    # compute basis functions
    B_r = Matrix(ders_bspline_interpolation_matrix(model.space[1], R, 1)[1])
    B_θ = Matrix(ders_bspline_interpolation_matrix(model.space[2], θ̂, 1)[1])
    B_z = Matrix(ders_bspline_interpolation_matrix(model.space[3], ẑ, 1)[1])
    
    # compute shape functions
    B = KroneckerProduct(B_z, B_θ, B_r)

    # compute and return interpolated moments
    return (B * model.data)'
end

function (model::QuadLearnModel)(R,θ̂,ẑ)

    # compute basis functions
    B_R = Matrix(ders_bspline_interpolation_matrix(model.space[1], R, 1)[1])
    B_θ̂ = Matrix(ders_bspline_interpolation_matrix(model.space[2], θ̂, 1)[1])
    B_ẑ = Matrix(ders_bspline_interpolation_matrix(model.space[3], ẑ, 1)[1])
    
    # compute shape functions
    B = KroneckerProduct(B_ẑ, B_θ̂, B_R)

    # compute and return interpolated moments
    return (B * model.data)'
end

(model::QuadLearnModel)(x::CartesianProduct{3}) = model(x.data[1], x.data[2], x.data[3])

function Base.empty!(cache::LFUDA) 
    for key in keys(cache)
        delete!(cache, key)
    end
end

# cached model loading with LFU cache (least frequently used) cache
function init_cached_model_loader(; maxsize = 75)
    @memoize LFUDA{Tuple{NTuple{3, Int64},@NamedTuple{case::Int64}},QuadLearnModel}(maxsize=maxsize) function get_model(R::Int, q::Int, p::Int; case::Int)
        filename_specifier = string("_R", R, "_q", q, "_p", p, "_case", case)
        FileIO.load(string("$(@__DIR__)/../models/data/modeldata", filename_specifier, ".jld2"), "model")
    end
end