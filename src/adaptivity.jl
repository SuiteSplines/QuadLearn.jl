import RegionTrees: AbstractRefinery, Cell, needs_refinement, refine_data, child_boundary, HyperRectangle, adaptivesampling!, allleaves

export get_weights, CurvatureTooLarge, QuadTreeRefinery, AdaptiveDataQuadRule, AdaptiveCutcellQuadratureRule

struct QuadTreeRefinery <: AbstractRefinery
    ϕ::AlgoimCallLevelSetFunction
    Rmin::Float64
    n::Int
    factor::Float64
end

hyperrecangle_to_element(r::HyperRectangle{2,T}) where {T} = Interval(r.origin[1], r.origin[1] + r.widths[1]) ⨱ Interval(r.origin[2], r.origin[2] + r.widths[2])

function needs_refinement(r::QuadTreeRefinery, cell)
    isinf(cell.data.R) && false
    cell.data.R / cell.boundary.widths[1] < r.Rmin * r.factor
end

function refine_data(r::QuadTreeRefinery, cell::Cell, indices)
    boundary = child_boundary(cell, indices)
    element = hyperrecangle_to_element(boundary)

    try
        R, a, b = levelset_least_squares_circle(r.ϕ, element; n=r.n, linear_max_constrain=true)
        return (R=R, a=a, b=b)
    catch e
        if isa(e, StrictlyPositiveLevelSet)
            return (R=Inf, a=Inf, b=Inf)
        elseif isa(e, StrictlyNegativeLevelSet)
            return (R=Inf, a=Inf, b=Inf)
        else
            rethrow(e)
        end
    end
end

struct CurvatureTooLarge <: Exception end 

function get_weights(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}, R::Float64, a::Float64, b::Float64; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, Rmax::Float64=30.0)
    try
        # width and height of element
        width, height = measure.(element.data)

        # elements must be square
        @assert isapprox(width, height, atol=10e-13) "element must be square"

        # compute scaling factor
        s = 1 / width

        # check condition for refinement
        (R*s < 1.0) && throw(CurvatureTooLarge())

        # compute reference cut parameters
        R_ref, θ_ref, z_ref = get_reference_cut_parameters(R, a, b, element)

        # get transformation rule
        θ_sym, transformation_rule = symmetry_transform(θ_ref)

        # get reference element domain in training range
        element_ref = get_cut_domain(1.0, θ_sym, z_ref)

        # get parametrization
        p = Parametrization(R_ref, 1.0)

        # get inverse parameterization
        p⁻¹ = inv(p)

        # get reference element case
        case = parameter_case(p, θ_sym, z_ref)

        # In an ideal world, where the cut contour is approximated
        # extremely well by a circle, `case` can never be == -1:
        # Least squares will throw an exception if element is not cut
        # and with that indicate in which phase the element is contained.
        #
        # BUT: the approximated contour is not exact and in some edge 
        # cases (elements that are barely cut), the approximated cut
        # contour might not intersect the element anymore and the cut case
        # test will return -1. This is more likely to happen on very fine
        # meshes but nonetheless it is quite rare.
        #
        # To handle that we can query the levelset value in element
        # center (furthest from the element boundary), check the sign
        # of the levelst at that point and throw StrictlyPositiveLevelSet or
        # StrictlyNegativeLevelSet exception. The error will be negligible.
        if case == -1 
            # center of element
            C = center(element)

            # sign of levelset value at element center
            test_sign = sign(ϕ(C))

            # levelset is positive
            test_sign == 1 && throw(StrictlyPositiveLevelSet())
            
            # levelset is negative
            test_sign == -1 && throw(StrictlyNegativeLevelSet())

            # levelset is zero? bug...
            error("cut case identification failed")
        end

        # compute reference element (θ̂, ẑ) normalized parameter pair
        local θ̂_ref, ẑ_ref
        try 
            θ̂_ref, ẑ_ref = p⁻¹(θ_sym, z_ref; case=case)
        catch e
            if isa(e, SingularTheta)
                θ̂_ref, ẑ_ref = 0.5, e.ẑ
            else
                rethrow(e)
            end
        end

        # get model
        model = loader(floor(Int, R_ref), order, interp_order; case=case)

        # evaluate model
        M_reference = reshape(model(R_ref, θ̂_ref, ẑ_ref), order+1, order+1)

        # apply transformation
        M_transformed = transformation_rule(M_reference)

        # compute Jacobian
        J = QuadLearnData.measure(element)

        # scale moments
        M_scaled = M_transformed * J

        # test complement condition
        complement = test_complement_condition(ϕ, element, R, a, b; phase=phase)

        # evaluate complement if needed
        M = complement ? integrate_bernstein_basis(order, element) .- M_scaled : M_scaled

        # compute linear map from weights to moments (memoized!)
        A = get_bernsteinfuns_at_gausspoints(; order=order)

        # compute test weights
        w = A \ M[:]
    catch e
        if isa(e, StrictlyPositiveLevelSet)
            w = (phase == 1) ? get_weights(element; order=order+1)[:] : zeros((order+1)^2)
        elseif isa(e, StrictlyNegativeLevelSet)
            w = (phase == -1) ? get_weights(element; order=order+1)[:] : zeros((order+1)^2)
        else
            rethrow(e)
        end
    end
end

struct AdaptiveDataQuadRule{Dim, X, W}
    x::X
    w::W
    phase::Int
    function AdaptiveDataQuadRule(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{Dim}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, Rmax::Float64=30.0, factor::Float64=1.0) where {Dim}
        x, w = fill_quad_data(ϕ, element; loader=loader, order=order, interp_order=interp_order, n=n, phase=phase, factor=factor)        
        X = typeof(x)
        W = typeof(w)
        return new{Dim,X,W}(x, w, phase) 
    end
end

function AdaptiveDataQuadRule(ϕ::AlgoimCallLevelSetFunction, element::Element{2}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, Rmax::Float64=30.0, factor::Float64=1.0)
    # get element domain
    Ωₑ = get_element_domain(element)

    # initialize
    AdaptiveDataQuadRule(ϕ, Ωₑ; loader=loader, order=order, interp_order=interp_order, n=n, phase=phase, Rmax=Rmax, factor=factor)
end

function fill_quad_data(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, Rmax::Float64=30.0, factor::Float64 = 1.0)
    # define adaptive refinery
    refinery = QuadTreeRefinery(ϕ, 1.0, n, factor)

    # define root cell
    origin = SVector(element.data[1][1], element.data[2][1])
    widths = SVector(measure.(element.data))
    local R, a, b
    try
        R, a, b = levelset_least_squares_circle(refinery.ϕ, element; n=refinery.n, linear_max_constrain=true)
    catch e
        if isa(e, StrictlyPositiveLevelSet)
            w = (phase == 1) ? get_weights(element; order=order+1)[:] : zeros((order+1)^2)
            x = get_stencil(element; order = order+1)[:]
            return x, w
        elseif isa(e, StrictlyNegativeLevelSet)
            w = (phase == -1) ? get_weights(element; order=order+1)[:] : zeros((order+1)^2)
            x = get_stencil(element; order = order+1)[:]
            return x, w
        else
            rethrow(e)
        end
    end
    root = Cell(origin, widths, (R=R,a=a,b=b))

    # perform adaptive sampling
    adaptivesampling!(root, refinery)

    # collect data
    subcells = collect(allleaves(root))
    nsubcells = length(subcells)
    npts = (order+1)^2
    x = Matrix{SVector{2,Float64}}(undef, npts, nsubcells)
    w = Matrix{Float64}(undef, npts, nsubcells)
    for (k,cell) in enumerate(subcells)
        e = QuadLearnData.hyperrecangle_to_element(cell.boundary)
        x[:,k] = SVector.(get_stencil(e; order = order+1))[:]
        w[:,k] = get_weights(ϕ, e; loader=loader, order=order, interp_order=interp_order, n=n, phase=phase)
    end

    # return a vector of quadrature points and corresponding weights
    return x[:], w[:]
end

Base.length(qr::AdaptiveDataQuadRule) = Base.length(qr.w)
Base.size(qr::AdaptiveDataQuadRule) = (length(qr),)





mutable struct AdaptiveCutcellQuadratureRule{Dim}
    partition::CartesianProduct{Dim}
    mapping::AlgoimCallLevelSetFunction
    npoints::Int
    loader::Function
    order::Int
    interp_order::Int
    n::Int
    Rmax::Float64
    factor::Float64
    element::Element
    quadrule::AdaptiveDataQuadRule
    function AdaptiveCutcellQuadratureRule(; partition::CartesianProduct{Dim}, mapping, npoints, loader, order, interp_order, n = 32, Rmax = 30.0, factor = 1.0) where {Dim}
        new{Dim}(partition, mapping, npoints, loader, order, interp_order, n, Rmax, factor)
    end
end

Base.isassigned(Q::AdaptiveCutcellQuadratureRule) = isdefined(Q, :element) && isdefined(Q, :quadrule)

function istodate(Q::AdaptiveCutcellQuadratureRule, e::Element, phase)::Bool
    
    # check if struct is already initialized
    if !isassigned(Q)
        return false
    end

    # check if the same quadrature rule is stored
    if Q.element==e && Q.quadrule.phase==phase
        return true
    end

    return false
end

function update!(Q::AdaptiveCutcellQuadratureRule, e::Element, phase::Int=-1)
    if !istodate(Q, e, phase)
        Q.element = e # element and quadrule are always updated simultaneously
        Q.quadrule = AdaptiveDataQuadRule(Q.mapping, e; loader=Q.loader, order=Q.order, interp_order=Q.interp_order, n=Q.n, phase=phase, factor=Q.factor)

    end
end

struct AdaptiveCutcellData{Dim, T, Q<:AdaptiveCutcellQuadratureRule}
    trialspace::TensorProduct{Dim,SplineSpace{T}}
    testspace::TensorProduct{Dim,SplineSpace{T}}
    quadrule::Q
end

Base.ndims(::AdaptiveCutcellData{Dim}) where Dim = Dim

function IgaFormation.ElementAccessorData(partition, trialspace, testspace, quadrule::AdaptiveCutcellQuadratureRule)
    return AdaptiveCutcellData(trialspace, testspace, quadrule)
end

function IgaFormation.QuadratureRule(data::AdaptiveCutcellData,  e::Element; phase::Int=-1)
    update!(data.quadrule, e, phase) # update algoim data 
    return data.quadrule.quadrule
end

function IgaFormation.QuadraturePoints(data::AdaptiveCutcellData,  e::Element; phase::Int=-1)
    Q = QuadratureRule(data, e; phase=phase)
    return Q.x
end

function IgaFormation.QuadratureWeights(data::AdaptiveCutcellData,  e::Element; phase::Int=-1)
    Q = QuadratureRule(data, e; phase=phase)
    return Q.w
end

function IgaFormation.TrialFunctions(data::AdaptiveCutcellData{Dim}, e::Element, i::Int; ders) where Dim
    x = single_quadrature_point(data, e, i) 
    eindex = IgaFormation.get_element_indices(e)
    funs = ntuple(k -> get_element_functions(data.trialspace[k], eindex[k], x[k], ders[k]), Dim)
    return KroneckerProduct(funs..., reverse=true)
end

function IgaFormation.TestFunctions(data::AdaptiveCutcellData{Dim}, e::Element, i::Int; ders) where Dim
    x = single_quadrature_point(data, e, i) 
    eindex = IgaFormation.get_element_indices(e)
    funs = ntuple(k -> get_element_functions(data.testspace[k], eindex[k], x[k], ders[k]), Dim)
    return KroneckerProduct(funs..., reverse=true)
end

function single_quadrature_point(data::AdaptiveCutcellData{Dim},  e::Element{Dim}, i::Int) where Dim
    @assert data.quadrule.element==e # check if data is up to date
    return data.quadrule.quadrule.x[i]
end

function single_quadrature_point(data::AdaptiveCutcellData{2},  e::Element{1}, i::Int)
    @assert data.quadrule.element==e # check if data is up to date
    k = IgaBase.find_singular_dimension(map(length, e.parent.indices)...)
    x = data.quadrule.quadrule.x[i]
    s = e[1][k]
    if k==1
        return (s, x[1])
    elseif k==2
        return (x[1], s)
    end
end

function single_quadrature_point(data::AdaptiveCutcellData{3},  e::Element{2}, i::Int)
    @assert data.quadrule.element==e # check if data is up to date
    k = IgaBase.find_singular_dimension(map(length, e.parent.indices)...)
    x = data.quadrule.quadrule.x[i]
    s = e[1,1][k]
    if k==1
        return (s, x[1], x[2])
    elseif k==2
        return (x[1], s, x[2])
    elseif k==3
        return (x[1], x[2], s)
    end
end

function get_element_functions(S::SplineSpace, e::Int, u, ders::Int)
    return dersbsplinebasisfuns(Degree(S), KnotVector(S), u, ders+1)[:,ders+1]
end
