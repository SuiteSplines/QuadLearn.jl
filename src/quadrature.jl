using UnivariateSplines

export get_stencil, get_weights, evaluate_distance_function
export get_patch_weights, get_patch_stencil

UnivariateSplines.GaussRule(method, n, I::Interval) = GaussRule(method, n, I.a, I.b)

function get_stencil(I::Interval; order)
    qx = GaussRule(Legendre, order, I)
    return qx.x
end

function get_stencil(element::CartesianProduct{2}; order)
    return get_stencil(element.data[1], order=order) ⨱ get_stencil(element.data[2], order=order)
end

function get_patch_stencil(partition::CartesianProduct{2}; order)
    univariate_patch_stencils = ntuple(d -> PatchRule(partition.data[d]; npoints=order).x[2:end-1], 2)
    CartesianProduct(univariate_patch_stencils...)
end

function get_weights(I::Interval; order)
    qx = GaussRule(Legendre, order, I)
    return qx.w
end

function get_weights(element::CartesianProduct{2}; order)
    wx, wy = get_weights.(element.data; order=order)
    wx * wy'
end

function get_weights(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, constrain::Bool=true)
    try
        # least square circular approximation of the cut
        R, a, b = levelset_least_squares_circle(ϕ, element; n=32, constrain=constrain)

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

function get_weights(ϕ::AlgoimCallLevelSetFunction, element::Element{2}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, constrain::Bool=true)
    # get element domain
    Ωₑ = get_element_domain(element)

    # get weights
    get_weights(ϕ, Ωₑ; loader=loader, order=order, interp_order=interp_order, n=n, phase=phase, constrain=constrain)
end

function get_patch_weights(ϕ::AlgoimCallLevelSetFunction, partition::CartesianProduct{2}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, constrain::Bool=true)
    # elements iterator
    elements = Elements(partition)

    # number of elements
    nelem = length(elements)

    # array of weights on each element (columnwise)
    weights = zeros((order+1)^2, nelem)

    # collect element weights
    for (lind, element) in enumerate(elements)
        # get weights on element
        w = get_weights(ϕ, element; loader=loader, order=order, interp_order=interp_order, n=n, phase=phase, constrain=constrain)

        # store weights (todo: implement get_weights!(v, ....))
        weights[:,lind] = w
    end

    # reshape to fit CartesianProduct quadrature points
    (n, m), q = size(elements), order + 1
    blocks = reshape(weights, q, q, n, m)
    permuted = permutedims(blocks, (1, 3, 2, 4))
    W = reshape(permuted, n*q, m*q)
end

function evaluate_distance_function(r, x, y)
    return x.^2 + y.^2 .- r^2;
end