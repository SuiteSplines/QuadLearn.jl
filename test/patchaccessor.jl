#function get_adaptive_weights(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, constrain::Bool=true)
#    try
#        # least square circular approximation of the cut
#        R, a, b = levelset_least_squares_circle(ϕ, element; n=32, constrain=constrain)
#
#        # compute reference cut parameters
#        R_ref, θ_ref, z_ref = get_reference_cut_parameters(R, a, b, element)
#
#        # get transformation rule
#        θ_sym, transformation_rule = symmetry_transform(θ_ref)
#
#        # get reference element domain in training range
#        element_ref = get_cut_domain(1.0, θ_sym, z_ref)
#
#        # get parametrization
#        p = Parametrization(R_ref, 1.0)
#
#        # get inverse parameterization
#        p⁻¹ = inv(p)
#
#        # get reference element case
#        case = parameter_case(p, θ_sym, z_ref)
#
#        # compute reference element (θ̂, ẑ) normalized parameter pair
#        local θ̂_ref, ẑ_ref
#        try 
#            θ̂_ref, ẑ_ref = p⁻¹(θ_sym, z_ref; case=case)
#        catch e
#            if isa(e, SingularTheta)
#                θ̂_ref, ẑ_ref = 0.5, e.ẑ
#            else
#                rethrow(e)
#            end
#        end
#
#        # get model
#        model = loader(floor(Int, R_ref), order, interp_order; case=case)
#
#        # evaluate model
#        M_reference = reshape(model(R_ref, θ̂_ref, ẑ_ref), order+1, order+1)
#
#        # apply transformation
#        M_transformed = transformation_rule(M_reference)
#
#        # compute Jacobian
#        J = QuadLearnData.measure(element)
#
#        # scale moments
#        M_scaled = M_transformed * J
#
#        # test complement condition
#        complement = test_complement_condition(ϕ, element, R, a, b; phase=phase)
#
#        # evaluate complement if needed
#        M = complement ? integrate_bernstein_basis(order, element) .- M_scaled : M_scaled
#
#        # compute linear map from weights to moments (memoized!)
#        A = get_bernsteinfuns_at_gausspoints(; order=order)
#
#        # compute test weights
#        w = A \ M[:]
#    catch e
#        if isa(e, StrictlyPositiveLevelSet)
#            w = (phase == 1) ? get_weights(element; order=order+1)[:] : zeros((order+1)^2)
#        elseif isa(e, StrictlyNegativeLevelSet)
#            w = (phase == -1) ? get_weights(element; order=order+1)[:] : zeros((order+1)^2)
#        else
#            rethrow(e)
#        end
#    end
#end
#
#function get_adaptive_weights(ϕ::AlgoimCallLevelSetFunction, element::Element{2}; loader::Function, order::Int, interp_order::Int, n::Int=32, phase::Int=-1, constrain::Bool=true)
#    # get element domain
#    Ωₑ = get_element_domain(element)
#
#    # get weights
#    get_weights(ϕ, Ωₑ; loader=loader, order=order, interp_order=interp_order, n=n, phase=phase, constrain=constrain)
#end