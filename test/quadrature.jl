module TestQuadrature

using Test
using QuadLearnData, SortedSequences, IgaFormation
using SortedSequences: IncreasingRange
using Algoim, AbstractMappings, IgaBase, ImmersedSplines
using LinearAlgebra

@testset "Test pipeline for area integration on exact circular cut" begin
    # circular cut radius
    radius = 2.0

    # embedding box width
    width = 4.0

    # levelset function for circle at origin with R = radius
    ϕ = implicit_circle_def(radius)

    # number of elements in each parametric direction
    nelem = 4 # (2*width / nelem) ≥ radius!

    # partition of [-width/2,width/2] × [-width/2,width/2]
    partition = IncreasingRange(-width/2, width/2, nelem) ⨱ IncreasingRange(-width/2, width/2, nelem)

    # model loader
    loader = init_cached_model_loader(; maxsize=75)

    # order of quadrature
    order = 2

    # order of interpolation model
    interp_order = 4

    # phase to integrate (1: outside, -1: inside)
    phase = 1

    # total area
    area = 0.0

    # element loop
    for element in Elements(partition)
        # get element domain
        Ωₑ = get_element_domain(element)

        try
            # least square circular approximation of the cut
            #R, a, b = levelset_least_squares_circle(ϕ, Ωₑ; n=32)
            R, a, b = radius, 0.0, 0.0 # in this case cut countour is exact!

            # compute reference cut parameters
            R_ref, θ_ref, z_ref = get_reference_cut_parameters(R, a, b, Ωₑ)

            # get transformation rule
            θ_sym, transformation_rule = symmetry_transform(θ_ref)

            # get reference element domain in training range
            Ωₑ_ref = get_cut_domain(1.0, θ_sym, z_ref)

            # get parametrization
            p = Parametrization(R_ref, 1.0)

            # get reference element case
            case = parameter_case(p, θ_sym, z_ref)

            # additional check because lsq approximation is not called!
            if case == -1
                if ϕ(center(Ωₑ)) < 0
                    throw(StrictlyNegativeLevelSet())
                elseif ϕ(center(Ωₑ)) > 0
                    throw(StrictlyPositiveLevelSet())
                else
                    error("case estimation failed unexpectedly")
                end
            end

            # get inverse parameterization
            p⁻¹ = inv(p)

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
            J = QuadLearnData.measure(Ωₑ)

            # scale moments
            M_scaled = M_transformed * J

            # test complement condition
            complement = test_complement_condition(ϕ, Ωₑ, R, a, b; phase=phase)

            # evaluate complement if needed
            M = complement ? integrate_bernstein_basis(order, Ωₑ) .- M_scaled : M_scaled

            # compute linear map from weights to moments (memoized!)
            A = get_bernsteinfuns_at_gausspoints(; order=order)

            # compute test weights
            w = A \ M[:]
            
            # add to area
            area += sum(w)
        catch e
            if isa(e, StrictlyPositiveLevelSet)
                area += (phase == 1)  ? QuadLearnData.measure(Ωₑ) : 0.0
            elseif isa(e, StrictlyNegativeLevelSet)
                area += (phase == -1) ? QuadLearnData.measure(Ωₑ) : 0.0
            else
                rethrow(e)
            end
        end
    end

    # error in integration only due to model accuracy
    @test isapprox(area, (width^2 - π*radius^2), atol=10e-9)
end

@testset "Test pipeline for area integration on approximate circular cut" begin
    # circular cut radius
    radius = 2.0

    # embedding box width
    width = 4.0

    # levelset function for circle at origin with R = radius
    ϕ = implicit_circle_def(radius)

    # number of elements in each parametric direction
    nelem = 32

    # partition of [-width/2,width/2] × [-width/2,width/2]
    partition = IncreasingRange(-width/2, width/2, nelem) ⨱ IncreasingRange(-width/2, width/2, nelem)

    # model loader
    loader = init_cached_model_loader(; maxsize=75)

    # order of quadrature
    order = 2

    # order of interpolation model
    interp_order = 4

    # phase to integrate (1: outside, -1: inside)
    phase = 1

    # total area
    area = 0.0

    # element loop
    for element in Elements(partition)
        # get element domain
        Ωₑ = get_element_domain(element)

        try
            # least square circular approximation of the cut
            R, a, b = levelset_least_squares_circle(ϕ, Ωₑ; n=32)

            # compute reference cut parameters
            R_ref, θ_ref, z_ref = get_reference_cut_parameters(R, a, b, Ωₑ)

            # get transformation rule
            θ_sym, transformation_rule = symmetry_transform(θ_ref)

            # get reference element domain in training range
            Ωₑ_ref = get_cut_domain(1.0, θ_sym, z_ref)

            # get parametrization
            p = Parametrization(R_ref, 1.0)

            # get inverse parameterization
            p⁻¹ = inv(p)

            # get reference element case
            case = parameter_case(p, θ_sym, z_ref)

            # case can never be case == -1...
            # Least squares will throw an exception if element is not cut
            # and indicated in which phase the element is contained

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
            J = QuadLearnData.measure(Ωₑ)

            # scale moments
            M_scaled = M_transformed * J

            # test complement condition
            complement = test_complement_condition(ϕ, Ωₑ, R, a, b; phase=phase)

            # evaluate complement if needed
            M = complement ? integrate_bernstein_basis(order, Ωₑ) .- M_scaled : M_scaled

            # compute linear map from weights to moments (memoized!)
            A = get_bernsteinfuns_at_gausspoints(; order=order)

            # compute test weights
            w = A \ M[:]
            
            # add to area
            area += sum(w)
        catch e
            if isa(e, StrictlyPositiveLevelSet)
                area += (phase == 1)  ? QuadLearnData.measure(Ωₑ) : 0.0
            elseif isa(e, StrictlyNegativeLevelSet)
                area += (phase == -1) ? QuadLearnData.measure(Ωₑ) : 0.0
            else
                rethrow(e)
            end
        end
    end

    # error in integration due to model accuracy AND cut contour approximation
    @test isapprox(area, (width^2 - π*radius^2), atol=10e-6)
end

@testset "Get weights element by element" begin
    # circular cut radius
    radius = 2.0

    # embedding box width
    width = 4.0

    # levelset function for circle at origin with R = radius
    ϕ = implicit_circle_def(radius)

    # number of elements in each parametric direction
    nelem = 32

    # partition of [-width/2,width/2] × [-width/2,width/2]
    partition = IncreasingRange(-width/2, width/2, nelem) ⨱ IncreasingRange(-width/2, width/2, nelem)

    # model loader
    loader = init_cached_model_loader(; maxsize=75)

    # order of quadrature
    order = 2

    # order of interpolation model
    interp_order = 4

    # phase to integrate (1: outside, -1: inside)
    phase = 1

    # total area
    area = 0.0

    # element loop
    for element in Elements(partition)
        # get weights
        w = get_weights(ϕ, element; loader=loader, order=order, interp_order=interp_order, phase=phase)

        # add to area
        area += sum(w)
    end

    # error in integration due to model accuracy AND cut contour approximation
    @test isapprox(area, (width^2 - π*radius^2), atol=10e-6)
end

@testset "Integration using trained models" begin
    # levelset (sinusoid squircle, (phase=-1) == inside)
    ϕ = Algoim.AlgoimCallLevelSetFunction(
        (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*x[1]/2) - sin(3*x[2]/2) - 1.5  ), 
        (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
    )

    # Bernstein basis order
    order = 2
    
    # interpolation order
    interp_order = 4

    # some cut element
    element = get_cut_domain(0.5, deg2rad(202), 2.54)

    # integrate exact momements (phase = -1)
    M_exact = integrate_bernstein_basis(order, ϕ, element)

    # least square circular approximation of the cut
    R, a, b = levelset_least_squares_circle(ϕ, element; n=512)

    # compute reference cut parameters (perform translation and scaling of lsq approximation)
    R_ref, θ_ref, z_ref = get_reference_cut_parameters(R, a, b, element)

    # get transformation rule (for evaluation in training range)
    θ_sym, transformation_rule = symmetry_transform(θ_ref)

    # get reference element in training range
    reference_element = get_cut_domain(1.0, θ_sym, z_ref)

    # compute test moments on reference cut and element
    M_reference_test = integrate_bernstein_basis(order, implicit_circle_def(R_ref), reference_element; n=2)

    # get parametrization
    p = Parametrization(R_ref, 1.0)

    # get inverse parameterization
    p⁻¹ = inv(p)

    # get reference element case
    case = parameter_case(p, θ_sym, z_ref)

    # compute reference element (θ̂, ẑ) normalized parameter pair
    θ̂_ref, ẑ_ref = p⁻¹(θ_sym, z_ref; case=case)

    # initialize lazy loader
    loader = init_cached_model_loader(; maxsize=1)

    # get model
    model = loader(floor(Int, R_ref), order, interp_order; case=case)

    # evaluate model
    M_reference = reshape(model(R_ref, θ̂_ref, ẑ_ref), order+1, order+1)

    # test against exactly integrated circular cut (without scaling)
    @test isapprox(M_reference, M_reference_test, atol=10e-10)

    # apply transformation
    M_transformed = transformation_rule(M_reference)

    # compute Jacobian
    J = QuadLearnData.measure(element)

    # scale moments
    M_scaled = M_transformed * J

    # compute moments on uncut element
    M_uncut = integrate_bernstein_basis(order, element)

    # evaluate complement if needed
    M_test = test_complement_condition(ϕ, element, R, a, b; phase=-1) ? M_uncut .- M_scaled : M_scaled

    # test complement condition
    @test test_complement_condition(ϕ, element, R, a, b; phase=-1) == true

    # test integration of moments on the approximate circular cut
    @test isapprox(M_exact, M_test, atol=10e-7)

    # compute linear map from weights to moments
    A = get_bernsteinfuns_at_gausspoints(; order=order)

    # compute exact weights
    w_exact = A \ M_exact[:]

    # compute test weights
    w_test = A \ M_test[:]

    # test weights
    @test isapprox(w_exact, w_test, atol=10e-6)

    # define quadrature rule
    x = get_stencil(element; order=order+1)
    w = reshape(w_test, order+1, order+1)
    Q = QuadratureRule(x, w)

    # integrate some polynomial on element
    f = ScalarFunction((x,y) -> 3*x^2 + 6*y^2)
    @evaluate y = f(x)
    ∫ = 0.0
    for k in eachindex(Q.x)
        ∫ += Q.w[k] * y[k]
    end

    # integrate polynomial on element using Algoim 
    f_exact = x -> f(x[1], x[2])
    Q_exact = ImmersedQuadRule(ϕ, element[1,1], element[2,2]; order=12, phase=-1)
    ∫_exact = sum(Q_exact.w .* f_exact.(Q_exact.x))

    # test integration of polynomial on element
    @test isapprox(∫, ∫_exact, atol=10-6)
end

end