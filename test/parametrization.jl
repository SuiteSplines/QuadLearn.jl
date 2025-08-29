module TestParametrization

using QuadLearnData, Test, Revise, Algoim

@testset "Initialization" begin
    p = Parametrization(4.75, 2.5)

    @test p.R == 4.75
    @test p.w == 2.5
    @test isapprox(p.z₁, 6.128825336065631, atol=1e-12)
    @test isapprox(p.z₂, 6.0, atol=1e-11)
    @test isapprox(p.z₄, 5.83257569495584, atol=1e-12)
    @test isapprox(p.z₅, 4.408798022137099, atol=1e-12)
    @test isapprox(p.z₆, 3.7165171868296265, atol=1e-12)
    @test isapprox(p.z₇, 3.5, atol=1e-12)
    @test isapprox(p.z₈, 3.33257569495584, atol=1e-12)
    @test isapprox(p.z₉, 2.982233047033631, atol=1e-12)
    @test isapprox(p.zmax, 6.517766952966369, atol=1e-12)
end

@testset "Parametrization (R,θ̂,ẑ; case) ↔ (R,θ,z)" begin
    p = Parametrization(4.75, 2.5)
    p⁻¹ = inv(p)
    parameter_range = LinRange(0, 1, 100)
    for case = 1:15
        for ẑ = parameter_range
            for θ̂ in parameter_range
                θ, z = p(θ̂, ẑ; case = case)
                try
                    θ̂_test, ẑ_test = p⁻¹(θ, z; case = case) 
                    @test isapprox(θ̂, θ̂_test, atol=10e-8)
                    @test isapprox(ẑ, ẑ_test, atol=10e-8)
                catch e_test 
                    if isa(e_test, SingularTheta) # limit cases where θ̂ is nearly arbitrary
                        @test isapprox(ẑ, e_test.ẑ, atol=10e-8) # ẑ is either close to 0 or 1
                    end
                end
            end
        end
    end
end

@testset "Parameter case identification and inverse parametrization" begin
    p = Parametrization(4.75, 2.5)
    p⁻¹ = inv(p)

    # test at midpoints of θ̂ and ẑ interval
    θ̂, ẑ = 0.5, 0.5
    for case = 1:15
        θ, z = p(θ̂, ẑ; case = case)
        case_test = parameter_case(p, θ, z)
        @test case == case_test
    end

    # test at the beginning of θ̂ and ẑ interval
    θ̂, ẑ = 0.0, 0.0
    @test parameter_case(p, p(θ̂, ẑ; case = 1)...) == 1
    @test parameter_case(p, p(θ̂, ẑ; case = 2)...) == 2
    @test parameter_case(p, p(θ̂, ẑ; case = 3)...) == 3
    @test parameter_case(p, p(θ̂, ẑ; case = 4)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 5)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 6)...) == 6
    @test parameter_case(p, p(θ̂, ẑ; case = 7)...) == 7
    @test parameter_case(p, p(θ̂, ẑ; case = 8)...) == 8
    @test parameter_case(p, p(θ̂, ẑ; case = 9)...) == 9
    @test parameter_case(p, p(θ̂, ẑ; case = 10)...) == 11
    @test parameter_case(p, p(θ̂, ẑ; case = 11)...) == 11
    @test parameter_case(p, p(θ̂, ẑ; case = 12)...) == 12
    @test parameter_case(p, p(θ̂, ẑ; case = 13)...) == 14
    @test parameter_case(p, p(θ̂, ẑ; case = 14)...) == 14
    @test parameter_case(p, p(θ̂, ẑ; case = 15)...) == 15


    # test at the end of θ̂ and ẑ interval
    θ̂, ẑ = 1.0, 1.0
    @test parameter_case(p, p(θ̂, ẑ; case = 1)...) == 1
    @test parameter_case(p, p(θ̂, ẑ; case = 2)...) == 1
    @test parameter_case(p, p(θ̂, ẑ; case = 3)...) == 1
    @test parameter_case(p, p(θ̂, ẑ; case = 4)...) == 3
    @test parameter_case(p, p(θ̂, ẑ; case = 5)...) == 3
    @test parameter_case(p, p(θ̂, ẑ; case = 6)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 7)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 8)...) == 7
    @test parameter_case(p, p(θ̂, ẑ; case = 9)...) == 7
    @test parameter_case(p, p(θ̂, ẑ; case = 10)...) == 8
    @test parameter_case(p, p(θ̂, ẑ; case = 11)...) == 9
    @test parameter_case(p, p(θ̂, ẑ; case = 12)...) == 9
    @test parameter_case(p, p(θ̂, ẑ; case = 13)...) == 12
    @test parameter_case(p, p(θ̂, ẑ; case = 14)...) == 12
    @test parameter_case(p, p(θ̂, ẑ; case = 15)...) == 14

    # test at the beginning of ẑ inverval and the end of θ̂ interval
    θ̂, ẑ = 1.0, 0.0
    @test parameter_case(p, p(θ̂, ẑ; case = 1)...) == 1
    @test parameter_case(p, p(θ̂, ẑ; case = 2)...) == 3
    @test parameter_case(p, p(θ̂, ẑ; case = 3)...) == 3
    @test parameter_case(p, p(θ̂, ẑ; case = 4)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 5)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 6)...) == 7
    @test parameter_case(p, p(θ̂, ẑ; case = 7)...) == 7
    @test parameter_case(p, p(θ̂, ẑ; case = 8)...) == 9
    @test parameter_case(p, p(θ̂, ẑ; case = 9)...) == 9
    @test parameter_case(p, p(θ̂, ẑ; case = 10)...) == 11
    @test parameter_case(p, p(θ̂, ẑ; case = 11)...) == 12
    @test parameter_case(p, p(θ̂, ẑ; case = 12)...) == 12
    @test parameter_case(p, p(θ̂, ẑ; case = 13)...) == 14
    @test parameter_case(p, p(θ̂, ẑ; case = 14)...) == 14
    @test parameter_case(p, p(θ̂, ẑ; case = 15)...) == 15

    # test at the beginning of θ̂ inverval and the end of ẑ interval
    θ̂, ẑ = 0.0, 1.0
    @test parameter_case(p, p(θ̂, ẑ; case = 1)...) == 1
    @test parameter_case(p, p(θ̂, ẑ; case = 2)...) == -1
    @test parameter_case(p, p(θ̂, ẑ; case = 3)...) == 1
    @test parameter_case(p, p(θ̂, ẑ; case = 4)...) == 2
    @test parameter_case(p, p(θ̂, ẑ; case = 5)...) == 3
    @test parameter_case(p, p(θ̂, ẑ; case = 6)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 7)...) == 5
    @test parameter_case(p, p(θ̂, ẑ; case = 8)...) == 6
    @test parameter_case(p, p(θ̂, ẑ; case = 9)...) == 7
    @test parameter_case(p, p(θ̂, ẑ; case = 10)...) == 8
    @test parameter_case(p, p(θ̂, ẑ; case = 11)...) == 8
    @test parameter_case(p, p(θ̂, ẑ; case = 12)...) == 9
    @test parameter_case(p, p(θ̂, ẑ; case = 13)...) == 11
    @test parameter_case(p, p(θ̂, ẑ; case = 14)...) == 12
    @test parameter_case(p, p(θ̂, ẑ; case = 15)...) == 14

    # greedy testing over the parameter space
    for θ = LinRange(0, π/4, 100)
        for z = LinRange(p.z₉, p.zmax, 100)
            case = parameter_case(p, θ, z)
            if case !== -1 # cut exists
                try
                    θ̂_test, ẑ_test = p⁻¹(θ, z; case = case)
                    θ_test, z_test = p(θ̂_test, ẑ_test; case = case)
                    @test isapprox(θ, θ_test, atol=1e-12)
                    @test isapprox(z, z_test, atol=1e-12)
                catch e_test 
                    if isa(e_test, SingularTheta) # θ is nearly arbitrary (singular)
                        θ_test, z_test = p(0.5, e_test.ẑ; case=case)
                        @test isapprox(θ, θ_test, atol=1e-12)
                        @test isapprox(z, z_test, atol=1e-12)
                    else
                        rethrow(e_test)
                    end
                end
            end
        end
    end
end

@testset "Symmetry transformations" begin
    # fixed parameters
    R, w, order = 2.5, 1.6, 2
    p = Parametrization(R, w)

    for z in LinRange(1.75, 3.25, 100) # inside → cut → outside
        for ϕ in 0:π/8:5π # angle ϕ inside and outside of training range [0, π/4]
            # element domain at angle ϕ
            element_ϕ = get_cut_domain(w, ϕ, z)

            # angle θ in parameter range [0, π/4] with analogous cut to the cut at ϕ 
            θ, transformation_rule = symmetry_transform(ϕ)

            # element domain in training parameter space
            element_θ = get_cut_domain(w, θ, z)

            # compute moments on element in training parameter space
            moments_θ = integrate_bernstein_basis(order, implicit_circle_def(R), element_θ; n=2)

            # compute test moments on element outside of training parameter space
            moments_ϕ = integrate_bernstein_basis(order, implicit_circle_def(R), element_ϕ; n=2)

            # test sum of moments
            @test isapprox(sum(moments_θ), sum(moments_ϕ), atol=10e-14)

            # test application of transformation rule
            @test isapprox(transformation_rule(moments_θ), moments_ϕ, atol=10e-14)
        end
    end
end

@testset "LevelSet contour coordinates" begin
    # levelset (sinusoid squircle, (phase=-1) == inside)
    ϕ = Algoim.AlgoimCallLevelSetFunction(
        (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*x[1]/2) - sin(3*x[2]/2) - 1.5  ), 
        (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
    )

    # some cut element fully inside, outside or cut
    element_inside = get_cut_domain(0.25, deg2rad(330), 2.0)
    element_outside = get_cut_domain(0.25, deg2rad(330), 3.0)
    element_cut = get_cut_domain(0.25, deg2rad(330), 2.5)

    # test if negative levelset is thrown
    @test_throws StrictlyNegativeLevelSet levelset_contour_coordinates(ϕ, element_inside)

    # test if positive levelset is thrown
    @test_throws StrictlyPositiveLevelSet levelset_contour_coordinates(ϕ, element_outside)

    # test if contour coordinates are returned
    x, y = levelset_contour_coordinates(ϕ, element_cut; n=32)
    for k in eachindex(x)
        X = (x[k], y[k])
        @test isapprox(ϕ(X), 0.0, atol=10e-5)
    end
end

@testset "Levelset least squares radius approximation" begin
    # levelset (off-center ellipsoid cut)
    ϕ = Algoim.AlgoimCallLevelSetFunction(
        (x) -> ((x[1]-3)/2)^2 + ((x[2]+1)/3)^2 - 6, 
        (x) -> [0.5 * (x[1] - 3), 2/9 * (x[2] + 1)]
    )

    # some elements which are in fact cut
    element₁ = get_cut_domain(1.0, deg2rad(30), 8.0) 
    element₂ = get_cut_domain(1.0, deg2rad(0), 7.5) 

    # radius of a circle centered at the origin which approximates the cut
    R₁ = levelset_least_squares_radius(ϕ, element₁; n=32)
    R₂ = levelset_least_squares_radius(ϕ, element₂; n=32)

    # test against visually obtained R_test
    @test isapprox(R₁, 7.71, atol=10e-3)
    @test isapprox(R₂, 7.85, atol=10e-3)

    # test moments for approximate circular cut on element₁ and element₂
    M₁ = integrate_bernstein_basis(4, implicit_circle_def(R₁), element₁)
    M₁_test = integrate_bernstein_basis(4, ϕ, element₁)

    M₂ = integrate_bernstein_basis(4, implicit_circle_def(R₂), element₂)
    M₂_test = integrate_bernstein_basis(4, ϕ, element₂)

    @test isapprox(M₁, M₁_test, atol=10e-4)
    @test isapprox(M₂, M₂_test, atol=10e-3)

    # test if for a smaller element the approximation improves
    @test !isapprox(M₁, M₁_test, atol=10e-5)
    element₁ = get_cut_domain(0.25, deg2rad(30), 7.7) 
    R₁ = levelset_least_squares_radius(ϕ, element₁; n=32)
    M₁ = integrate_bernstein_basis(4, implicit_circle_def(R₁), element₁)
    M₁_test = integrate_bernstein_basis(4, ϕ, element₁)
    @test isapprox(M₁, M₁_test, atol=10e-5)
end

@testset "Levelset least squares circle approximation" begin
    # levelset (sinusoid squircle)
    ϕ = Algoim.AlgoimCallLevelSetFunction(
        (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*x[1]/2) - sin(3*x[2]/2) - 1.5  ), 
        (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
    )

    # some elements which are in fact cut
    element₁ = get_cut_domain(1.0, deg2rad(44), 2.8) 
    element₂ = get_cut_domain(1.0, deg2rad(164), 2.8) 

    # circle with radius R centered at (a,b) that approximates the cut
    R₁, a₁, b₁ = levelset_least_squares_circle(ϕ, element₁; n=32)
    R₂, a₂, b₂ = levelset_least_squares_circle(ϕ, element₂; n=32)

    # test moments for approximate circular cut on element₁ and element₂
    M₁ = integrate_bernstein_basis(4, implicit_circle_def(R₁, a₁, b₁), element₁)
    M₁_test = integrate_bernstein_basis(4, ϕ, element₁)

    M₂ = integrate_bernstein_basis(4, implicit_circle_def(R₂, a₂, b₂), element₂)
    M₂_test = integrate_bernstein_basis(4, ϕ, element₂)

    @test isapprox(M₁, M₁_test, atol=10e-4)
    @test isapprox(M₁, M₁_test, atol=10e-4)

    # test if for a smaller element the approximation improves
    element₁ = get_cut_domain(0.25, deg2rad(44), 2.9) 
    element₂ = get_cut_domain(0.25, deg2rad(164), 2.8) 
    R₁, a₁, b₁ = levelset_least_squares_circle(ϕ, element₁; n=32)
    R₂, a₂, b₂ = levelset_least_squares_circle(ϕ, element₂; n=32)

    M₁ = integrate_bernstein_basis(4, implicit_circle_def(R₁, a₁, b₁), element₁)
    M₁_test = integrate_bernstein_basis(4, ϕ, element₁)

    M₂ = integrate_bernstein_basis(4, implicit_circle_def(R₂, a₂, b₂), element₂)
    M₂_test = integrate_bernstein_basis(4, ϕ, element₂)

    @test isapprox(M₁, M₁_test, atol=10e-7)
    @test isapprox(M₂, M₂_test, atol=10e-7)
end

@testset "Reference cut parameters and transformation" begin
    # levelset (sinusoid squircle)
    ϕ = Algoim.AlgoimCallLevelSetFunction(
        (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*x[1]/2) - sin(3*x[2]/2) - 1.5  ), 
        (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
    )

    # some element which is in fact cut and not unit square
    element = get_cut_domain(0.8, deg2rad(44), 2.8)

    # integrate exact momements
    M_exact = integrate_bernstein_basis(4, ϕ, element)

    # least square circular approximation of the cut
    R, a, b = levelset_least_squares_circle(ϕ, element; n=32)

    # integrate momements on approximated cut
    M = integrate_bernstein_basis(4, implicit_circle_def(R, a, b), element)

    # compute reference cut parameters
    R_ref, θ_ref, z_ref = get_reference_cut_parameters(R, a, b, element)

    # get transformation rule
    θ, transformation_rule = symmetry_transform(θ_ref)

    # get reference element in training range
    reference_element = get_cut_domain(1.0, θ, z_ref)

    # compute test moments on reference cut and element
    M_ref = integrate_bernstein_basis(4, implicit_circle_def(R_ref), reference_element; n=2)

    # apply transformation
    M_test = transformation_rule(M_ref)

    # compute Jacobian
    J = measure(element)

    # test integration of moments
    @test isapprox(M, (M_test * J), atol=10e-14)
    @test isapprox(M_exact, (M_test * J), atol=10e-5)
end

@testset "Complement condition test" begin
    # levelset (sinusoid squircle)
    ϕ = Algoim.AlgoimCallLevelSetFunction(
        (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*x[1]/2) - sin(3*x[2]/2) - 1.5  ), 
        (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
    )

    element = get_cut_domain(0.8, deg2rad(44), 2.8)
    R, a, b = levelset_least_squares_circle(ϕ, element; n=32)
    @test test_complement_condition(ϕ, element, R, a, b; phase=1) == true
    @test test_complement_condition(ϕ, element, R, a, b; phase=-1) == false

    element = get_cut_domain(0.8, deg2rad(90), 2.0)
    R, a, b = levelset_least_squares_circle(ϕ, element; n=32)
    @test test_complement_condition(ϕ, element, R, a, b; phase=1) == false
    @test test_complement_condition(ϕ, element, R, a, b; phase=-1) == true

    element = get_cut_domain(0.8, deg2rad(160), 2.8)
    R, a, b = levelset_least_squares_circle(ϕ, element; n=32)
    @test test_complement_condition(ϕ, element, R, a, b; phase=1) == true
    @test test_complement_condition(ϕ, element, R, a, b; phase=-1) == false

    element = get_cut_domain(0.8, deg2rad(270), 0.5)
    R, a, b = levelset_least_squares_circle(ϕ, element; n=32)
    @test test_complement_condition(ϕ, element, R, a, b; phase=1) == false
    @test test_complement_condition(ϕ, element, R, a, b; phase=-1) == true
end


@testset "Circular cut approximation and integration" begin
    # levelset (sinusoid squircle, (phase=-1) == inside)
    ϕ = Algoim.AlgoimCallLevelSetFunction(
        (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*x[1]/2) - sin(3*x[2]/2) - 1.5  ), 
        (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
    )

    # some cut element
    element = get_cut_domain(0.25, deg2rad(239), 0.88)

    # integrate exact momements (phase = -1)
    M_exact = integrate_bernstein_basis(4, ϕ, element)

    # least square circular approximation of the cut (falls outside of phase=-1)
    R, a, b = levelset_least_squares_circle(ϕ, element; n=512)

    # compute reference cut parameters
    R_ref, θ_ref, z_ref = get_reference_cut_parameters(R, a, b, element)

    # get transformation rule
    θ, transformation_rule = symmetry_transform(θ_ref)

    # get reference element in training range
    reference_element = get_cut_domain(1.0, θ, z_ref)

    # compute test moments on reference cut and element
    M_reference = integrate_bernstein_basis(4, implicit_circle_def(R_ref), reference_element; n=2)

    # apply transformation
    M_transformed = transformation_rule(M_reference)

    # compute Jacobian
    J = measure(element)

    # compute moments on uncut element
    M_uncut = integrate_bernstein_basis(4, element)

    # test that the complement is actually needed
    @test test_complement_condition(ϕ, element, R, a, b; phase=-1) == true

    # test integration of moments on the complement of circular cut
    @test isapprox(M_exact, (M_uncut .- M_transformed * J), atol=10e-7)
end

end # module TestParametrization