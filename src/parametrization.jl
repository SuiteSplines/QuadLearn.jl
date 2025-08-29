export Parametrization, InverseParametrization, parameter_case_test, parameter_case
export symmetry_case, symmetry_case_test, symmetry_transform
export SingularTheta
export levelset_least_squares_radius
export levelset_least_squares_circle
export test_levelset_sign_change
export fit_constrained_circle
export get_reference_cut_parameters
export test_complement_condition
export levelset_contour_coordinates, levelset_contours
export StrictlyNegativeLevelSet, StrictlyPositiveLevelSet
export NotSimplyConnectedCutDomain

struct StrictlyNegativeLevelSet <: Exception end
struct StrictlyPositiveLevelSet <: Exception end
struct NotSimplyConnectedCutDomain <: Exception end 

struct SingularTheta{T} <: Exception
    ẑ::T
end
Base.showerror(io::IO, e::SingularTheta) = print(io, "θ is almost arbitrary
for ẑ=$(e.ẑ). The value of ẑ is stored in this exception.")

struct Parametrization{T}
    R::T
    w::T
    z₁::T
    z₂::T
    z₄::T
    z₅::T
    z₆::T
    z₇::T
    z₈::T
    z₉::T
    zmax ::T
    function Parametrization(R::T, w::T) where {T<:Real}
        # R must be larger than w
        @assert R ≥ w "Parametrization is undefined for R < w (R=$R, w=$w)" 

        # angle at which for z = z₁ the element corner E₁
        # just touches the circular cut at exactly at (R,0)
        θ₁ = atan(w / (2 * (R + w / 2)))

        # angle at which for z = z₆ the element corner E₂
        # just touches the circular cut at exactly at (R,0)
        θ₅ = atan(w / (2 * (R - w / 2)))

        # maximum z for which there exists a cut in θ ∈ [0, π/4]
        z₁ = w / (2 * sin(θ₁))

        # length of the element diagonal
        d = sqrt(2) * w

        # maximum z for which there is no cut for any θ ∈ [0, π/4]
        zmax = R + d / 2

        # characteristic z values defining smooth cut regions for range θ ⊂ [0, pi/4]
        z₂ = R + w / 2
        z₄ = sqrt(R^2 - w^2 / 4) + w / 2
        z₅ = sqrt(R^2 - d^2 / 4)
        z₆ = w / (2 * sin(θ₅))
        z₇ = R - w / 2
        z₈ = sqrt(R^2 - w^2 / 4) - w / 2
        z₉ = R - w / sqrt(2)

        # initialize parametrization for a fixed (R,w)
        new{T}(R, w, z₁, z₂, z₄, z₅, z₆, z₇, z₈, z₉, zmax)
    end
end

struct InverseParametrization{T}
    p::Parametrization{T}
    function InverseParametrization(R::T, w::T) where {T}
        new{T}(Parametrization(R,w))
    end
end
Base.getproperty(p::InverseParametrization, s::Symbol) = Base.getfield(Base.getfield(p, :p), s)
Base.propertynames(p::InverseParametrization) = Base.propertynames(Base.getfield(p, :p))

# given parametrization, (θ̂, ẑ) and case return (θ,z)
(p::Parametrization)(θ̂, ẑ; case::Int) = p(θ̂, ẑ, Val(case))

# given inverse parametrization, (θ, z) and case return (θ̂,ẑ)
(p::InverseParametrization)(θ, z; case::Int) = p(θ, z, Val(case))
Base.inv(p::Parametrization) = InverseParametrization(p.R, p.w)

# given parametrization and (θ,z) return case
function parameter_case(p::Parametrization{T}, θ::T, z::T) where {T}
    # perform a greedy check (on purpose):
    #   run all test and check if none or only one is true...
    #   if there is more then one positive test, that's certainly a bug!
    tests = findall(k -> parameter_case_test(p, θ, z; case = k), Base.OneTo(15))
    n = length(tests)

    @assert n <= 1 "More then one case test returned true. \
                    This indicates a bug in parameter_case_test()!"

    # return cut case index or -1 if the element is not cut
    (n == 1) ? tests[1] : -1
end

# return true if (θ,z) corresponds to case and false otherwise
function parameter_case_test(p::Parametrization{T}, θ::T, z::T; case::Int) where {T}
    try
        return parameter_case_test(p, θ, z, Val(case))
    catch e
        isa(e, DomainError) && return false
        rethrow(e)
    end
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{1}) where {T}
    R, w, z₁, zmax = p.R, p.w, p.z₁, p.zmax
    θ₀ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    return (z ≤ p.zmax) & (z ≥ z₁) & (θ ≤ π/4) & (θ ≥ θ₀)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{2}) where {T}
    R, w, z₂, z₁ = p.R, p.w, p.z₂, p.z₁
    θ₂ = acos((R + w / 2) / z)
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    return (z < z₁) & (z ≥ z₂) & (θ < θ₃) & (θ ≥ θ₂)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{3}) where {T}
    R, w, z₂, z₁ = p.R, p.w, p.z₂, p.z₁
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    return (z < z₁) & (z ≥ z₂) & (θ ≤ π/4) & (θ ≥ θ₃)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{4}) where {T}
    R, w, z₄, z₂ = p.R, p.w, p.z₄, p.z₂
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    return (z < z₂) & (z ≥ z₄) & (θ ≥ 0) & (θ < θ₃)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{5}) where {T}
    R, w, z₄, z₂ = p.R, p.w, p.z₄, p.z₂
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    return (z < z₂) & (z ≥ z₄) & (θ ≥ θ₃) & (θ ≤ π/4)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{6}) where {T}
    R, w, z₅, z₄ = p.R, p.w, p.z₅, p.z₄
    θ₄ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) + π / 4

    return (z < z₄) & (z ≥ z₅) & (θ ≥ 0) & (θ < θ₄)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{7}) where {T}
    R, w, z₅, z₄ = p.R, p.w, p.z₅, p.z₄
    θ₄ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) + π / 4

    return (z < z₄) & (z ≥ z₅) & (θ ≥ θ₄) & (θ ≤ π/4)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{8}) where {T}
    R, w, z₆, z₅ = p.R, p.w, p.z₆, p.z₅
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    return (z < z₅) & (z ≥ z₆) & (θ ≥ 0) & (θ < θ₆)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{9}) where {T}
    R, w, z₆, z₅ = p.R, p.w, p.z₆, p.z₅
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    return (z < z₅) & (z ≥ z₆) & (θ ≥ θ₆) & (θ ≤ π/4)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{10}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆
    θ₇ = acos((R - w / 2) / z)

    return (z ≥ z₇) & (z < z₆) & (θ ≥ 0) & (θ < θ₇)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{11}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆
    θ₇ = acos((R - w / 2) / z)
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    return (z ≥ z₇) & (z < z₆) & (θ ≥ θ₇) & (θ < θ₆)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{12}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    return (z ≥ z₇) & (z < z₆) & (θ ≥ θ₆) & (θ ≤ π/4)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{13}) where {T}
    R, w, z₈, z₇ = p.R, p.w, p.z₈, p.z₇
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    return (z ≥ z₈) & (z < z₇) & (θ ≥ 0) & (θ < θ₆)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{14}) where {T}
    R, w, z₈, z₇ = p.R, p.w, p.z₈, p.z₇
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    return (z ≥ z₈) & (z < z₇) & (θ ≥ θ₆) & (θ ≤ π/4)
end

function parameter_case_test(p::Parametrization{T}, θ::T, z::T, ::Val{15}) where {T}
    R, w, z₉, z₈ = p.R, p.w, p.z₉, p.z₈
    θ₈ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    return (z ≥ z₉) & (z < z₈) & (θ ≥ θ₈) & (θ ≤ π/4)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{1}) where {T}
    R, w, z₁, zmax = p.R, p.w, p.z₁, p.zmax

    z = z₁ * (1 - ẑ) + zmax * ẑ

    # in edge cases float precision kicks in... this handles the imprecision
    # still, if the argument is much smaller than 0 or much larger than 1
    # arg will be returned back and asin will throw DomainError
    arg = (z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)
    arg = (arg < 0 && abs(arg) < sqrt(eps(T)) ) ? 0.0 : arg
    arg = (arg > 1 && abs(1 - arg) < sqrt(eps(T)) ) ? 1.0 : arg

    θ₀ = asin(arg) - π / 4
    θ = θ₀ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{2}) where {T}
    R, w, z₂, z₁ = p.R, p.w, p.z₂, p.z₁

    z = z₂ * (1 - ẑ) + z₁ * ẑ
    θ₂ = acos((R + w / 2) / z)
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    θ = θ₂ * (1 - θ̂) + θ₃ * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{3}) where {T}
    R, w, z₂, z₁ = p.R, p.w, p.z₂, p.z₁

    z = z₂ * (1 - ẑ) + z₁ * ẑ
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₃ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{4}) where {T}
    R, w, z₄, z₂ = p.R, p.w, p.z₄, p.z₂

    z = z₄ * (1 - ẑ) + z₂ * ẑ
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₃ * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{5}) where {T}
    R, w, z₄, z₂ = p.R, p.w, p.z₄, p.z₂

    z = z₄ * (1 - ẑ) + z₂ * ẑ
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₃ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{6}) where {T}
    R, w, z₅, z₄ = p.R, p.w, p.z₅, p.z₄

    z = z₅ * (1 - ẑ) + z₄ * ẑ
    θ₄ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) + π / 4
    θ = θ₄ * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{7}) where {T}
    R, w, z₅, z₄ = p.R, p.w, p.z₅, p.z₄

    z = z₅ * (1 - ẑ) + z₄ * ẑ
    θ₄ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) + π / 4
    θ = θ₄ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{8}) where {T}
    R, w, z₆, z₅ = p.R, p.w, p.z₆, p.z₅

    z = z₆ * (1 - ẑ) + z₅ * ẑ
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₆ * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{9}) where {T}
    R, w, z₆, z₅ = p.R, p.w, p.z₆, p.z₅

    z = z₆ * (1 - ẑ) + z₅ * ẑ
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₆ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{10}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆

    z = z₇ * (1 - ẑ) + z₆ * ẑ
    θ₇ = acos((R - w / 2) / z)
    θ = θ₇ * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{11}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆

    z = z₇ * (1 - ẑ) + z₆ * ẑ
    θ₇ = acos((R - w / 2) / z)
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₇ * (1 - θ̂) + θ₆ * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{12}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆

    z = z₇ * (1 - ẑ) + z₆ * ẑ
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₆ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{13}) where {T}
    R, w, z₈, z₇ = p.R, p.w, p.z₈, p.z₇

    z = z₈ * (1 - ẑ) + z₇ * ẑ
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₆ * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{14}) where {T}
    R, w, z₈, z₇ = p.R, p.w, p.z₈, p.z₇

    z = z₈ * (1 - ẑ) + z₇ * ẑ
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4
    θ = θ₆ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::Parametrization{T})(θ̂::T, ẑ::T, ::Val{15}) where {T}
    R, w, z₉, z₈ = p.R, p.w, p.z₉, p.z₈

    z = z₉ * (1 - ẑ) + z₈ * ẑ

    # in edge cases float precision kicks in... this handles the imprecision
    # still, if the argument is much smaller than 0 or much larger than 1
    # arg will be returned back and asin will throw DomainError
    arg = (R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)
    arg = (arg < 0 && abs(arg) < sqrt(eps(T)) ) ? 0.0 : arg
    arg = (arg > 1 && abs(1 - arg) < sqrt(eps(T)) ) ? 1.0 : arg

    θ₈ = asin(arg) - π / 4
    θ = θ₈ * (1 - θ̂) + π/4 * θ̂

    return (θ = θ, z = z)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{1}) where {T}
    R, w, z₁, zmax = p.R, p.w, p.z₁, p.zmax
    θ₀ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₁) / (zmax - z₁)
    θ̂ = (θ - θ₀) / (π/4 - θ₀)

    (ẑ > (1 - 10e-5)) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{2}) where {T}
    R, w, z₂, z₁ = p.R, p.w, p.z₂, p.z₁
    θ₂ = acos((R + w / 2) / z)
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₂) / (z₁ - z₂)
    θ̂ = (θ - θ₂) / (θ₃ - θ₂)
    
    (ẑ > (1 - 10e-5)) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{3}) where {T}
    R, w, z₂, z₁ = p.R, p.w, p.z₂, p.z₁
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₂) / (z₁ - z₂)
    θ̂ = (θ - θ₃) / (π/4 - θ₃)

    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{4}) where {T}
    R, w, z₄, z₂ = p.R, p.w, p.z₄, p.z₂
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₄) / (z₂ - z₄)
    θ̂ = θ / θ₃

    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{5}) where {T}
    R, w, z₄, z₂ = p.R, p.w, p.z₄, p.z₂
    θ₃ = asin((z^2 + w^2 / 2 - R^2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₄) / (z₂ - z₄)
    θ̂ = (θ - θ₃) / (π/4 - θ₃)

    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{6}) where {T}
    R, w, z₅, z₄ = p.R, p.w, p.z₅, p.z₄
    θ₄ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) + π / 4

    ẑ = (z - z₅) / (z₄ - z₅)
    θ̂ = θ / θ₄

    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{7}) where {T}
    R, w, z₅, z₄ = p.R, p.w, p.z₅, p.z₄
    θ₄ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) + π / 4

    ẑ = (z - z₅) / (z₄ - z₅)
    θ̂ = (θ - θ₄) / (π/4 - θ₄)

    (ẑ < 10e-5) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{8}) where {T}
    R, w, z₆, z₅ = p.R, p.w, p.z₆, p.z₅
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₆) / (z₅ - z₆)
    θ̂ = θ / θ₆

    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{9}) where {T}
    R, w, z₆, z₅ = p.R, p.w, p.z₆, p.z₅
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₆) / (z₅ - z₆)
    θ̂ = (θ - θ₆) / (π/4 - θ₆)

    (ẑ > (1 - 10e-5)) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{10}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆
    θ₇ = acos((R - w / 2) / z)

    ẑ = (z - z₇) / (z₆ - z₇)
    θ̂ = θ / θ₇

    (ẑ < 10e-5) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{11}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆
    θ₇ = acos((R - w / 2) / z)
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₇) / (z₆ - z₇)
    θ̂ = (θ - θ₇) / (θ₆ - θ₇)

    (ẑ > (1 - 10e-5)) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{12}) where {T}
    R, w, z₇, z₆ = p.R, p.w, p.z₇, p.z₆
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₇) / (z₆ - z₇)
    θ̂ = (θ - θ₆) / (π/4 - θ₆)

    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{13}) where {T}
    R, w, z₈, z₇ = p.R, p.w, p.z₈, p.z₇
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₈) / (z₇ - z₈)
    θ̂ = θ / θ₆

    (ẑ < 10e-5) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{14}) where {T}
    R, w, z₈, z₇ = p.R, p.w, p.z₈, p.z₇
    θ₆ = acos((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₈) / (z₇ - z₈)
    θ̂ = (θ - θ₆) / (π/4 - θ₆)

    return (θ̂ = θ̂, ẑ = ẑ)
end

function (p::InverseParametrization{T})(θ::T, z::T, ::Val{15}) where {T}
    R, w, z₉, z₈ = p.R, p.w, p.z₉, p.z₈
    θ₈ = asin((R^2 - z^2 - w^2 / 2) / (sqrt(2) * z * w)) - π / 4

    ẑ = (z - z₉) / (z₈ - z₉)
    θ̂ = (θ - θ₈) / (π/4 - θ₈)

    (ẑ < 10e-5) && throw(SingularTheta(ẑ))
    return (θ̂ = θ̂, ẑ = ẑ)
end

function get_cut_domain(p::Parametrization{T}, θ̂::T, ẑ::T; case::Int) where {T}
    # evaluate parametrization
    θ, z = p(θ̂, ẑ; case = case)

    # physical element width
    w = p.w

    # unit vector in direction of translation
    v = SVector(cos(θ), sin(θ))

    # center of physical element
    C = z * v

    # return physical element
    return Interval(C[1] - w/2, C[1] + w/2) ⨱ Interval(C[2] - w/2, C[2] + w/2)
end

function get_cut_domain(w::T, θ::T, z::T) where {T}
    # unit vector in direction of translation
    v = SVector(cos(θ), sin(θ))

    # center of physical element
    C = z * v

    # return physical element
    return Interval(C[1] - w/2, C[1] + w/2) ⨱ Interval(C[2] - w/2, C[2] + w/2)
end


## symmetry handling
function symmetry_case(θ::T) where {T}
    # sanitize input s.t. θ ∈ [0, 2π)
    θ = mod(θ, 2π)

    # perform a greedy check (on purpose):
    #   run all tests and check if none or only one is true...
    #   if there is more then one positive test, that's certainly a bug!
    tests = findall(k -> symmetry_case_test(θ; case = k), Base.OneTo(8))
    n = length(tests)

    @assert n <= 1 "More then one case test returned true. \
                    This indicates a bug in symmetry_case_test()!"

    # return cut case index or -1 if the element is not cut
    (n == 1) ? tests[1] : -1
end
symmetry_case_test(θ::T; case::Int) where {T} = symmetry_case_test(θ, Val(case))
symmetry_case_test(θ::T, ::Val{1}) where {T} = ((θ ≥ 0)    && (θ ≤ π/4))  ? true : false
symmetry_case_test(θ::T, ::Val{2}) where {T} = ((θ > π/4)  && (θ ≤ π/2))  ? true : false
symmetry_case_test(θ::T, ::Val{3}) where {T} = ((θ > π/2)  && (θ ≤ 3π/4)) ? true : false
symmetry_case_test(θ::T, ::Val{4}) where {T} = ((θ > 3π/4) && (θ ≤ π))    ? true : false
symmetry_case_test(θ::T, ::Val{5}) where {T} = ((θ > π)    && (θ ≤ 5π/4)) ? true : false
symmetry_case_test(θ::T, ::Val{6}) where {T} = ((θ > 5π/4) && (θ ≤ 3π/2)) ? true : false
symmetry_case_test(θ::T, ::Val{7}) where {T} = ((θ > 3π/2) && (θ ≤ 7π/4)) ? true : false
symmetry_case_test(θ::T, ::Val{8}) where {T} = ((θ > 7π/4) && (θ < 2π))   ? true : false

reverse1(A) = reverse(A, dims=1)
reverse2(A) = reverse(A, dims=2)

function symmetry_transform(θ::T) where {T}
    # sanitize input s.t. θ ∈ [0, 2π)
    θ = mod(θ, 2π) # [0, 2π]
    θ = isapprox(θ, 2π, atol=sqrt(eps(T))) ? 0.0 : θ # [0, 2π)

    # find case
    case = symmetry_case(θ)
    @assert case ≥ 1

    # return named tuple with reference θ and transformation rule
    symmetry_transform(θ, Val(case))
end
symmetry_transform(θ::T, ::Val{1}) where {T} = (θ = θ, t = identity)
symmetry_transform(θ::T, ::Val{2}) where {T} = (θ = (π/2 - θ), t = transpose)
symmetry_transform(θ::T, ::Val{3}) where {T} = (θ = (θ - π/2), t = reverse1 ∘ transpose)
symmetry_transform(θ::T, ::Val{4}) where {T} = (θ = (π - θ), t = reverse1)
symmetry_transform(θ::T, ::Val{5}) where {T} = (θ = (θ - π), t = reverse)
symmetry_transform(θ::T, ::Val{6}) where {T} = (θ = (3π/2 - θ), t = transpose ∘ reverse)
symmetry_transform(θ::T, ::Val{7}) where {T} = (θ = (θ - 3π/2), t = reverse2 ∘ transpose)
symmetry_transform(θ::T, ::Val{8}) where {T} = (θ = (2π - θ), t = reverse2)

function levelset_contours(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}; n::Int=32, level::T=0.0) where {T}
    # x and y interval of element sides
    Ix, Iy = element.data[1], element.data[2]

    # test points
    X = IncreasingRange(Ix.a, Ix.b, n) ⨱ IncreasingRange(Iy.a, Iy.b, n)

    # evalaute levelset at each test point
    z = ϕ.(X)

    # compute contour using marching squares
    c = Contour.contour(X.data[1], X.data[2], z, level)

    # get contour lines
    l = lines(c)

    # return contour lines
    coordinates.(l)
end

function levelset_contour_coordinates(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}; n::Int=32, level::T=0.0) where {T}
    # x and y interval of element sides
    Ix, Iy = element.data[1], element.data[2]

    # test points
    X = IncreasingRange(Ix.a, Ix.b, n) ⨱ IncreasingRange(Iy.a, Iy.b, n)

    # evalaute levelset at each test point
    z = ϕ.(X)

    # compute contour using marching squares
    c = Contour.contour(X.data[1], X.data[2], z, level)

    # get contour lines
    l = lines(c)

    # if no contour is not found, element is not cut
    # thus ϕ must be either strictly positive or negative on element
    if length(l) == 0
        if all(z .> 0)
            throw(StrictlyPositiveLevelSet())
        elseif all(z .< 0)
            throw(StrictlyNegativeLevelSet())
        else
            error("levelset contour estimation failed!")
        end
    end
    
    # collect contour lines coordinates
    coords = coordinates.(l)

    # collect all x and y coordinates
    x = vcat(getindex.(coords, 1)...)
    y = vcat(getindex.(coords, 2)...)
    
    # the contour touches the boundary only at a single point
    # thus ϕ must be either strictly positive or negative on element
    if length(x) == 2
        if all(z .≥ -10e-3)
            throw(StrictlyPositiveLevelSet())
        elseif all(z .≤ 10e-3)
            throw(StrictlyNegativeLevelSet())
        else
            error("levelset contour estimation failed!")
        end
    end

    return x, y
end

function levelset_least_squares_radius(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}; n=32)
    # compute zero contour line coordinates using marching squares
    x, y = levelset_contour_coordinates(ϕ, element; n=n)

    # compute least squares R = ∑ᵢ ( √(xᵢ^2 + yᵢ^2) ) / N
    R = sum(norm.(eachrow([x y]))) / length(x)
end

fit_circle_objective(x::Vector{T}, y::Vector{T}, R::T, a::T, b::T) where {T} = sum(((x .- a).^2 .+ (y .- b).^2 .- R^2).^2)

function fit_constrained_circle(x::Vector{T}, y::Vector{T}; R::T) where {T}
    # initial guess
    a₀, b₀ = mean(x), mean(y)

    # bounds
    lower_bounds = [-Inf, -Inf]
    upper_bounds = [Inf,  Inf]

    # optimize
    opt = optimize(
        p -> fit_circle_objective(x, y, R, p[1], p[2]),
        lower_bounds,
        upper_bounds,
        [a₀, b₀],
        Fminbox(),
    )
    a, b = Optim.minimizer(opt)

    return (R=R, a=a, b=b)

    # TODO: error handling / residual check / ...
end

function levelset_least_squares_circle(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}; n=32, constrain::Bool = false, linear_max_constrain::Bool = false, Rmax::T=30.0, Rmin::T=1.0) where {T}
    # compute zero contour line coordinates using marching squares
    x, y = levelset_contour_coordinates(ϕ, element; n=n)
    ncoords = length(x)

    # solve least squares for circle (x - a)² + (y - b)² = R² that best fits the curve described by ϕ = 0
    A = [x y ones(ncoords)]
    r = -(x.^2 + y.^2)
    lsq = A \ r
    
    # recover parameters
    a = -lsq[1] / 2
    b = -lsq[2] / 2
    R = sqrt(a^2 + b^2 - lsq[3])

    if linear_max_constrain
        # width and height of element
        w, h = measure.(element.data)

        # elements must be square
        @assert isapprox(w, h, atol=10e-13) "element must be square"

        # compute scaling factor
        s = 1 / w

        # reference radius
        Rref = R * s

        # constrained circle if necessary
        if Rref > Rmax
            # two points on circle contour
            p₁ = SVector(first(x), first(y))
            p₂ = SVector(last(x), last(y))

            # circle through p₁ and p₂ that approximates too large (R,a,b)
            v = p₂ - p₁
            m = 0.5 * (p₁ + p₂)
            d = norm(v)

            h = sqrt((w*Rmax)^2 - d^2/4)
            u = 1/d * v
            n = SVector(-u[2], u[1])
            O = m + h * n

            # this is not *the* best fit but very close for large training Rmax
            # in addition, if the cut is a curve between p₁ and p₂ and interesects
            # the cell boundaries only twice then the contour will be continouous
            # an adjacent cells
            return (R=Rmax, a=O[1], b=O[2])
        end
    end


    # enforce constraints (nonlinear optimization)
    # this will become unnecessary as soon as
    # straight cut are implemented
    if constrain
        # width and height of element
        w, h = measure.(element.data)

        # elements must be square
        @assert isapprox(w, h, atol=10e-13) "element must be square"

        # compute scaling factor
        s = 1 / w

        # reference radius
        Rref = R * s

        # perfrom constrained optimization if necessary
        if Rref > Rmax
            #@warn "Running constrained optimization"
            R, a, b = fit_constrained_circle(x, y; R=w*Rmax)
        elseif Rref < Rmin
            #@warn "Running constrained optimization"
            R, a, b = fit_constrained_circle(x, y; R=w*Rmin)
        end
    end


    return (R=R, a=a, b=b)
end

function get_reference_cut_parameters(R::T, a::T, b::T, element::CartesianProduct{2}; unitscale::Bool=true) where {T}
    # width and height of element
    w, h = measure.(element.data)

    # elements must be square
    @assert isapprox(w, h, atol=10e-13) "element must be square"

    # compute scaling factor
    s = unitscale ? 1.0 / w : 1.0

    # center of physical element
    C = center(element)

    # center of reference element
    C_ref = C .- (a, b)

    # angle parameter for reference element
    θ_ref = angle(complex(C_ref...))

    # translation parameter for reference element
    z_ref = norm(C_ref)

    # return parameter triplet
    return (R = s * R, θ = θ_ref, z = s * z_ref)
end

function is_point_inside_circle(x::Tuple{T,T}, R::T, a::T, b::T) where {T}
    # vector from circle center to point x
    Δ = x .- (a, b)

    # true if length of Δ is smaller or equal to circle radius
    norm(Δ) ≤ R
end

test_complement_condition(ϕ_point::T, point_inside_circle::Val{true},  phase::Int) where {T} = sign(ϕ_point) == phase ? false : true
test_complement_condition(ϕ_point::T, point_inside_circle::Val{false}, phase::Int) where {T} = sign(ϕ_point) == phase ? true : false 
function test_complement_condition(ϕ::AlgoimCallLevelSetFunction, element::CartesianProduct{2}, R::T, a::T, b::T; phase::Int=-1, n::Int=8) where {T}
    # Note: this test logic is reasonably robust but might fail for
    # cases for which the cut contour is not approximated well enough
    # by the circle. These cases are nonsense anyways because the
    # quadrature rule would be garbage too...

    # x and y interval of element sides
    Ix, Iy = element.data[1], element.data[2]

    # test points
    x = IncreasingRange(Ix.a, Ix.b, n) ⨱ IncreasingRange(Iy.a, Iy.b, n)

    # evalaute levelset at each test point
    z = ϕ.(x)

    # index of test point for which |ϕ| is maximal (cut accuracy safeguard!)
    k = argmax(abs.(z))

    # test if test point is inside circle
    inside_circle = is_point_inside_circle(x[k], R, a, b)

    # test complement condition and return
    test_complement_condition(z[k], Val(inside_circle), phase)
end