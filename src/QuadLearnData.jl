module QuadLearnData

using LinearAlgebra, StaticArrays, StatsBase
import Statistics: minmax, middle

using Algoim, ImmersedSplines, Contour
using UnivariateSplines, CartesianProducts, KroneckerProducts, IgaFormation, SortedSequences
using SpecialSpaces

using Optim

using IgaBase: QuadratureRule

using ProgressMeter

using LFUDACache, Memoization, FileIO

import AbstractMappings: Interval

export implicit_circle_def
export Interval
export CartesianProducts, ⨱

implicit_circle_def(r::Real) = Algoim.AlgoimCallLevelSetFunction(
                (x) -> x[1]*x[1] + x[2]*x[2] - r^2, 
                (x) -> [2.0 * x[1], 2.0 * x[2]]
            )

implicit_circle_def(r::Real, a::Real, b::Real) = Algoim.AlgoimCallLevelSetFunction(
                (x) -> (x[1] - a)^2 + (x[2] - b)^2 - r^2, 
                (x) -> [2.0 * (x[1] - a), 2.0 * (x[2] - b)]
            )

include("utils.jl")
include("bernstein.jl")
include("parametrization.jl")
include("quadrature.jl")
include("trainingdata.jl")
include("models.jl")
include("adaptivity.jl")

end # module 