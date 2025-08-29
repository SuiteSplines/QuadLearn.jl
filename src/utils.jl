# given a uniform knot vector [0,1] generate knot vectors
# s.t. at intervals (l,r) the knots are clustered at 
# the left or the right edge.
#
# The clustering is defined by (p,l,r) = (Int, Bool, Bool)

export refined_uniform_unit_linear_range
export map_refined_range
export refined_increasing_vector
export measure, center

function refined_uniform_unit_linear_range(p_l, p_r, n)
    v = LinRange(0.0, 1.0, n)
    w = map(x -> x^p_l / (x^p_l + (1 - x)^p_r), v)
end

function map_refined_range(w, a, b)
    @assert b > a
    d = b - a
    v = (d .* w) .+ a
end

function refined_increasing_vector(; n::NTuple{N,Int}, p_l::NTuple{N,T}, p_r::NTuple{N,T}, breakpoints::NTuple{M,T}) where {N,M,T}
    @assert M == (N + 1)
    v = T[]
    for k in 1:N
        w = refined_uniform_unit_linear_range(p_l[k], p_r[k], n[k])
        push!(v, map_refined_range(w, breakpoints[k], breakpoints[k+1])...)
    end
    return IncreasingVector(unique(v))
end

# center of Interval and element
center(I::Interval) = 0.5 * (I.a+I.b)
center(I...) = map(center,I...)

# measure intervals
measure(I::Interval) = I.b-I.a

# element area
measure(element) = measure(element.data[1]) * measure(element.data[2])
