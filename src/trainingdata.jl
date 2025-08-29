export compute_moments_bernstein_basis, integrate_bernstein_basis, get_cut_domain, area_cut_element
export generate_training_data, getdata, get_bernsteinfuns_at_gausspoints

function integrate_bernstein_basis(q, map, element; n = 16)#32)

    # domain of cut element
    Ix, Iy = element.data

    # compute bernstein polynomials
    bx = x -> bernstein(Ix, q, x)
    by = y -> bernstein(Iy, q, y)

    # allocate array for moments
    A = zeros(q+1,q+1)

    # get quadrature rule
    partition = IncreasingRange(Ix.a, Ix.b, n+1) ⨱ IncreasingRange(Iy.a, Iy.b, n+1)
    
    # loop over integration mesh
    for e in Elements(partition)
        Q = ImmersedQuadRule(map, e; order = 12)
        integrate_bernstein_basis!(A, Q, bx, by)
    end
    
    return A
end

function integrate_bernstein_basis!(A::Matrix, Q::ImmersedQuadRule, Bx::Function, By::Function)

    @assert size(A,1) == size(A,2) "Matrix is not square."
    q = size(A,2)-1

    # loop over quadrature points
    for l in 1:length(Q)

        # get quadrature point
        x, y, w = Q.x[l][1], Q.x[l][2], Q.w[l]

        # compute bernstein polynomials
        bx = Bx(x)
        by = By(y)

        # loop over functions
        for j in 1:q+1
            for i in 1:q+1
                A[i,j] += bx[i] * by[j] * w
            end
        end
    end
end

function integrate_bernstein_basis(q, element)

    # domain of cut element
    Ix, Iy = element.data

    # quadrature rules
    Qx, Qy = GaussRule(Legendre, q+1, Ix), GaussRule(Legendre, q+1, Iy)
    
    # compute moments of Bernstein functions
    return  bernstein(Ix, q, Qx.x) * Qx.w * (bernstein(Iy, q, Qy.x) * Qy.w)'
end

@memoize function get_bernsteinfuns_at_gausspoints(; order)
    x = get_stencil(Interval(0.0,1.0), order = order+1)
    b = bernstein(Interval(0.0,1.0), order, x)
    return KroneckerProduct(b, b)
end

struct TrainingSet{S<:CartesianProduct, A<:KroneckerProduct, T<:Real}
    order::Int
    input::S
    operator::A
    x::Array{T}
    y::Array{T}
end

function generate_training_data(; order, data, case, w=1.0)

    # arrays to store training data
    xtrain = zeros((order+1)^2, length(data))
    ytrain = zeros((order+1)^2, length(data))

    # operator mapping weights to moments
    A = get_bernsteinfuns_at_gausspoints(order=order)

    # data-point
    @showprogress dt=1 desc="Computing..." for k in 1:length(data)

        # current data-point
        R, θ̂, ẑ = data[k]

        # get parametrization for element of width = w
        p = Parametrization(R, 1.0)

        # get cut-element
        element = get_cut_domain(p, θ̂, ẑ, case=case)

        # compute distance function at Gauss-Legendre stencil
        x, y = collect(get_stencil(element, order = order+1))
        distance = evaluate_distance_function(R, x, y)

        # compute moments
        moments = integrate_bernstein_basis(order, implicit_circle_def(R), element; n=1)

        # save training data
        xtrain[:,k] = distance[:]
        ytrain[:,k] = moments[:]
    end

    return TrainingSet(order, data, A, xtrain, ytrain)
end

function getdata(; order, R, θ̂, ẑ, case)
    traindata = generate_training_data(order=order, data=CartesianProduct(R, θ̂, ẑ), case=case)
    testdata = generate_training_data(order=order, data=CartesianProduct(midpoints(R), midpoints(θ̂), midpoints(ẑ)), case=case)
    return traindata, testdata
end