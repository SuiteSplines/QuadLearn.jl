using ImmersedSplines, Printf, Algoim

function algoim_compute_area(ϕ, partition, order)
    # area
    A = 0.0

    # npts
    n = 0

    # compute area
    for element in Elements(partition)
        Q = ImmersedQuadRule(ϕ, element; order = order, phase=-1)
        #@show typeof(Q.x)
        A += sum(Q.w)
        n += length(Q.w)
    end

    return A
end

function quadlearn_compute_area(ϕ, partition, order)
    # model loader
    loader = init_cached_model_loader(; maxsize=75)

    # order of quadrature
    order = order

    # order of interpolation model
    interp_order = order * order

    # phase to integrate (1: outside, -1: inside)
    phase = -1

    # get patch weights
    weights = get_patch_weights(ϕ, partition; loader=loader, order=order, interp_order=interp_order, n=32, phase=phase)

    # area and npts
    sum(weights)
end


# geometry
offset = [-3.0, -1.5, -0.5, -0.1, 0.1, 0.5]
ϕ = Algoim.AlgoimCallLevelSetFunction(
    (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*(x[1])/2) - sin(3*(x[2])/2) + offset[2]  ), 
    (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
)
#ϕ = implicit_circle_def(2.0)

# embedding box width
width = 6.0

# generate partition
ref_partition = IncreasingRange(-width/2, width/2, 1024) ⨱ IncreasingRange(-width/2, width/2, 1024)
#ref_partition = IncreasingRange(0.0, 2.0, 512) ⨱ IncreasingRange(0.0, 2.0, 512)

# reference area
A_ref = algoim_compute_area(ϕ, ref_partition, 3)

# order
order = 3

for nelem = [4,8,16,32,64,128] .+ 1
    @info "Nelem = $(nelem-1)"

    # generate partition
    partition = IncreasingRange(-width/2, width/2, nelem) ⨱ IncreasingRange(-width/2, width/2, nelem)
    #partition = IncreasingRange(0.0, 2.0, nelem) ⨱ IncreasingRange(0.0, 2.0, nelem)

    # algoim
    A_algoim = algoim_compute_area(ϕ, partition, order)
    x = abs(A_ref - A_algoim)
    sci_str = @sprintf("%.4e", x)
    latex_str = replace(sci_str, r"e([+-]?\d+)" => s" \\cdot 10^{\1}")
    println("abs(A_ref - A_algoim) = \$$latex_str\$")

    # quadlearn
    A_quadlearn = quadlearn_compute_area(ϕ, partition, order)
    x = abs(A_ref - A_quadlearn)
    sci_str = @sprintf("%.4e", x)
    latex_str = replace(sci_str, r"e([+-]?\d+)" => s" \\cdot 10^{\1}")
    println("abs(A_ref - A_quadlearn) = \$$latex_str\$")
end

