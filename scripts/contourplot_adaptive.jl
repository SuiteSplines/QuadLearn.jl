using QuadLearnData, SortedSequences, IgaFormation
using GLMakie, Makie.GeometryBasics
using Algoim
using StatsBase
using ImplicitGeometries
using StaticArrays
using SpecialSpaces
using RegionTrees

#import CairoMakie

GLMakie.activate!()

# get figure
f = Figure()

# set figure size
screen = display(f)
resize!(screen, 1200, 600)

# zoom out axis
ax1 = Axis(f[1,1], aspect = DataAspect(), xgridvisible = false, ygridvisible = false)

# zoom in axis
ax2 = Axis(f[1,2], aspect = DataAspect(), xgridvisible = false, ygridvisible = false)

# hide clutter
hidedecorations!.([ax1, ax2])
hidespines!(ax1)
hidespines!(ax2)

# flag for quadrature points plots
plot_quad_points = true


# empty axes
empty!.([ax1, ax2])

# embedding box width and height
width = 4.0

# set initial limits 
limits!(ax1, -width/2, width/2, -width/2, width/2)
limits!(ax2, -width/2, width/2, -width/2, width/2)

# number of elements in each direction
nelem = 16

# generate partition
partition = IncreasingRange(-width/2, width/2, nelem+1) ⨱ IncreasingRange(-width/2, width/2, nelem+1)
h = width/(nelem-1)
constrain = true
contour_constrain = constrain

# get embedding domain
domain = get_cut_domain(width, 0.0, 0.0)

# levelset function function (interesting offsets: -1.5 / -0.5 / -0.1 / +0.1)
offset = [-3.0, -1.5, -0.5, -0.1, 0.1, 0.5]
ϕ = Algoim.AlgoimCallLevelSetFunction(
    (x) -> ( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*(x[1])/2) - sin(3*(x[2])/2) + offset[2]  ), 
    (x) -> [ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
)

#ϕ = Algoim.AlgoimCallLevelSetFunction(
#    (x) -> (sqrt(0.25x[1]^2 + 0.25x[2]^2) - (1.25 - 0.5 * cos(2*angle(x[1] + im * x[2]) + pi/3)^4 ) ), 
#    (x) -> [0.0,0.0]
#)

#ϕ = implicit_circle_def(h + 0h, -0.5, -0.7)
#ϕ = AlgoimCallLevelSetFunction(
#    x -> -x[1]*x[1] - x[2]*x[2] + (1/3)^2,  # ϕ
#    x -> [-2.0*x[1], -2.0*x[2]]         # ∇ϕ
#)

#constrain = true
contour_constrain = false
# !nelem = 128

#polygon = Ring(Polygon(;
#    v = [
#        SVector(-1.0, -1.0),
#        SVector( 1.0, -1.0),
#        SVector( 1.0,  1.0),
#        SVector( 0.0,  2.0),
#        SVector(-1.0,  1.0),
#    ]
#), r=0.1)
#
#
#geometry = SmoothSubtraction(polygon, Translation(Circle(; r=0.5), dx=-0.5, dy=0.75); k=1.0)
#geometry = SmoothSubtraction(geometry, Translation(Circle(; r=0.5), dx=0.95, dy=0.75); k=0.5)
#geometry = SmoothSubtraction(geometry, Translation(Circle(; r=0.5), dx=0.25, dy=-0.7); k=0.8)
#geometry = SmoothSubtraction(geometry, Translation(Circle(; r=0.38), dx=-0.45, dy=-0.7); k=0.05)
#geometry = Scal(geometry; s=0.8)
#gradient = ImplicitGeometries.Gradient(geometry)
#ϕ = AlgoimCallLevelSetFunction(
#    x -> geometry(SVector(x)),  # ϕ
#    x -> gradient(SVector(x))   # ∇ϕ
#)

a = 0.5
b = 0.25
tooth = Translation(Scaling(QuadraticBezier(; v=[
    SVector(-b, 0.0), SVector(-a, 1.0), SVector(0.0, 1.0),
    SVector(0.0, 1.0), SVector(a, 1.0), SVector(b, 0.0),
]), s=0.5), dy = 0.95)

geometry = ImplicitGeometries.Circle(; r=1.1)
Δα = π/6
for α = 0:Δα:2π - Δα
    geometry = SmoothUnion(geometry, Rotation(tooth; θ=α); k=0.05)
end
geometry = geometry \ ImplicitGeometries.Circle(; r=0.2)
geometry = Scaling(geometry; s=2/3)
ϕ = AlgoimCallLevelSetFunction(
    x -> geometry(SVector(x[1], x[2])),  # ϕ
    x -> gradient(SVector(x[1], x[2]))   # ∇ϕ
)


#r = 1.0
#dd = -2.0
#th = 0.45
#geometry = Translation(ImplicitGeometries.Circle(; r=r); dx = -dd, dy = dd)
#geometry = geometry ∪ Translation(ImplicitGeometries.Circle(r=r); dx = -dd, dy = -dd)
#geometry = geometry ∪ Translation(ImplicitGeometries.Circle(r=r); dx = dd, dy = dd)
#geometry = geometry ∪ Translation(ImplicitGeometries.Circle(r=r); dx = dd, dy = -dd)
#geometry = SmoothUnion(geometry, Translation(Rectangle(; w = 4width, h = th); dy = dd); k=0.3)
#geometry = SmoothUnion(geometry, Translation(Rectangle(; w = 4width, h = th); dy = -dd); k=0.3)
#geometry = SmoothUnion(geometry, Translation(Rectangle(; w = th, h = 4width); dx = dd); k=0.3)
#geometry = SmoothUnion(geometry, Translation(Rectangle(; w = th, h = 4width); dx = -dd); k=0.3)
#gradient = ImplicitGeometries.Gradient(geometry)
#ϕ = AlgoimCallLevelSetFunction(
#    x -> geometry(SVector(x)),  # ϕ
#    x -> gradient(SVector(x))   # ∇ϕ
#)

#ww, hh = 2/6, 1/6
#ϕ = AlgoimCallLevelSetFunction(
#    x -> x[1]*x[1] / ww^2 + x[2]*x[2] / hh^2 - 1.0,  # ϕ
#    x -> [2.0*x[1]/ww^2, 2.0*x[2]/hh^2]         # ∇ϕ
#)


# compute zero contour on domain
cs = levelset_contours(ϕ, domain; n=1024)
for c in cs
    lines!(ax1, c[1], c[2], color=:black, linewidth=5)
    lines!(ax2, c[1], c[2], color=:black, linewidth=5)
end

# plot contour on both axes

# grid lines (mesh)
vlines!(ax1, partition.data[1], color = (:gray, 0.25))
hlines!(ax1, partition.data[2], color = (:gray, 0.25))
#vlines!(ax2, partition.data[1], color = (:gray, 0.25))
#hlines!(ax2, partition.data[2], color = (:gray, 0.25))

# plot element
element_corners = Observable(Point2f[])
polygon = @lift(Makie.Polygon($element_corners))
poly!(ax1, polygon, color = (:gray, 0.5))
poly!(ax2, polygon, color = (:gray, 0.5))

# element loop
for (k,element) in enumerate(Elements(partition))
    # get element domain
    Ωₑ = get_element_domain(element)

    # update element plot
    element_corners[] = Point2f[Ωₑ[1], Ωₑ[2], Ωₑ[4], Ωₑ[3]] 

    # update zoom in axis limits
    C = QuadLearnData.center(Ωₑ)
    limits!(ax2, C[1] - h, C[1] + h, C[2] - h, C[2] + h)

    try
        # get contour coordinates for approximate cut and plot
        #color = (mod(k, 2) == 0) ? RGBf(0.749, 0.141, 0.157) : RGBf(0.973, 0.71, 0.251)
        color = RGBf(0.973, 0.71, 0.251)
        for c in cs
            #lines!(ax1, c[1], c[2], color=color, linewidth=6)
            #lines!(ax2, c[1], c[2], color=color, linewidth=6)
            #lines!(ax1, c[1], c[2], color=color, linewidth=4)
            #lines!(ax2, c[1], c[2], color=color, linewidth=4)
        end

        # wait
        #sleep(0.15)
    catch e
        if isa(e, StrictlyPositiveLevelSet)
            #sleep(0.15)
            nothing
        elseif isa(e, StrictlyNegativeLevelSet)
            #sleep(0.15)
            nothing
        else
            rethrow(e)
        end
    end

end

# get elements
elements = Elements(partition)

# dummy space
space = ScalarSplineSpace(1, partition)

# model loader
loader = init_cached_model_loader(; maxsize=75)

# init quadrule
order = 2
quadrule = AdaptiveCutcellQuadratureRule(;
    partition = partition, 
    mapping = ϕ,
    npoints = (order+1)^2,
    loader = loader,
    order = order,
    interp_order = order^2,
    n = 32,
    Rmax = 30,
    factor = 1.0
)

# accessor
acc = ElementAccessor(testspace=space, trialspace=space, quadrule=quadrule);

# compute quadrature
quadpoints, quadweights = [], []
for e in Elements(partition)
    Ωₑ = get_element_domain(e)
    refinery = QuadTreeRefinery(ϕ, 1.0, 32, 1.0)
    origin = SVector(Ωₑ.data[1][1], Ωₑ.data[2][1])
    widths = SVector(QuadLearnData.measure.(Ωₑ.data))
    try 
        RR, aa, bb = levelset_least_squares_circle(refinery.ϕ, Ωₑ; n=refinery.n, linear_max_constrain=true)
        root = Cell(origin, widths, (R=RR,a=aa,b=bb))
        adaptivesampling!(root, refinery)
        for leaf in allleaves(root)
            ψ = implicit_circle_def(leaf.data.R, leaf.data.a, leaf.data.b)
            v = hcat(collect(RegionTrees.vertices(leaf.boundary))...)
            lines!(ax2, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]], color = (:grey, 0.25), linewidth=1.0)
            ee = QuadLearnData.hyperrecangle_to_element(leaf.boundary)
            color = RGBf(0.973, 0.71, 0.251)
            cs = levelset_contours(ψ, ee; n=256)
            for c in cs
                lines!(ax1, c[1], c[2], color=color, linewidth=2)
                lines!(ax2, c[1], c[2], color=color, linewidth=2)
            end
        end
    catch ex
        lines!(ax2, [e[1,1][1], e[1,2][1], e[2,2][1], e[1,2][1]], [e[1,1][2], e[1,2][2], e[2,2][2], e[1,2][2]], color = (:grey, 0.25), linewidth=1.0)
        if !isa(ex, StrictlyPositiveLevelSet) && !isa(ex, StrictlyNegativeLevelSet)
            rethrow(ex)
        end 
    end

    Q = QuadratureRule(acc, e; phase=-1)
    push!(quadpoints, Q.x...)
    push!(quadweights, Q.w...)
end

# plot quadrature points from models
if plot_quad_points == true
    # collect quadrature points
    x = first.(quadpoints)
    y = last.(quadpoints)
    quadrature_points = Point2f.(x, y)

    # plot quadrature points 
    markersize = abs.(quadweights[:]).^(1/4)
    markersize = markersize / maximum(markersize) * 10
    colormap = :roma
    #scatter!(ax1, quadrature_points, colormap = colormap, color = -weights[:], markersize = markersize)
    #scatter!(ax2, quadrature_points, colormap = colormap, color = -weights[:], markersize = markersize * 2)
    scatter!(ax1, quadrature_points, color=RGBf(0.3,0.3,0.3), markersize = markersize * 0.8)
    scatter!(ax2, quadrature_points, color=RGBf(0.3,0.3,0.3), markersize = markersize * 1.2)
end

function update_zoom_in(pos)
    eind = closest_element(partition, pos)
    Ωₑ = get_element_domain(elements[eind...])
    C = QuadLearnData.center(Ωₑ)
    limits!(ax2, C[1] - h, C[1] + h, C[2] - h, C[2] + h)
    element_corners[] = Point2f[Ωₑ[1], Ωₑ[2], Ωₑ[4], Ωₑ[3]] 
end

function closest_element(partition, pos)
    x, y = pos[1], pos[2]
    xmid = midpoints(partition.data[1])
    ymid = midpoints(partition.data[2])
    xind = argmin(abs.(xmid .- x))
    yind = argmin(abs.(xmid .- y))
    return xind, yind
end


# set zoom in at center
update_zoom_in((0.0,h))

# disable rectangle zoom on ax1 and ax2
deregister_interaction!(ax1, :rectanglezoom)
deregister_interaction!(ax2, :rectanglezoom)

# add callback for element selection
on(events(ax1).mousebutton) do event
    if event.button == Mouse.left
        if event.action == Mouse.press || event.action == Mouse.release
            mp = mouseposition_px(ax1.scene)
            pos = to_world(ax1.scene, Point2f(mp))
            update_zoom_in(pos)
        end
    end
end

# add callback for arrow keys control
on(events(f).keyboardbutton) do event
    if event.action == Keyboard.press || event.action == Keyboard.repeat
        if event.key == Keyboard.up
            update_zoom_in(sum(element_corners[])/4 .+ [0, h])
        elseif event.key == Keyboard.down
            update_zoom_in(sum(element_corners[])/4 .+ [0, -h])
        elseif event.key == Keyboard.left
            update_zoom_in(sum(element_corners[])/4 .+ [-h, 0])
        elseif event.key == Keyboard.right
            update_zoom_in(sum(element_corners[])/4 .+ [h, 0])
        end
    end
end

display(screen)
resize!(screen, 1200, 600)