using QuadLearnData, SortedSequences, IgaFormation
using GLMakie, Makie.GeometryBasics
using Algoim
using StatsBase
using ImplicitGeometries
using StaticArrays

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
width = 1.0

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
    (x) -> -( (x[1]/2)^4 + (x[2]/2)^4 + cos(3*(x[1])/2) - sin(3*(x[2])/2) + offset[2]  ), 
    (x) -> -[ (x[1]^3/8 - 3 / 2 * sin(3*x[1]/2)), (x[2]^3 / 8 - 3/2 * cos(3*x[2]/2)) ]
)

#ϕ = Algoim.AlgoimCallLevelSetFunction(
#    (x) -> (sqrt(0.25x[1]^2 + 0.25x[2]^2) - (1.25 - 0.5 * cos(2*angle(x[1] + im * x[2]) + pi/3)^4 ) ), 
#    (x) -> [0.0,0.0]
#)

#ϕ = implicit_circle_def(h + 0h, -0.5, -0.7)
#ϕ = AlgoimCallLevelSetFunction(
#    x -> -(-x[1]*x[1] - x[2]*x[2] + (1/3)^2),  # ϕ
#    x -> -[-2.0*x[1], -2.0*x[2]]         # ∇ϕ
#)

w = 1/3
h = 1/6
ϕ = AlgoimCallLevelSetFunction(
    x -> -(-x[1]*x[1] / w^2 - x[2]*x[2] / h^2 + 1.0),  # ϕ
    x -> [-2.0*x[1]/w^2, -2.0*x[2]/h^2]         # ∇ϕ
)

#airfoil = [
#    SVector(-1.0,  0.0), SVector( 0.91070,  0.39097), SVector( 1.0,  0.0),
#    SVector( 1.0,  0.0), SVector( 1.06054, -0.26505), SVector( 0.4, -0.1),
#    SVector( 0.4, -0.1), SVector(-0.20900,  0.05217), SVector(-1.0,  0.0)
#];
#geometry = Rotation(Ring(QuadraticBezier(; v=airfoil); r=0.08); θ=deg2rad(15))
#gradient = ImplicitGeometries.Gradient(geometry)
#ϕ = AlgoimCallLevelSetFunction(
#    x -> geometry(SVector(x)),  # ϕ
#    x -> gradient(SVector(x))   # ∇ϕ
#)
##constrain = true
#contour_constrain = false
## !nelem = 128

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
vlines!(ax2, partition.data[1], color = (:gray, 0.25))
hlines!(ax2, partition.data[2], color = (:gray, 0.25))

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
        # least square circular approximation of the cut
        R, a, b = levelset_least_squares_circle(ϕ, Ωₑ; n=32, constrain=contour_constrain)

        # define exact levelset for lsq approximation
        ψ = implicit_circle_def(R, a, b)

        # get contour coordinates for approximate cut and plot
        #color = (mod(k, 2) == 0) ? RGBf(0.749, 0.141, 0.157) : RGBf(0.973, 0.71, 0.251)
        color = RGBf(0.973, 0.71, 0.251)
        cs = levelset_contours(ψ, Ωₑ; n=256)
        for c in cs
            #lines!(ax1, c[1], c[2], color=color, linewidth=6)
            #lines!(ax2, c[1], c[2], color=color, linewidth=6)
            lines!(ax1, c[1], c[2], color=color, linewidth=2)
            lines!(ax2, c[1], c[2], color=color, linewidth=2)
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

# plot quadrature points from models
if plot_quad_points
    # model loader
    loader = init_cached_model_loader(; maxsize=75)

    # order of quadrature
    order = 3

    # order of interpolation model
    interp_order = order^2

    # phase to integrate (1: outside, -1: inside)
    phase = 1

    # get patch stencil
    stencil = get_patch_stencil(partition; order=order+1)

    # get patch weights
    weights = get_patch_weights(ϕ, partition; loader=loader, order=order, interp_order=interp_order, n=32, phase=phase, constrain=constrain)

    # collect quadrature points
    x = map(k -> first(stencil[k]), 1:length(stencil))
    y = map(k -> last(stencil[k]), 1:length(stencil))
    quadrature_points = Point2f.(x, y)

    # plot quadrature points 
    markersize = abs.(weights[:]).^(1/4)
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