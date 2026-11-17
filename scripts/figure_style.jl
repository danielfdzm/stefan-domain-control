# Shared visual style of all figures: palette, theme and export helpers.
using CairoMakie, Printf
# Curves: LaTeX blue and red first, then black and green; further colors follow xcolor names.
const BLUE = "#0000FF"; const RED = "#FF0000"; const BLACK = "#000000"; const GREEN = "#008000"
const ORANGE = "#FF8000"; const VIOLET = "#800080"; const BROWN = "#BF8040"
const GRAY = "#808080"; const INK = "#1A2742"
const CURVES = [BLUE, RED, BLACK, GREEN, ORANGE, VIOLET, BROWN]
# Signed fields use cividis with a symmetric range; nonnegative fields use Blues.
const FIELD = :cividis
const SEQUENTIAL = :Blues
const FAN_TARGETS = (0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15)
fancolor(Rf) = CURVES[findfirst(==(Rf), FAN_TARGETS)]
"Blue for contraction, red for expansion, when only two targets are compared."
paircolor(Rf) = Rf < 1 ? BLUE : RED
decimal2(x) = @sprintf("%.2f", x)
decimal3(x) = @sprintf("%.3f", x)

# At the manuscript's 5.9-inch figure width, 28 units give approximately
# 10-point labels for the 1180-unit, two-column figures.
set_theme!(Theme(fontsize=28, fonts=(;regular="TeX Gyre Heros", bold="TeX Gyre Heros Bold"),
    textcolor=INK, backgroundcolor=:white,
    Axis=(;titlealign=:center, titlefont=:bold, spinewidth=1.1,
        titlesize=28, xlabelsize=28, ylabelsize=28, xticklabelsize=26, yticklabelsize=26,
        xticks=WilkinsonTicks(4), yticks=WilkinsonTicks(4),
        xgridcolor=(:gray,0.13), ygridcolor=(:gray,0.13),
        topspinevisible=false, rightspinevisible=false),
    Axis3=(;titlealign=:center, titlefont=:bold),
    Colorbar=(;labelsize=26, ticklabelsize=24),
    Legend=(;labelsize=24),
    Lines=(;linewidth=3.2), Scatter=(;markersize=12)))

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RUNS = joinpath(ROOT, "data", "runs")
const PAPER = joinpath(ROOT, "figures", "paper")
const PREVIEWS = joinpath(ROOT, "figures", "previews")

function exportfigure(f, name)
    mkpath(PAPER); mkpath(PREVIEWS)
    save(joinpath(PAPER, name*".pdf"), f; pt_per_unit=0.75)
    save(joinpath(PREVIEWS, name*".png"), f; px_per_unit=2)
    println("Saved ", name); flush(stdout)
end

"Range of the finite values of one or several arrays, for a shared color scale."
fieldrange(arrays...) = extrema(Iterators.flatten(filter(isfinite, vec(A)) for A in arrays))

"Symmetric range of signed fields, preserving zero at the colormap midpoint."
function signedrange(arrays...)
    lo, hi = fieldrange(arrays...)
    bound = max(abs(lo), abs(hi), eps(Float64))
    (-bound, bound)
end

"Linear interpolation of nodal values on the uniform grid of [0,1] at ρ ∈ [0,1]."
function radial_interp(values, ρ)
    ξ = clamp(ρ, 0, 1)*(length(values)-1)
    j = min(floor(Int, ξ)+1, length(values)-1); θ = ξ-(j-1)
    (1-θ)*values[j] + θ*values[j+1]
end

"Temperature on a fixed physical grid of r; NaN outside the moving ball."
function physical_field(result; nr=320, rmax=maximum(result.R))
    Y = physical_state(result); p = result.p
    rr = collect(range(0, rmax, length=nr))
    F = fill(NaN, nr, p.M+1)
    for n in 1:p.M+1, (j,r) in enumerate(rr)
        r <= result.R[n] || continue
        F[j,n] = radial_interp(view(Y,:,n), r/result.R[n])
    end
    rr, F
end

"Hemispherical cutaway of the ball of radius R with the radial temperature on a diametral disk."
function cutaway!(ax, values, R, b; colorrange, nθ=65, nρ=100, nφ=160)
    θ = range(0, pi, length=nθ); φ = range(0, pi, length=nθ)
    X = [R*sin(v)*cos(u) for u in θ, v in φ]
    YY = [R*sin(v)*sin(u) for u in θ, v in φ]
    Z = [R*cos(v) for u in θ, v in φ]
    surface!(ax, X, YY, Z, color=zeros(size(X)), colormap=FIELD, colorrange=colorrange,
        shading=NoShading, rasterize=3)
    for latitude in (-pi/3, -pi/6, 0, pi/6, pi/3)
        lines!(ax, R*cos(latitude).*cos.(θ), R*cos(latitude).*sin.(θ),
            fill(R*sin(latitude), length(θ)), color=(:white,0.35), linewidth=0.7)
    end
    ρ = range(0, 1, length=nρ); ψ = range(0, 2pi, length=nφ)
    Xd = [R*v*cos(u) for v in ρ, u in ψ]
    Zd = [R*v*sin(u) for v in ρ, u in ψ]
    colors = [radial_interp(values, v) for v in ρ, u in ψ]
    surface!(ax, Xd, fill(-0.001, size(Xd)), Zd, color=colors, colormap=FIELD,
        colorrange=colorrange, shading=NoShading, rasterize=3)
    lines!(ax, b.*cos.(ψ), fill(-0.006, length(ψ)), b.*sin.(ψ), color=BLACK, linestyle=:dash, linewidth=2.2)
    lines!(ax, R.*cos.(ψ), fill(-0.004, length(ψ)), R.*sin.(ψ), color=BLACK, linewidth=1.4)
    nothing
end

function ball_axis(where; title="", limit=1.2)
    Axis3(where, title=title, aspect=:data, azimuth=-0.25pi, elevation=0.13pi, perspectiveness=0.15,
        limits=(-limit,limit,-limit,limit,-limit,limit),
        xlabel="x₁", ylabel="x₂", zlabel="x₃", xlabeloffset=42, ylabeloffset=42, zlabeloffset=42,
        xlabelsize=34, ylabelsize=34, zlabelsize=34,
        xticklabelsize=30, yticklabelsize=30, zticklabelsize=30,
        xticks=[-1,0,1], yticks=[-1,0,1], zticks=[-1,0,1])
end

"Spherical L2 norms in time for a run."
timenorms(result) = diagnostics(result).norms
"L2(B_b) norm of the physical control on each time interval."
function controlnorms(result)
    p = result.p; Mq = mass_control(p)
    [sqrt(dot(result.C[:,n], Mq*result.C[:,n])) for n in 1:size(result.C,2)]
end
