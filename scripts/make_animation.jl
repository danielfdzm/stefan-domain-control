#!/usr/bin/env julia
# Animated hemispherical cutaways: temperature-only control with free terminal radius
# against prescribed contraction and expansion, with the radius traces drawn in time.
using StefanDomainControl, CairoMakie, LinearAlgebra, Printf
include(joinpath(@__DIR__, "figure_style.jl"))
const ANIM = joinpath(ROOT, "figures", "animation"); mkpath(ANIM)
frames = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 90
frames >= 2 || error("The animation needs at least two frames")
const ANIMATION_RUNS = joinpath(ROOT, "data", "animation")
# Keep the README animation consistent with the public blue-gray-yellow palette.
const ANIMATION_FIELD = :cividis

cases = ((load_run(ANIMATION_RUNS, "animation_Rf_0.50"), L"prescribed $R_f = 0.5$"),
         (load_run(ANIMATION_RUNS, "animation_free_radius"), L"temperature only, free $R(T)$"),
         (load_run(ANIMATION_RUNS, "animation_Rf_1.50"), L"prescribed $R_f = 1.5$"))
p = cases[1][1].p; t = collect(range(0, p.T, length=p.M+1))
fields = [physical_state(run) for (run,_) in cases]
# Independent symmetric ranges put zero at the gray midpoint of cividis.
# Each amplitude covers the full trajectory and stays fixed throughout playback.
amplitudes = [maximum(abs, field) for field in fields]
cranges = [(-max(a, eps(Float64)), max(a, eps(Float64))) for a in amplitudes]
radius_limit = maximum(run.p.Rmax for (run, _) in cases)
colors = [BLUE, BLACK, RED]

"Cutaway whose radius and disk colors are updated through observables."
function animated_cutaway!(ax, b; colorrange, nθ=61, nρ=70, nφ=120)
    θ = range(0, pi, length=nθ); φ = range(0, pi, length=nθ)
    us = [sin(v)*cos(u) for u in θ, v in φ]; vs = [sin(v)*sin(u) for u in θ, v in φ]; ws = [cos(v) for u in θ, v in φ]
    X = Observable(copy(us)); Y = Observable(copy(vs)); Z = Observable(copy(ws))
    surface!(ax, X, Y, Z, color=zeros(size(us)), colormap=ANIMATION_FIELD, colorrange=colorrange, shading=NoShading)
    ρ = range(0, 1, length=nρ); ψ = range(0, 2pi, length=nφ)
    dx = [v*cos(u) for v in ρ, u in ψ]; dz = [v*sin(u) for v in ρ, u in ψ]
    DX = Observable(copy(dx)); DZ = Observable(copy(dz)); DC = Observable(zeros(size(dx)))
    surface!(ax, DX, fill(-0.001, size(dx)), DZ, color=DC, colormap=ANIMATION_FIELD, colorrange=colorrange, shading=NoShading)
    lines!(ax, b.*cos.(ψ), fill(-0.006, length(ψ)), b.*sin.(ψ), color=BLACK, linestyle=:dash, linewidth=2.2)
    BX = Observable(cos.(ψ)); BZ = Observable(sin.(ψ))
    lines!(ax, BX, fill(-0.004, length(ψ)), BZ, color=BLACK, linewidth=1.5)
    EX = Observable(cos.(θ)); EY = Observable(sin.(θ))
    lines!(ax, EX, EY, zeros(length(θ)), color=(:white,0.4), linewidth=0.8)
    function update!(values, R)
        X[] = R.*us; Y[] = R.*vs; Z[] = R.*ws
        DX[] = R.*dx; DZ[] = R.*dz
        DC[] = [radial_interp(values, v) for v in ρ, u in ψ]
        BX[] = R.*cos.(ψ); BZ[] = R.*sin.(ψ)
        EX[] = R.*cos.(θ); EY[] = R.*sin.(θ)
    end
    update!
end

f = Figure(size=(1600, 1000), figure_padding=(20, 65, 20, 20))
clock = Observable(L"t / T = %$(decimal2(0.0))")
Label(f[0,1:3], clock, fontsize=32, font=:bold, color=INK, tellwidth=false)
updaters = []
for (j, (run, label)) in enumerate(cases)
    ax = ball_axis(f[1,j]; title=label, limit=radius_limit)
    push!(updaters, animated_cutaway!(ax, run.p.b; colorrange=cranges[j]))
    Colorbar(f[2,j], colormap=ANIMATION_FIELD, limits=cranges[j], vertical=false,
        label=L"temperature $y$", height=14, width=Relative(0.8),
        ticks=([cranges[j][1], 0.0, cranges[j][2]],
               [@sprintf("%.3g", cranges[j][1]), "0", @sprintf("%.3g", cranges[j][2])]),
        tellwidth=false, labelsize=28)
end
ax = Axis(f[3,1:3], xlabel=L"time $t$", ylabel=L"radius $R(t)$", height=250)
hlines!(ax, [p.Rmin, p.Rmax], color=(GRAY,0.8), linestyle=:dash, linewidth=1.2)
traces = []
for (j, (run, label)) in enumerate(cases)
    tx = Observable(t[1:1]); ry = Observable(run.R[1:1])
    lines!(ax, tx, ry, color=colors[j], linewidth=3, label=label)
    head = Observable(Point2f[(t[1], run.R[1])])
    scatter!(ax, head, color=colors[j], markersize=12)
    run.p.constrained && scatter!(ax, [p.T], [run.p.Rf], color=:white, strokecolor=colors[j], strokewidth=2.6, markersize=14)
    push!(traces, (tx, ry, head))
end
xlims!(ax, 0, p.T*1.02); ylims!(ax, p.Rmin-0.02, p.Rmax+0.02)
axislegend(ax, position=:lt, framevisible=false, labelsize=25, nbanks=3)
rowgap!(f.layout, 24)

inds = unique(round.(Int, range(1, p.M+1, length=frames)))
inds = vcat(inds, fill(inds[end], 12))   # hold the final state
# Makie rounds fractional frame rates to integers. Repeat each frame at the
# original 15 fps to halve playback speed, including the final pause: 13.6 s.
inds = repeat(inds; inner=2)
record(f, joinpath(ANIM, "stefan_domain_control.gif"), inds; framerate=15, px_per_unit=1) do n
    clock[] = L"t / T = %$(decimal2(t[n]/p.T))"
    for (j, (run, _)) in enumerate(cases)
        updaters[j](fields[j][:,n], run.R[n])
        tx, ry, head = traces[j]
        tx[] = t[1:n]; ry[] = run.R[1:n]; head[] = Point2f[(t[n], run.R[n])]
    end
end
save(joinpath(ANIM, "stefan_domain_control_still.png"), f; px_per_unit=1)
println("Saved animation with ", length(inds), " frames and a still of the final state")
