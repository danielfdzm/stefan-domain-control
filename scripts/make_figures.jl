#!/usr/bin/env julia
# Publication figures from the stored records in data/runs and data/summary.csv.
# Optional arguments select figures: fan pure balls landscape structure surfaces.
using StefanDomainControl, CairoMakie, DelimitedFiles, LinearAlgebra, Printf
include(joinpath(@__DIR__, "figure_style.jl"))
const ONLY = ARGS
want(name) = isempty(ONLY) || name in ONLY
readrun(name) = load_run(RUNS, name)

function summary_table()
    data, header = readdlm(joinpath(ROOT, "data", "summary.csv"), ',', Any; header=true)
    cols = Dict(string(k) => i for (i,k) in enumerate(vec(header)))
    rows = Dict(string(data[i,1]) => i for i in axes(data,1))
    (data=data, cols=cols, rows=rows)
end
value(S, name, key) = S.data[S.rows[name], S.cols[key]]

# ---------------------------------------------------------------------------------
# 1. One initial state, seven prescribed terminal radii.
if want("fan"); let
    runs = Dict(Rf => readrun(@sprintf("fan_Rf_%.2f", Rf)) for Rf in FAN_TARGETS)
    free = readrun("free_radius")
    base = runs[0.85]; p = base.p
    t = collect(range(0, p.T, length=p.M+1))
    passive = base.passive
    f = Figure(size=(1180, 980))
    # Plot the signed physical temperature at the centre, without norm or normalization.
    centre_temperature(run) = vec(physical_state(run)[1,:])
    temperatures = [centre_temperature(run) for run in (passive, free, (runs[Rf] for Rf in FAN_TARGETS)...)]
    temperature_min, temperature_max = fieldrange(temperatures...)
    temperature_padding = 0.10max(temperature_max-temperature_min, eps(Float64))
    ax = Axis(f[1,1], title=L"\text{a) Temperature extinction}", xlabel="time t", xticks=[0.0,0.1,0.2],
        ylabel=L"y(0,t)", yticks=WilkinsonTicks(5; k_min=4))
    hlines!(ax, [0.0], color=(GRAY,0.5), linewidth=1.2)
    lines!(ax, t, temperatures[1], color=GRAY, linestyle=:dash, label=L"\text{uncontrolled}")
    lines!(ax, t, temperatures[2], color=GRAY, linestyle=:dot, linewidth=3.2, label=L"free $R(T)$")
    for Rf in FAN_TARGETS
        lines!(ax, t, centre_temperature(runs[Rf]), color=fancolor(Rf), label=L"R_f = %$(decimal2(Rf))")
    end
    ylims!(ax, temperature_min-temperature_padding, temperature_max+temperature_padding)
    Legend(f[3,1:2], ax, orientation=:horizontal, framevisible=false, labelsize=24,
        nbanks=2, patchsize=(24,10), colgap=18, rowgap=5, tellheight=true)
    ax = Axis(f[1,2], title=L"\text{b) Prescribed terminal radii}", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="radius R(t)")
    hlines!(ax, [p.Rmin, p.Rmax], color=(GRAY,0.8), linestyle=:dash, linewidth=1.2)
    lines!(ax, t, passive.R, color=GRAY, linestyle=:dash, label="uncontrolled")
    lines!(ax, t, free.R, color=GRAY, linestyle=:dot, linewidth=3.2, label=L"free $R(T)$")
    for Rf in FAN_TARGETS
        lines!(ax, t, runs[Rf].R, color=fancolor(Rf))
        scatter!(ax, [p.T], [Rf], color=:white, strokecolor=fancolor(Rf), strokewidth=2.6, markersize=13)
    end
    text!(ax, p.T-0.012, 1.175, text=L"targets $R_f$", align=(:right,:center), fontsize=24, color=GRAY)
    ylims!(ax, 0.79, 1.21)
    gleft = GridLayout(f[2,1]); gright = GridLayout(f[2,2])
    rr, F = physical_field(base; rmax=1.0)
    ax = Axis(gleft[1,1], title=L"c) Temperature, $R_f = %$(decimal2(base.p.Rf))$",
        xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="physical radius r")
    m = heatmap!(ax, t, rr, permutedims(F), colormap=FIELD, colorrange=signedrange(F),
        nan_color=:transparent, rasterize=3)
    lines!(ax, t, base.R, color=BLACK, linewidth=2)
    hlines!(ax, [p.b], color=BLACK, linestyle=:dash, linewidth=2)
    text!(ax, 0.004, p.b, text="actuator radius b", align=(:left,:bottom), fontsize=24, color=BLACK)
    text!(ax, 0.004, 0.89, text="free boundary R(t)", align=(:left,:top), fontsize=24, color=BLACK)
    Colorbar(gleft[1,2], m, label="temperature y", width=16)
    rc, V = physical_control(base)
    ax = Axis(gright[1,1], title=L"\text{d) Control in the fixed ball}", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="physical radius r")
    m = heatmap!(ax, (collect(1:p.M).-0.5).*p.T/p.M, rc, permutedims(V), colormap=FIELD,
        colorrange=signedrange(V), rasterize=3)
    Colorbar(gright[1,2], m, label="control u", width=16)
    rowgap!(f.layout, 34); colgap!(f.layout, 35)
    exportfigure(f, "fig_fan")
end; end

# ---------------------------------------------------------------------------------
# 2. Pure domain control from zero temperature.
if want("pure"); let
    grow = readrun("pure_Rf_1.10"); shrink = readrun("pure_Rf_0.90")
    p = grow.p; t = collect(range(0, p.T, length=p.M+1))
    f = Figure(size=(1180, 900))
    fields = Dict(:grow => physical_field(grow; rmax=1.1), :shrink => physical_field(shrink; rmax=1.1))
    crange = signedrange(fields[:grow][2], fields[:shrink][2])
    for (j, (key, run, title)) in enumerate(((:grow, grow, L"a) Expansion, $R_f = 1.10$"),
                                             (:shrink, shrink, L"b) Contraction, $R_f = 0.90$")))
        rr, F = fields[key]
        g = GridLayout(f[1,j])
        ax = Axis(g[1,1], title=title, xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="physical radius r")
        m = heatmap!(ax, t, rr, permutedims(F), colormap=FIELD, colorrange=crange,
            nan_color=:transparent, rasterize=3)
        lines!(ax, t, run.R, color=BLACK, linewidth=2)
        hlines!(ax, [p.b], color=BLACK, linestyle=:dash, linewidth=2)
        hlines!(ax, [run.p.Rf], color=paircolor(run.p.Rf), linestyle=:dot, linewidth=2)
        ylims!(ax, 0, 1.12)
        j == 2 && Colorbar(g[1,2], m, label="temperature y", width=16)
    end
    ax = Axis(f[2,1], title=L"\text{c) Radius trajectories}", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="radius R(t)")
    for run in (grow, shrink)
        lines!(ax, t, run.R, color=paircolor(run.p.Rf), label=L"R_f = %$(decimal2(run.p.Rf))")
        scatter!(ax, [p.T], [run.p.Rf], color=:white, strokecolor=paircolor(run.p.Rf), strokewidth=2.6, markersize=13)
    end
    hlines!(ax, [p.R0], color=(GRAY,0.8), linestyle=:dash, linewidth=1.2)
    ax = Axis(f[2,2], title=L"\text{d) Temperature pulse}", xlabel="time t", xticks=[0.0,0.1,0.2],
        ylabel="3D L² temperature norm")
    for run in (grow, shrink)
        lines!(ax, t, timenorms(run), color=paircolor(run.p.Rf), label=L"R_f = %$(decimal2(run.p.Rf))")
    end
    # One shared legend outside the data axes prevents labels from covering the pulse.
    Legend(f[3,1:2], ax, orientation=:horizontal, framevisible=false, labelsize=26, tellheight=true)
    rowgap!(f.layout, 34); colgap!(f.layout, 35)
    exportfigure(f, "fig_pure_domain")
end; end

# ---------------------------------------------------------------------------------
# 3. Hemispherical cutaways: the ball visibly shrinks or grows while the temperature vanishes.
if want("balls"); let
    cases = ((readrun("fan_Rf_0.85"), L"R_f = 0.85"), (readrun("fan_Rf_1.15"), L"R_f = 1.15"))
    fractions = (0.0, 0.25, 0.5, 1.0)
    snapshots(run) = [round(Int, 1 + q*run.p.M) for q in fractions]
    crange = signedrange((physical_state(run)[:, snapshots(run)] for (run,_) in cases)...)
    f = Figure(size=(1500, 820))
    for (j, q) in enumerate(fractions)
        Label(f[0,j], L"t / T = %$(decimal2(q))", fontsize=34, font=:bold, tellwidth=false)
    end
    for (i, (run, label)) in enumerate(cases)
        Label(f[i,0], label, fontsize=34, font=:bold, rotation=pi/2, tellheight=false)
        p = run.p; Y = physical_state(run); t = range(0, p.T, length=p.M+1)
        for (j, n) in enumerate(snapshots(run))
            ax = ball_axis(f[i,j])
            cutaway!(ax, Y[:,n], run.R[n], p.b; colorrange=crange)
        end
    end
    Colorbar(f[3,1:4], colormap=FIELD, limits=crange, vertical=false,
        label="temperature y", labelsize=32, ticklabelsize=30,
        height=17, width=Relative(0.6), tellwidth=false)
    Label(f[4,1:4], "Dashed circle: actuator radius b     •     Solid circle: free boundary R(t)",
        fontsize=30, tellwidth=false)
    # The narrow row-label column and non-sizing headers leave the four axes
    # free to use the entire available figure width.
    colsize!(f.layout, 0, Fixed(46))
    rowgap!(f.layout, 15); colgap!(f.layout, 22)
    exportfigure(f, "fig_balls")
end; end

# ---------------------------------------------------------------------------------
# 4. Local regime: cost landscape over amplitude and prescribed radius change.
if want("landscape"); let
    S = summary_table()
    names = [n for n in keys(S.rows) if startswith(n, "landscape")]
    amps = sort(unique([value(S, n, "amplitude") for n in names]))
    dRs = sort(unique([round(value(S, n, "Rf") - value(S, n, "R0"); digits=6) for n in names]))
    cost = fill(NaN, length(dRs), length(amps)); share = copy(cost)
    for n in names
        i = findfirst(==(round(value(S,n,"Rf")-value(S,n,"R0"); digits=6)), dRs)
        j = findfirst(==(value(S,n,"amplitude")), amps)
        cost[i,j] = value(S,n,"energy")
        share[i,j] = value(S,n,"energy") > 0 ? value(S,n,"energy_parallel")/value(S,n,"energy") : NaN
    end
    f = Figure(size=(1680, 520))
    ax = Axis(f[1,1], title=L"Control cost $\Vert u\Vert$", xlabel=L"radius change $R_f - R_0$", ylabel="initial amplitude a")
    m = contourf!(ax, dRs, amps, cost, levels=range(0, maximum(cost), length=16), colormap=SEQUENTIAL)
    contour!(ax, dRs, amps, cost, levels=0.5:0.5:3.5, color=(BLACK,0.35), linewidth=0.8)
    Colorbar(f[1,2], m, width=13)
    ax = Axis(f[1,3], title=L"Cost against $\delta_0$", xlabel=L"\delta_0 = \Vert\nabla y_0\Vert + |R_0 - R_f|", ylabel=L"control cost $\Vert u\Vert$")
    δs = [value(S,n,"delta0") for n in names]; costs = [value(S,n,"energy") for n in names]
    as = [value(S,n,"amplitude") for n in names]
    slope = maximum(filter(isfinite, costs ./ δs))
    lines!(ax, [0, maximum(δs)], [0, slope*maximum(δs)], color=BLACK, linestyle=:dash, linewidth=1.6,
        label=L"\Vert u\Vert = %$(decimal2(slope))\,\delta_0")
    sc = scatter!(ax, δs, costs, color=as, colormap=SEQUENTIAL, strokecolor=BLACK, strokewidth=0.6, markersize=10)
    axislegend(ax, position=:lt, framevisible=false, labelsize=24)
    Colorbar(f[1,4], sc, width=13, label="initial amplitude a", labelsize=24)
    ax = Axis(f[1,5], title=L"Share of $u_{\parallel}$", xlabel=L"radius change $R_f - R_0$", ylabel="initial amplitude a")
    m = heatmap!(ax, dRs, amps, share, colormap=SEQUENTIAL, colorrange=(0,1), nan_color=:transparent)
    Colorbar(f[1,6], m, width=13, label=L"\Vert u_{\parallel}\Vert / \Vert u\Vert", labelsize=24)
    colgap!(f.layout, 18)
    exportfigure(f, "fig_landscape")
end; end

# ---------------------------------------------------------------------------------
# 5. Structure of the control: auxiliary adjoint and the two orthogonal components.
if want("structure"); let
    base = readrun("fan_Rf_0.95"); p = base.p
    t = collect(range(0, p.T, length=p.M+1)); tm = (collect(1:p.M).-0.5).*p.T/p.M
    s = grid(p)
    # The decomposition is evaluated at the converged frozen coefficients.
    _, _, _, psi, info = frozen_control(base.W, base.R, p)
    dt = p.T/p.M; Uc = cholesky(mass_control(p)).U
    xpar = (dot(info.gvec, info.x)/info.psinorm^2) .* info.gvec
    tocontrol(x) = (Uc \ reshape(x, p.K, p.M)) ./ sqrt(dt)
    Cpar = tocontrol(xpar); Cort = tocontrol(info.x .- xpar)
    rc, Vpar = physical_control((p=p, C=Cpar)); _, Vort = physical_control((p=p, C=Cort))
    f = Figure(size=(1180, 900), figure_padding=(20,20,20,34))
    # Normalize using the same right-endpoint radii as the frozen state operators.
    # The last spatial node carries the one-sided flux-stencil artifact; explicitly
    # omit it from the continuum comparison, keeping all remaining values unclipped.
    psi_unit = normalized_auxiliary_adjoint(psi, base.R, p)
    psi_interior = psi_unit[1:end-1,:]
    ax = Axis(f[1,1], title=L"a) Auxiliary adjoint $\psi$", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="reference radius s = r / R(t)")
    m = heatmap!(ax, t[2:end], s[1:end-1], permutedims(psi_interior), colormap=SEQUENTIAL,
        colorrange=fieldrange(psi_interior), rasterize=3)
    lines!(ax, t, p.b./base.R, color=BLACK, linestyle=:dash, linewidth=2)
    text!(ax, 0.004, p.b/base.R[1]-0.02, text="s < b / R(t)", align=(:left,:top), fontsize=24, color=BLACK)
    Colorbar(f[1,2], m, width=16)
    crange = signedrange(Vpar, Vort)
    ax = Axis(f[1,3], title=L"b) Volume-adjusting part $u_{\parallel}$", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="physical radius r")
    m = heatmap!(ax, tm, rc, permutedims(Vpar), colormap=FIELD, colorrange=crange, rasterize=3)
    Colorbar(f[1,4], m, width=16)
    ax = Axis(f[2,1], title=L"c) Orthogonal part $u_{\perp}$", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="physical radius r")
    m = heatmap!(ax, tm, rc, permutedims(Vort), colormap=FIELD, colorrange=crange, rasterize=3)
    Colorbar(f[2,2], m, width=16)
    ax = Axis(f[2,3], title=L"\text{d) Control norms}", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel=L"$L^2(B_b)$ norm")
    lines!(ax, tm, controlnorms(base), color=BLACK, label=L"\text{total control}")
    lines!(ax, tm, controlnorms((p=p, C=Cpar)), color=BLUE, label=L"volume-adjusting $u_{\parallel}$")
    lines!(ax, tm, controlnorms((p=p, C=Cort)), color=RED, label=L"volume-neutral $u_{\perp}$")
    Legend(f[3,1:4], ax, orientation=:horizontal, framevisible=false, labelsize=24, tellheight=true)
    Label(f[0,1:4], L"R_f = %$(decimal2(p.Rf)):\quad \Vert u\Vert = %$(decimal3(info.energy)),\quad \Vert u_{\parallel}\Vert = %$(decimal3(info.energy_parallel)),\quad \Vert u_{\perp}\Vert = %$(decimal3(info.energy_orthogonal))",
        fontsize=26, color=(INK,0.8), tellwidth=false)
    rowgap!(f.layout, 26); colgap!(f.layout, 18)
    exportfigure(f, "fig_structure")
end; end

# ---------------------------------------------------------------------------------
# 6. Space-time surfaces on the moving ball and inside the fixed actuator.
if want("surfaces"); let
    run = readrun("fan_Rf_0.85"); p = run.p
    t = collect(range(0, p.T, length=p.M+1)); s = collect(range(0, 1, length=p.N+2))
    Y = physical_state(run); rc, V = physical_control(run)
    f = Figure(size=(1320, 575), figure_padding=(75, 20, 20, 20))
    inds = unique(round.(Int, range(1, p.M+1, length=100)))
    TT = repeat(reshape(t[inds], 1, :), length(s), 1)
    RR = s*transpose(run.R[inds])
    ax = Axis3(f[1,1], title=L"Temperature on the ball shrinking to $R_f = %$(decimal2(p.Rf))$", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="radius r",
        zlabel="y", azimuth=0.30pi, elevation=0.20pi, perspectiveness=0.25,
        xlabeloffset=48, ylabeloffset=48, zlabeloffset=75, aspect=(1.3,1,0.8))
    surface!(ax, TT, RR, Y[:,inds], colormap=FIELD, colorrange=signedrange(Y), shading=NoShading, rasterize=3)
    lines!(ax, t, run.R, zeros(length(t)), color=BLACK, linewidth=2.5)
    ax = Axis3(f[1,2], title=L"Control inside the fixed ball $B_b$", xlabel="time t", xticks=[0.0,0.1,0.2], ylabel="radius r",
        zlabel="u", azimuth=0.30pi, elevation=0.20pi, perspectiveness=0.25,
        xlabeloffset=48, ylabeloffset=48, zlabeloffset=75, aspect=(1.3,1,0.8))
    ii = unique(round.(Int, range(1, p.M, length=100)))
    surface!(ax, (collect(1:p.M)[ii].-0.5).*p.T/p.M, rc, permutedims(V[:,ii]), colormap=FIELD,
        colorrange=signedrange(V), shading=NoShading, rasterize=3)
    exportfigure(f, "fig_surfaces")
end; end
