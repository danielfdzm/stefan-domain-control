#!/usr/bin/env julia
# Deterministic experiment driver. Every run is recomputed from scratch and overwrites
# its records in data/runs; the summary table is rewritten after every run.
# Optional arguments select groups: fan pure mesh eta horizon oscillatory landscape.
using StefanDomainControl, LinearAlgebra, Printf, TOML, Dates
BLAS.set_num_threads(1)
root = normpath(joinpath(@__DIR__, ".."))
const DATA = joinpath(root, "data"); const RUNS = joinpath(DATA, "runs")
mkpath(RUNS)
const ONLY = ARGS
group(name) = isempty(ONLY) || name in ONLY
rows = NamedTuple[]

extinction(d) = d.terminal / maximum(d.norms)

function write_summary()
    open(joinpath(DATA, "summary.csv"), "w") do io
        println(io, join(string.(keys(first(rows))), ","))
        for row in rows; println(io, join(values(row), ",")); end
    end
end

function experiment(name, p; refined=true, save=true, psi=false)
    print("Running ", name, " ... "); flush(stdout)
    t0 = time()
    r = solve_control(p; verbose=false)
    d = diagnostics(r)
    # Fresh nonlinear forward solve with the stored control: nothing is prescribed.
    same = replay(r)
    @assert maximum(abs, same.W - r.W) < 1e-9
    @assert abs(same.R[end] - r.R[end]) < 1e-9
    fine = refined ? replay(r; N=2p.N+1, M=2p.M) : same
    df = diagnostics(fine)
    passive = replay(r; uncontrolled=true)
    dp = diagnostics(passive)
    save && save_run(RUNS, name, r; psi=psi, passive=passive)
    row = (name=name, N=p.N, M=p.M, K=p.K, amplitude=p.amplitude, profile=string(p.profile),
           T=p.T, b=p.b, R0=p.R0, Rf=p.Rf, eta=p.eta, constrained=p.constrained,
           h1=d.h1, delta0=d.delta0, peak=maximum(d.norms), terminal=d.terminal,
           relative=d.relative, extinction=extinction(d),
           replay_terminal=df.terminal, replay_relative=df.relative, replay_extinction=extinction(df),
           uncontrolled_relative=dp.relative,
           energy=d.energy, energy_parallel=r.info.energy_parallel,
           energy_orthogonal=r.info.energy_orthogonal, psinorm=r.info.psinorm, lambda=r.info.lambda,
           excursion=d.excursion, Rlo=d.Rlo, Rhi=d.Rhi, Rfinal=d.Rfinal, radius_error=d.radius_error,
           replay_Rfinal=df.Rfinal, replay_radius_error=df.radius_error,
           uncontrolled_Rfinal=dp.Rfinal, radius_h1=d.radius_h1,
           pde=d.pde, stefan=d.stefan, iterations=length(r.history),
           defect=r.history[end].defect, optimality=r.info.optimality,
           constraint=r.info.constraint, seconds=round(time()-t0, digits=2))
    push!(rows, row)
    write_summary()
    @printf("%.1fs | extinction %.2e (replay %.2e) | R(T)-Rf %.1e (replay %.1e) | cost %.4f | %d iterations\n",
        row.seconds, row.extinction, row.replay_extinction, d.Rfinal-p.Rf,
        df.Rfinal-p.Rf, d.energy, row.iterations)
    flush(stdout)
    r
end

const FINE = (N=127, M=512, K=32)

# Prescribed terminal radii from one initial state; the free-radius control is the
# temperature-only reference.
if group("fan")
    for Rf in (0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15)
        experiment(@sprintf("fan_Rf_%.2f", Rf), Config(; FINE..., Rf=Rf); psi=(Rf == 0.95))
    end
    experiment("free_radius", Config(; FINE..., constrained=false))
end

# Pure domain control: zero initial temperature, prescribed terminal radius.
if group("pure")
    for Rf in (0.90, 1.10)
        experiment(@sprintf("pure_Rf_%.2f", Rf), Config(; FINE..., Rf=Rf, amplitude=0.0); psi=true)
    end
end

# Joint refinement of state grid, time grid and control space at fixed regularization.
if group("mesh")
    for (N, M, K) in ((31,128,8), (63,256,16), (127,512,32), (255,1024,64))
        experiment("mesh_$N", Config(N=N, M=M, K=K, Rf=0.95))
    end
end

# Regularization study: the volume constraint is exact for every eta.
if group("eta")
    for eta in (1e-4, 1e-6, 1e-8, 1e-12)
        experiment("eta_$(round(Int, log10(eta)))", Config(Rf=0.95, eta=eta))
    end
end

# Cost of shorter horizons and smaller fixed actuators.
if group("horizon")
    for T in (0.15, 0.20, 0.35); experiment("time_$T", Config(Rf=0.95, T=T)); end
    for b in (0.40, 0.475, 0.65); experiment("support_$b", Config(Rf=0.95, b=b)); end
end

# Sign-changing initial temperature with contraction, return and expansion targets.
if group("oscillatory")
    for Rf in (0.97, 1.00, 1.03)
        experiment(@sprintf("oscillatory_Rf_%.2f", Rf), Config(profile=:oscillatory, Rf=Rf))
    end
end

# Local regime: control cost over the plane of initial amplitude and radius change.
if group("landscape")
    for a in 0:0.005:0.04, dR in -0.05:0.01:0.05
        experiment(@sprintf("landscape_a%.3f_dR%+.2f", a, dR), Config(amplitude=a, Rf=1+dR);
                   refined=false, save=false)
    end
end

open(joinpath(DATA, "provenance.toml"), "w") do io
    TOML.print(io, Dict(
        "julia_version" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
        "solver" => "backward Euler ALE in s=r/R(t), volume Stefan law p'=-W_s(1,t), constrained Tikhonov HUM by SVD, relaxed Picard iteration",
        "reaction" => "f(y)=-4*tanh(y)", "fixed_point_tolerance" => 1e-9,
        "forward_iteration_tolerance" => 2e-14, "relaxation" => 0.8,
        "experiment_count" => length(rows), "date" => string(Dates.today())))
end
println("Completed ", length(rows), " experiments.")
