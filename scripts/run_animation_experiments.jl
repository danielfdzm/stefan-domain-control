#!/usr/bin/env julia
# Dedicated README animation experiments; the paper records in data/runs and
# data/summary.csv are independent and are not overwritten by this driver.
using StefanDomainControl, LinearAlgebra, Printf, TOML, Dates

BLAS.set_num_threads(1)
const ROOT = normpath(joinpath(@__DIR__, ".."))
const ANIMATION_DATA = joinpath(ROOT, "data", "animation")
const ANIMATION_CONFIG = (N=127, M=512, K=32, b=0.4, Rmin=0.45, Rmax=1.6)
mkpath(ANIMATION_DATA)

extinction(d) = d.terminal / maximum(d.norms)

const REFINED_REPLAY_TOLERANCE = 1e-12

function animation_experiment(name, p)
    println("Running ", name); flush(stdout)
    t0 = time()
    r = solve_control(p; verbose=true)
    d = diagnostics(r)

    # Replay the computed control without imposing a terminal condition.
    same = replay(r)
    same_state_error = maximum(abs, same.W - r.W)
    same_radius_error = maximum(abs, same.R - r.R)
    @assert same_state_error < 1e-9
    @assert same_radius_error < 1e-9
    @assert !p.constrained || abs(r.R[end] - p.Rf) < 1e-10
    @assert all(p.b .< r.R)
    passive = replay(r; uncontrolled=true)
    dp = diagnostics(passive)
    save_run(ANIMATION_DATA, name, r; passive=passive)
    # Confirm the archived controls and trajectories reproduce the solved run.
    stored = load_run(ANIMATION_DATA, name)
    @assert stored.W == r.W && stored.R == r.R && stored.C == r.C
    # The large expansion produces transformed temperatures of order 20. On the
    # refined grid, the default absolute 2e-14 threshold falls below roundoff.
    fine = replay(r; N=2p.N+1, M=2p.M, tol=REFINED_REPLAY_TOLERANCE)
    df = diagnostics(fine)
    @assert df.pde < 2e-7 && df.stefan < 2e-7

    row = (name=name, N=p.N, M=p.M, K=p.K, amplitude=p.amplitude,
           profile=string(p.profile), T=p.T, b=p.b, R0=p.R0, Rf=p.Rf,
           Rmin=p.Rmin, Rmax=p.Rmax, eta=p.eta, constrained=p.constrained,
           h1=d.h1, delta0=d.delta0, peak=maximum(d.norms), terminal=d.terminal,
           relative=d.relative, extinction=extinction(d),
           replay_terminal=df.terminal, replay_relative=df.relative,
           replay_extinction=extinction(df), uncontrolled_relative=dp.relative,
           replay_pde=df.pde, replay_stefan=df.stefan,
           energy=d.energy, energy_parallel=r.info.energy_parallel,
           energy_orthogonal=r.info.energy_orthogonal, psinorm=r.info.psinorm,
           lambda=r.info.lambda, excursion=d.excursion, Rlo=d.Rlo, Rhi=d.Rhi,
           Rfinal=d.Rfinal, radius_error=d.radius_error,
           replay_Rfinal=df.Rfinal, replay_radius_error=df.radius_error,
           uncontrolled_Rfinal=dp.Rfinal, radius_h1=d.radius_h1,
           pde=d.pde, stefan=d.stefan, same_state_error=same_state_error,
           same_radius_error=same_radius_error, iterations=length(r.history),
           defect=r.history[end].defect, optimality=r.info.optimality,
           constraint=r.info.constraint, seconds=round(time()-t0, digits=2))
    @printf("%s: %.1fs | extinction %.2e (replay %.2e) | R(T)=%.12f (replay %.12f) | %d iterations\n",
            name, row.seconds, row.extinction, row.replay_extinction,
            row.Rfinal, row.replay_Rfinal, row.iterations)
    flush(stdout)
    row
end

# All three cases share one initial temperature and one fixed actuator. The
# actuator must lie strictly inside the smallest prescribed terminal ball.
cases = (("animation_Rf_0.50", Config(; ANIMATION_CONFIG..., Rf=0.5)),
         ("animation_free_radius", Config(; ANIMATION_CONFIG..., Rf=1.0, constrained=false)),
         ("animation_Rf_1.50", Config(; ANIMATION_CONFIG..., Rf=1.5)))
rows = [animation_experiment(name, p) for (name, p) in cases]
open(joinpath(ANIMATION_DATA, "summary.csv"), "w") do io
    println(io, join(string.(keys(first(rows))), ","))
    for row in rows
        println(io, join(values(row), ","))
    end
end
open(joinpath(ANIMATION_DATA, "provenance.toml"), "w") do io
    TOML.print(io, Dict(
        "julia_version" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
        "solver" => "backward Euler ALE in s=r/R(t), volume Stefan law p'=-W_s(1,t), constrained Tikhonov HUM by SVD, relaxed Picard iteration",
        "reaction" => "f(y)=-4*tanh(y)", "fixed_point_tolerance" => 1e-9,
        "forward_iteration_tolerance" => 2e-14, "relaxation" => 0.8,
        "refined_replay_iteration_tolerance" => REFINED_REPLAY_TOLERANCE,
        "experiment_count" => length(rows), "date" => string(Dates.today()),
        "purpose" => "README animation with prescribed terminal radii 0.5 and 1.5 and a temperature-only reference",
        "replay_state_nodes" => 2ANIMATION_CONFIG.N+1,
        "replay_time_steps" => 2ANIMATION_CONFIG.M))
end
println("Completed ", length(rows), " animation experiments in ", ANIMATION_DATA)
