#!/usr/bin/env julia
# Verify the stored records, figures and tables against the checksum manifest, check the
# recorded tolerances of every run, and recompute two small experiments from scratch.
# `--write` regenerates the manifest instead of verifying it.
using StefanDomainControl, SHA, DelimitedFiles, LinearAlgebra, Printf, TOML
BLAS.set_num_threads(1)
root = normpath(joinpath(@__DIR__, ".."))
digest(path) = open(io -> bytes2hex(sha256(io)), path)
tracked = String[]
for dir in ("data", "figures", "tables"), (base, _, files) in walkdir(joinpath(root, dir)), f in files
    startswith(f, ".") || push!(tracked, relpath(joinpath(base, f), root))
end
sort!(tracked)
manifest = joinpath(root, "checksums.sha256")
if "--write" in ARGS
    open(manifest, "w") do io
        for rel in tracked; println(io, digest(joinpath(root, rel)), "  ", rel); end
    end
    println("Wrote ", manifest, " with ", length(tracked), " entries.")
else
    listed = Dict{String,String}()
    for line in eachline(manifest)
        isempty(strip(line)) && continue
        h, rel = split(line; limit=2)
        listed[strip(rel)] = h
    end
    for rel in tracked
        haskey(listed, rel) || error("File not in the manifest: $rel")
        digest(joinpath(root, rel)) == listed[rel] || error("Checksum mismatch: $rel")
    end
    for rel in keys(listed)
        rel in tracked || error("File listed in the manifest is missing: $rel")
    end
    println("All ", length(tracked), " records, figures and tables match the checksum manifest.")
end

data, header = readdlm(joinpath(root, "data", "summary.csv"), ',', Any; header=true)
cols = Dict(string(k) => i for (i,k) in enumerate(vec(header)))
column(key) = data[:, cols[key]]
names = string.(column("name"))
constrained = column("constrained") .== true
all(Float64.(column("pde")) .< 2e-7) || error("PDE residual above tolerance")
all(Float64.(column("stefan")) .< 2e-7) || error("Stefan residual above tolerance")
all(Float64.(column("defect")) .< 1e-9) || error("Fixed-point defect above tolerance")
# The relative KKT residual is normalized by the larger of the initial gradient and the
# constraint force. With zero initial temperature both are of the order of eta, so double
# precision leaves a relative residual near 1e-6 there and near 1e-12 elsewhere.
all(Float64.(column("optimality")) .< 1e-5) || error("KKT residual above tolerance")
amplitude = Float64.(column("amplitude"))
all(Float64.(column("optimality"))[amplitude .> 0] .< 1e-10) || error("KKT residual above tolerance for nonzero initial data")
all(Float64.(column("radius_error"))[constrained] .< 1e-12) || error("Prescribed terminal radius not reached on the design grid")
all(Float64.(column("constraint"))[constrained] .< 1e-12) || error("Volume constraint residual above tolerance")
all(Float64.(column("replay_radius_error"))[constrained] .< 2e-3) || error("Replay radius error above tolerance")
energy = Float64.(column("energy"))
parallel = Float64.(column("energy_parallel")); orthogonal = Float64.(column("energy_orthogonal"))
all(abs.(energy.^2 .- parallel.^2 .- orthogonal.^2) .< 1e-8 .* max.(energy.^2, 1e-12)) || error("Control decomposition is not orthogonal")
prov = TOML.parsefile(joinpath(root, "data", "provenance.toml"))
prov["experiment_count"] == length(names) == 126 || error("Unexpected number of experiments")
println("All ", length(names), " runs satisfy the recorded tolerances; every prescribed terminal radius is met to rounding on the design grid.")

# The README animation has its own actuator, radius bounds and records. Its
# extreme targets do not inherit the paper suite's finer-grid error thresholds.
animation_dir = joinpath(root, "data", "animation")
animation_prov = TOML.parsefile(joinpath(animation_dir, "provenance.toml"))
animation_prov["experiment_count"] == 3 || error("Unexpected number of animation experiments")
for (name, target) in (("animation_Rf_0.50", 0.5),
                       ("animation_free_radius", nothing),
                       ("animation_Rf_1.50", 1.5))
    stored = load_run(animation_dir, name)
    d = diagnostics(stored)
    stored.p.b == 0.4 || error("Unexpected animation actuator")
    stored.p.constrained == (target !== nothing) || error("Incorrect animation constraint mode")
    all(stored.p.b .< stored.R) || error("Animation actuator leaves the moving ball")
    all(stored.p.Rmin .< stored.R .< stored.p.Rmax) || error("Animation radius out of bounds")
    d.pde < 2e-7 && d.stefan < 2e-7 || error("Animation equations fail verification")
    if target !== nothing
        stored.p.Rf == target || error("Incorrect animation target")
        abs(d.Rfinal-target) < 1e-10 || error("Animation target not reached")
    end
    same = replay(stored)
    maximum(abs, same.W-stored.W) < 1e-9 || error("Animation state differs from independent replay")
    maximum(abs, same.R-stored.R) < 1e-9 || error("Animation radius differs from independent replay")
    @printf("Verified %s: terminal radius %.12f, extinction ratio %.3e\n",
        name, d.Rfinal, d.terminal/maximum(d.norms))
end

for (name, cfg) in (("mesh_31", Config(N=31, M=128, K=8, Rf=0.95)),
                    ("oscillatory_Rf_1.00", Config(profile=:oscillatory, Rf=1.00)))
    i = findfirst(==(name), names)
    r = solve_control(cfg; verbose=false); d = diagnostics(r)
    isapprox(d.energy, data[i, cols["energy"]]; rtol=1e-8) || error("Recomputed control norm differs for $name")
    abs(d.terminal - data[i, cols["terminal"]]) < 1e-9 || error("Recomputed terminal norm differs for $name")
    stored = load_run(joinpath(root, "data", "runs"), name)
    maximum(abs, stored.W - r.W) < 1e-9 || error("Recomputed state differs from the stored record for $name")
    maximum(abs, stored.R - r.R) < 1e-9 || error("Recomputed radius differs from the stored record for $name")
    @printf("Recomputed %s: control norm %.6f, terminal radius %.12f, extinction ratio %.3e\n",
        name, d.energy, d.Rfinal, d.terminal/maximum(d.norms))
end
println("Check passed.")
