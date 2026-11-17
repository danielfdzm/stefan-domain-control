#!/usr/bin/env julia
# LaTeX and CSV tables from data/summary.csv.
using StefanDomainControl, DelimitedFiles, Printf, LinearAlgebra
root = normpath(joinpath(@__DIR__, ".."))
out = joinpath(root, "tables"); mkpath(out)
data, header = readdlm(joinpath(root, "data", "summary.csv"), ',', Any; header=true)
cols = Dict(string(k) => i for (i,k) in enumerate(vec(header)))
rows = Dict(string(data[i,1]) => i for i in axes(data,1))
value(name, key) = data[rows[name], cols[key]]

# Reconstruct the initial norm with exactly the spherical quadrature used by
# diagnostics. Landscape replays use the design mesh; all other paper runs use
# the doubled mesh, whose initial norm differs slightly by quadrature.
function initialnorm(name; replay=false)
    N = Int(value(name, "N"))
    replay && !startswith(name, "landscape_") && (N = 2N+1)
    p = Config(N=N, R0=Float64(value(name, "R0")),
               amplitude=Float64(value(name, "amplitude")),
               profile=Symbol(value(name, "profile")))
    sqrt(4pi*p.R0/(N+1))*norm(initial_state(p))
end

function sci(x)
    x == 0 && return "\$0\$"
    e = floor(Int, log10(abs(x)))
    @sprintf("\$%.2f\\times10^{%d}\$", x/10.0^e, e)
end
function table(path, caption, label, colspec, head, body, csvhead, csvbody)
    open(joinpath(out, path*".tex"), "w") do io
        println(io, "\\begin{table}[tbp]")
        println(io, "\\caption{", caption, "}")
        println(io, "\\label{", label, "}")
        println(io, "\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}")
        println(io, "\\begin{tabular}{", colspec, "}")
        println(io, "\\toprule")
        println(io, head, "\\\\")
        println(io, "\\midrule")
        for line in body; println(io, line, " \\\\"); end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
    end
    open(joinpath(out, path*".csv"), "w") do io
        println(io, csvhead)
        for line in csvbody; println(io, line); end
    end
    println("Wrote ", path)
end

# Prescribed terminal radii from one initial state (fine mesh), with the free-radius reference.
body = String[]; csv = String[]
for name in ["fan_Rf_0.85","fan_Rf_0.90","fan_Rf_0.95","fan_Rf_1.00","fan_Rf_1.05","fan_Rf_1.10","fan_Rf_1.15","free_radius"]
    Rf = value(name,"constrained") == true ? @sprintf("%.2f", value(name,"Rf")) : "free"
    push!(body, @sprintf("%s & %.4f & %.4f & %.4f & %s & %s & %s", Rf,
        value(name,"energy"), value(name,"energy_parallel"), value(name,"energy_orthogonal"),
        sci(value(name,"extinction")), sci(value(name,"replay_extinction")),
        value(name,"constrained") == true ? sci(value(name,"replay_radius_error")) : @sprintf("\\(R(T)=%.5f\\)", value(name,"replay_Rfinal"))))
    push!(csv, join([Rf, value(name,"energy"), value(name,"energy_parallel"), value(name,"energy_orthogonal"),
        value(name,"extinction"), value(name,"replay_extinction"), value(name,"replay_Rfinal"), value(name,"replay_radius_error")], ","))
end
table("targets",
    "Prescribed radii for \$a=0.02\$ with the baseline discretization. Replay uses the doubled mesh. The last row drops the volume constraint and reports the replay terminal radius instead of its error.",
    "tab:num_targets", "ccccccc",
    "\$R_f\$ & \$\\|u\\|_{L^2}\$ & \$\\|u_\\parallel\\|_{L^2}\$ & \$\\|u_\\perp\\|_{L^2}\$ & Extinction & Replay extinction & Replay \$|R(T)-R_f|\$",
    body, "Rf,energy,energy_parallel,energy_orthogonal,extinction,replay_extinction,replay_Rfinal,replay_radius_error", csv)

# Absolute errors and initial-relative accuracy, alongside peak-normalized
# extinction. Undefined zero-initial-data ratios are explicit, never set to zero.
body = String[]; csv = String[]
for name in ["fan_Rf_0.85","fan_Rf_0.90","fan_Rf_0.95","fan_Rf_1.00",
             "fan_Rf_1.05","fan_Rf_1.10","fan_Rf_1.15","free_radius",
             "pure_Rf_0.90","pure_Rf_1.10"]
    Rf = value(name,"constrained") == true ? @sprintf("%.2f", value(name,"Rf")) : "free"
    ini = initialnorm(name); ini_replay = initialnorm(name; replay=true)
    rel = ini_replay > 0 ? value(name,"replay_terminal")/ini_replay : NaN
    reltext = isfinite(rel) ? @sprintf("%.3g", 100rel) : "---"
    push!(body, @sprintf("%.2f & %s & %s & %s & %s & %s", value(name,"amplitude"), Rf,
        sci(ini), sci(value(name,"terminal")), sci(value(name,"replay_terminal")), reltext))
    push!(csv, join([name, value(name,"amplitude"), Rf, ini, ini_replay,
        value(name,"terminal"), value(name,"replay_terminal"), rel], ","))
end
table("terminal_accuracy",
    "Absolute spherical temperature norms and initial-relative error for the prescribed-radius and zero-temperature experiments. Replay uses the doubled mesh. The last column is a percentage, using the initial norm on the replay mesh; it is undefined for zero initial temperature.",
    "tab:num_terminal_accuracy", "cccccc",
    "\$a\$ & \$R_f\$ & \$\\|y_0\\|_{L^2}\$ & Design \$\\|y(T)\\|_{L^2}\$ & Replay \$\\|y(T)\\|_{L^2}\$ & Replay/initial (\\%)",
    body, "name,amplitude,Rf,initial,replay_initial,terminal,replay_terminal,replay_initial_relative", csv)

# Machine-readable accuracy for all 126 cases, including the same-grid
# landscape replays. Keep the summary archive unchanged.
open(joinpath(out, "terminal_accuracy_all.csv"), "w") do io
    println(io, "name,replay_grid,N,replay_N,initial,replay_initial,terminal,initial_relative,replay_terminal,replay_initial_relative,extinction,replay_extinction")
    for name in string.(data[:,1])
        N = Int(value(name,"N")); doubled = !startswith(name,"landscape_")
        ini = initialnorm(name); ini_replay = initialnorm(name; replay=true)
        rel = ini > 0 ? value(name,"terminal")/ini : NaN
        replayrel = ini_replay > 0 ? value(name,"replay_terminal")/ini_replay : NaN
        @assert isequal(isnan(rel), isnan(value(name,"relative")))
        @assert !isfinite(rel) || isapprox(rel, value(name,"relative"); rtol=1e-12)
        @assert !isfinite(replayrel) || isapprox(replayrel, value(name,"replay_relative"); rtol=1e-12)
        println(io, join([name, doubled ? "doubled" : "design", N, doubled ? 2N+1 : N,
            ini, ini_replay, value(name,"terminal"), rel, value(name,"replay_terminal"),
            replayrel, value(name,"extinction"), value(name,"replay_extinction")], ","))
    end
end

# Mesh refinement at Rf = 0.95.
body = String[]; csv = String[]
for name in ["mesh_31","mesh_63","mesh_127","mesh_255"]
    push!(body, @sprintf("%d & %d & %d & %s & %s & %s & %.4f & %.5f", value(name,"N"), value(name,"M"), value(name,"K"),
        sci(value(name,"extinction")), sci(value(name,"replay_extinction")), sci(value(name,"replay_radius_error")),
        value(name,"energy"), value(name,"psinorm")))
    push!(csv, join([value(name,"N"), value(name,"M"), value(name,"K"), value(name,"extinction"), value(name,"replay_extinction"),
        value(name,"replay_radius_error"), value(name,"energy"), value(name,"psinorm")], ","))
end
table("mesh_refinement",
    "Joint refinement for \$R_f=0.95\$, \$a=0.02\$ and \$\\eta=10^{-10}\$, with replay on the doubled mesh. The design radius error is at rounding level; the last column is the norm \$|g|\$ of the discrete volume gradient in the physical control metric.",
    "tab:num_mesh", "rrrccccc",
    "\$N\$ & \$M\$ & \$K\$ & Extinction & Replay extinction & Replay \$|R(T)-R_f|\$ & \$\\|u\\|_{L^2}\$ & \$|g|\$",
    body, "N,M,K,extinction,replay_extinction,replay_radius_error,energy,volume_gradient_norm", csv)

# Regularization at Rf = 0.95 on the intermediate mesh.
body = String[]; csv = String[]
for (name, e) in (("eta_-4",-4),("eta_-6",-6),("eta_-8",-8),("mesh_63",-10),("eta_-12",-12))
    push!(body, @sprintf("\$10^{%d}\$ & %s & %s & %s & %s & %.4f", e, sci(value(name,"extinction")),
        sci(value(name,"replay_extinction")), sci(value(name,"radius_error")), sci(value(name,"replay_radius_error")), value(name,"energy")))
    push!(csv, join([10.0^e, value(name,"extinction"), value(name,"replay_extinction"), value(name,"radius_error"),
        value(name,"replay_radius_error"), value(name,"energy")], ","))
end
table("regularization",
    "Regularization for \$R_f=0.95\$ on \$(63,256,16)\$, with replay on the doubled mesh. The terminal volume constraint holds for every \$\\eta\$.",
    "tab:num_eta", "cccccc",
    "\$\\eta\$ & Extinction & Replay extinction & \$|R(T)-R_f|\$ & Replay \$|R(T)-R_f|\$ & \$\\|u\\|_{L^2}\$",
    body, "eta,extinction,replay_extinction,radius_error,replay_radius_error,energy", csv)

# Horizon and actuator sweeps at Rf = 0.95.
body = String[]; csv = String[]
for (label, name) in (("\$T=0.15\$","time_0.15"),("\$T=0.20\$","time_0.2"),("\$T=0.25\$","mesh_63"),("\$T=0.35\$","time_0.35"),
                      ("\$b=0.40\$","support_0.4"),("\$b=0.475\$","support_0.475"),("\$b=0.55\$","mesh_63"),("\$b=0.65\$","support_0.65"))
    push!(body, @sprintf("%s & %.4f & %.4f & %.4f & %s & %s", label, value(name,"energy"), value(name,"energy_parallel"),
        value(name,"energy_orthogonal"), sci(value(name,"replay_extinction")), sci(value(name,"replay_radius_error"))))
    push!(csv, join([replace(label, "\$"=>""), value(name,"energy"), value(name,"energy_parallel"), value(name,"energy_orthogonal"),
        value(name,"replay_extinction"), value(name,"replay_radius_error")], ","))
end
table("horizon_actuator",
    "Horizon and actuator sensitivity for \$R_f=0.95\$, \$a=0.02\$ on \$(63,256,16)\$, with \$b=0.55\$ in the horizon study and \$T=0.25\$ in the actuator study.",
    "tab:num_horizon", "lccccc",
    "Parameter & \$\\|u\\|_{L^2}\$ & \$\\|u_\\parallel\\|_{L^2}\$ & \$\\|u_\\perp\\|_{L^2}\$ & Replay extinction & Replay \$|R(T)-R_f|\$",
    body, "parameter,energy,energy_parallel,energy_orthogonal,replay_extinction,replay_radius_error", csv)
