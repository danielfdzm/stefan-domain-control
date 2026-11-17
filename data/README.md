# Records of the experiments

All records are plain text and describe deterministic experiments. The solver,
recorded parameters and generation scripts are public. See
[the reproduction guide](../REPRODUCING.md) for commands that regenerate the
figures and tables from these records or recompute the experiments.

| File | Content |
| --- | --- |
| `summary.csv` | One row per experiment with the parameters, the measured norms and errors, the two control components, the discrete volume-gradient norm in the physical control metric, the multiplier, the residuals of the discrete equations and the runtime |
| `provenance.toml` | Julia version, BLAS threads, solver description, tolerances, experiment count and date |
| `runs/NAME_config.toml` | The `Config` of the run |
| `runs/NAME_state.csv` | No header; columns `t, R, W_1, ..., W_N` with `W(s,t) = R(t) s z(R(t) s, t)` at the interior nodes `s_i = i/(N+1)` |
| `runs/NAME_control.csv` | No header; columns are the right endpoint `t_n` of each time interval followed by the hat-function values of `q = r u` at `r_j = b j/(K+1)`, held constant on `(t_{n-1}, t_n]` |
| `runs/NAME_history.csv` | Outer-iteration defects, optimality residual, constraint residual, multiplier and control norm |
| `runs/NAME_passive_state.csv` | The uncontrolled trajectory from the same initial data on the same grid |
| `runs/NAME_psi.csv` | No header; rows `t_n, Lambda_1^n, ..., Lambda_N^n` contain the raw discrete adjoint of the integrated boundary flux for `fan_Rf_0.95` and the pure domain-control runs |

The 99 `landscape_*` runs appear only in `summary.csv`. Their `replay_*`
diagnostics use a same-grid independent replay. The other 27 paper runs use
a doubled-grid replay for these columns.

The adjoint CSV files store the raw discrete adjoint $\Lambda$. The normalized
continuous auxiliary adjoint is approximated by
$-\Lambda_i^n/[hR(t_n)^2]$, with $h=1/(N+1)$ and $R(t_n)$ from the matching
state file. The last spatial node is affected by the one-sided boundary-flux
stencil and is omitted only in the adjoint visualization; the archived arrays
and the control calculation retain it.

## README animation records

`animation/` contains three additional runs: `animation_Rf_0.50`,
`animation_free_radius` and `animation_Rf_1.50`. All three use actuator radius
`b=0.4` and admissible radius interval `(0.45, 1.6)`, with the same initial data
and mesh `(N,M,K)=(127,512,32)`. The constrained runs prescribe terminal radii
`0.5` and `1.5`; the middle comparison leaves the terminal radius free.

The `NAME_config.toml`, `NAME_state.csv`, `NAME_control.csv`,
`NAME_history.csv` and `NAME_passive_state.csv` files use the formats above.
`animation/summary.csv` records the design and replay diagnostics, and
`animation/provenance.toml` records the environment and experiment settings.
In the animation summary, `Rmin` and `Rmax` are the admissible bounds, and
`same_state_error` and `same_radius_error` are the maximum differences between
the optimized trajectory and its independent replay on the design mesh.
These records are separate from the 126 paper experiments in `summary.csv`.
Load an animation run with
`load_run("data/animation", "animation_Rf_0.50")`.

## Columns of `summary.csv`

`name, N, M, K, amplitude, profile, T, b, R0, Rf, eta, constrained` are the
parameters. `h1` is the gradient seminorm of the initial temperature on the
ball of radius `R0`, and `delta0 = h1 + |R0 - Rf|` is the smallness parameter
of the main theorem. `peak` and `terminal` are the maximal and terminal
spherical `L2` norms of the temperature on the design grid, `relative` divides
the terminal norm by the initial norm (undefined for zero initial temperature)
and `extinction` divides it by `peak`. The `replay_*` columns measure the same
quantities after a fresh nonlinear forward solve with the stored control on the
grid `(2N+1, 2M)` for the 27 non-landscape runs and on the design grid for the
99 landscape runs. Initial-relative errors use the initial norm on the
corresponding grid. `uncontrolled_*` refers to the zero control. `energy` is
the control norm in `L2(B_b x (0,T))`, split into `energy_parallel` (the
component along the discrete volume gradient $g$ that adjusts the volume) and
`energy_orthogonal` (volume-neutral part), with
`energy^2 = energy_parallel^2 + energy_orthogonal^2`. `psinorm` is the legacy
column name for $|g|$, the Euclidean norm of the integrated-flux gradient in the weighted physical-control coordinates;
it is not the unweighted spatial norm of the continuous auxiliary adjoint.
`lambda` is the multiplier of the volume constraint, `excursion = max |R(t) - R0|`,
`Rlo` and `Rhi` the extreme radii, `Rfinal` and `radius_error = |R(T) - Rf|` the terminal radius on the
design grid, `replay_Rfinal` and `replay_radius_error` the same after the
corresponding replay, `radius_h1` the discrete `H1(0,T)` norm of `R - Rf`, `pde` and
`stefan` the maximal absolute residuals of the discrete equations evaluated
with the computed radius, `iterations` the number of outer iterations,
`defect` the final fixed-point defect, `optimality` the relative KKT residual,
`constraint` the residual of the linear volume constraint and `seconds` the
runtime of the run including its replays.

## Derived terminal-accuracy tables

`scripts/make_tables.jl` writes `tables/terminal_accuracy.tex` and
`tables/terminal_accuracy.csv` for the seven targets, the free-radius reference
and the two zero-initial-temperature cases. The tables report absolute initial
and terminal spherical norms and replay errors relative to the initial norm.
Undefined initial-relative errors for zero initial temperature are marked
explicitly. `tables/terminal_accuracy_all.csv` contains these diagnostics for
all 126 runs and labels each replay as same-grid or doubled-grid. These
artifacts are derived from `summary.csv` and the configured initial profiles;
they do not replace the archived raw diagnostics.

## Loading a record

```julia
using StefanDomainControl
run = load_run("data/runs", "fan_Rf_0.85")
run.p            # Config
run.R            # radius at the M+1 time nodes
run.W            # N x (M+1) transformed state
run.C            # K x M control coefficients
Y = physical_state(run)          # temperature at the s-nodes including the endpoints
d = diagnostics(run)             # norms, cost, residuals, terminal radius error
```
