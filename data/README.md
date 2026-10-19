# Records of the experiments

All records are plain text and describe deterministic experiments. The solver
and recorded parameters are public; the experiment-generation scripts are
maintained privately.

| File | Content |
| --- | --- |
| `summary.csv` | One row per experiment with the parameters, the measured norms and errors, the two control components, the discrete norm of the auxiliary adjoint, the multiplier, the residuals of the discrete equations and the runtime |
| `provenance.toml` | Julia version, BLAS threads, solver description, tolerances, experiment count and date |
| `runs/NAME_config.toml` | The `Config` of the run |
| `runs/NAME_state.csv` | No header; columns `t, R, W_1, ..., W_N` with `W(s,t) = R(t) s z(R(t) s, t)` at the interior nodes `s_i = i/(N+1)` |
| `runs/NAME_control.csv` | No header; columns are the right endpoint `t_n` of each time interval followed by the hat-function values of `q = r u` at `r_j = b j/(K+1)`, held constant on `(t_{n-1}, t_n]` |
| `runs/NAME_history.csv` | Outer-iteration defects, optimality residual, constraint residual, multiplier and control norm |
| `runs/NAME_passive_state.csv` | The uncontrolled trajectory from the same initial data on the same grid |
| `runs/NAME_psi.csv` | The discrete adjoint of the terminal-volume functional, one column per time interval, for `fan_Rf_0.95` and the pure domain-control runs |

The 99 `landscape_*` runs appear only in `summary.csv`.

## Columns of `summary.csv`

`name, N, M, K, amplitude, profile, T, b, R0, Rf, eta, constrained` are the
parameters. `h1` is the gradient seminorm of the initial temperature on the
ball of radius `R0`, and `delta0 = h1 + |R0 - Rf|` is the smallness parameter
of the main theorem. `peak` and `terminal` are the maximal and terminal
spherical `L2` norms of the temperature on the design grid, `relative` divides
the terminal norm by the initial norm (undefined for zero initial temperature)
and `extinction` divides it by `peak`. The `replay_*` columns measure the same
quantities after a fresh nonlinear forward solve with the stored control on the
grid `(2N+1, 2M)`, and `uncontrolled_*` refer to the zero control. `energy` is
the control norm in `L2(B_b x (0,T))`, split into `energy_parallel` (the
multiple of the auxiliary adjoint that adjusts the volume) and
`energy_orthogonal` (volume-neutral part), with
`energy^2 = energy_parallel^2 + energy_orthogonal^2`. `psinorm` is the norm of
the auxiliary adjoint in the control space, `lambda` the multiplier of the
volume constraint, `excursion = max |R(t) - R0|`, `Rlo` and `Rhi` the extreme
radii, `Rfinal` and `radius_error = |R(T) - Rf|` the terminal radius on the
design grid, `replay_Rfinal` and `replay_radius_error` the same after the fine
replay, `radius_h1` the discrete `H1(0,T)` norm of `R - Rf`, `pde` and
`stefan` the maximal absolute residuals of the discrete equations evaluated
with the computed radius, `iterations` the number of outer iterations,
`defect` the final fixed-point defect, `optimality` the relative KKT residual,
`constraint` the residual of the linear volume constraint and `seconds` the
runtime of the run including its replays.

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
