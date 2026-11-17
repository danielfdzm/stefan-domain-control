# Methods and verification

The repository includes the solver, tests, experiment and visualization
scripts, archived records, figures, tables and animation. All commands below
run from the repository root and use the checked-in Julia environment.

## Environment and included records

Use the Julia environment pinned in `Project.toml` and `Manifest.toml`. The
recorded environment is Julia 1.13.0 with CairoMakie 0.15.15 and one BLAS
thread. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
julia --project=. scripts/check.jl
```

The tests cover the spherical weighting of the control, the physical support,
the sign of the Stefan law, manufactured consistency of the ALE scheme, the
discrete forward–adjoint identity, the affine structure of the volume
functional and its gradient, the KKT conditions and the orthogonal
decomposition of the control, the coupled solve with prescribed radius, the
independent replay, pure domain control from zero temperature, the zero-data
invariance, the plain-text record round trip and the continuum normalization
of the auxiliary adjoint at a nonunit radius.

`scripts/check.jl` verifies all files in `data/`, `figures/` and `tables/`
against `checksums.sha256`, checks the recorded tolerances for the 126 paper
runs, independently replays the three animation controls on their design
meshes, and recomputes two small paper experiments. The test suite separately
exercises the solver and its numerical identities. To check only archive
integrity, use `shasum -a 256 -c checksums.sha256` (or
`sha256sum -c checksums.sha256` where available).

## Regenerate figures and tables from the archive

The archived records are sufficient to regenerate the paper artifacts without
rerunning the full experiment suite:

```bash
julia --project=. scripts/make_figures.jl
julia --project=. scripts/make_tables.jl
```

Figures are written to `figures/paper/` as PDFs and `figures/previews/` as PNGs.
Tables are written to `tables/` as LaTeX and CSV files, including
`terminal_accuracy.tex`, `terminal_accuracy.csv` and the all-experiment
`terminal_accuracy_all.csv`. The figure script accepts optional selections:

```bash
julia --project=. scripts/make_figures.jl structure
```

Available figure names are `fan`, `pure`, `balls`, `landscape`, `structure`
and `surfaces`. The plotting scripts may evaluate the frozen control
construction to recover components from an archived trajectory, but do not
rerun the full nonlinear experiment suite. Shared plotting style is in
`scripts/figure_style.jl`.

## Recompute the paper experiments

```bash
julia --project=. scripts/run_experiments.jl
julia --project=. scripts/make_figures.jl
julia --project=. scripts/make_tables.jl
julia --project=. test/runtests.jl
```

The driver runs all 126 experiments with one BLAS thread and overwrites
`data/summary.csv`, `data/provenance.toml` and the corresponding trajectories
in `data/runs/`. It checks a fresh same-grid nonlinear replay for every run.
The 27 runs outside the landscape sweep also receive a doubled-grid replay;
the 99 landscape runs retain only their summary and a same-grid replay.

Optional driver arguments select groups: `fan`, `pure`, `mesh`, `eta`,
`horizon`, `oscillatory` and `landscape`. A selected-group run replaces the
summary and provenance with that subset, so run it in a separate checkout if
you want to preserve the complete archive. The table and figure generators,
and the verification script, expect the complete 126-run summary.

The computation is deterministic up to floating-point and library/platform
differences. Timings and provenance dates change on rerun; PDF metadata and
font rendering can also change artifact bytes. The original checksum manifest
verifies the supplied archive, not byte identity of a newly generated archive.
After deliberately regenerating and reviewing the outputs, refresh and verify
the manifest with:

```bash
julia --project=. scripts/check.jl --write
julia --project=. scripts/check.jl
```

The `--write` command also performs the numerical checks described above.

## Figures, tables and their records

The 126 deterministic experiments are documented in `data/summary.csv`,
`data/provenance.toml` and `data/runs/`. The following table maps each generated
figure and table to its underlying records.

| Figure or table | Records |
| --- | --- |
| `fig_fan.pdf` | `fan_Rf_0.85` … `fan_Rf_1.15`, `free_radius` and the uncontrolled replay |
| `fig_balls.pdf` | `fan_Rf_0.85`, `fan_Rf_1.15` |
| `fig_pure_domain.pdf` | `pure_Rf_0.90`, `pure_Rf_1.10` |
| `fig_structure.pdf` | `fan_Rf_0.95`; the decomposition is evaluated at the converged frozen coefficients |
| `fig_landscape.pdf` | the 99 `landscape_*` rows of `summary.csv` |
| `fig_surfaces.pdf` | `fan_Rf_0.85` |
| `stefan_domain_control.gif` | `animation_Rf_0.50`, `animation_free_radius`, `animation_Rf_1.50` in `data/animation/` |
| `targets.tex` | fan and `free_radius` |
| `terminal_accuracy.tex`, `terminal_accuracy.csv` | fan, `free_radius` and the two `pure_Rf_*` cases |
| `terminal_accuracy_all.csv` | all 126 rows of `summary.csv`, with the replay grid identified |
| `mesh_refinement.tex` | `mesh_31`, `mesh_63`, `mesh_127`, `mesh_255` |
| `regularization.tex` | `eta_-4`, `eta_-6`, `eta_-8`, `mesh_63`, `eta_-12` |
| `horizon_actuator.tex` | `time_*`, `support_*`, `mesh_63` |

In the paper figures, signed temperature and control fields use the
blue–gray–yellow `cividis` palette with limits $[-A,A]$, where $A$ is the largest
absolute value in the fields being compared, so zero is gray. Compared panels
share these limits. Nonnegative diagnostics retain the sequential `Blues`
palette, with ranges shown on their colorbars. The animation uses a separate
colorbar for each evolution, as described below.

## README animation experiment

Regenerate the GIF and still image from the archived animation records:

```bash
julia --project=. scripts/make_animation.jl
```

An optional integer sets the number of animation frames, for example
`julia --project=. scripts/make_animation.jl 90`. To recompute the three
animation experiments first:

```bash
julia --project=. scripts/run_animation_experiments.jl
julia --project=. scripts/make_animation.jl
```

The animation driver updates only `data/animation/`. Rendering writes
`figures/animation/stefan_domain_control.gif` and
`figures/animation/stefan_domain_control_still.png`.


The animation compares prescribed targets $R_f=0.5$ and $1.5$ with a
recomputed temperature-only control on the same mesh $(N,M,K)=(127,512,32)$.
All three runs use $R_0=1$, $T=0.25$, $b=0.4$, $(R_*,E)=(0.45,1.6)$,
$\eta=10^{-10}$ and the same positive initial profile with amplitude $a=0.02$.
The smaller actuator stays inside the contracting ball. The three trajectories
share one spatial scale and the cividis palette. Each evolution $j$ has its
own symmetric colorbar with limits $[-A_j,A_j]$, where $A_j$ is the largest
absolute temperature over that run’s full trajectory. Zero is the gray
midpoint, with blue below zero and yellow above. These independent ranges stay
fixed across frames. Playback lasts 13.6 seconds at the default frame count, including the final-state hold;
playback timing does not change the simulation horizon.

Their parameters, states, controls and iteration histories are stored in
`data/animation/`, alongside a separate `summary.csv` and `provenance.toml`.
The summary records the terminal errors and the independent replays on both
the design mesh and a doubled mesh. These three runs are separate from the
126 paper experiments in `data/summary.csv`. The doubled-mesh animation
replays use an iteration tolerance of $10^{-12}$ to accommodate floating-point
roundoff in the larger temperature pulses; design-mesh replays retain
$2\times10^{-14}$.

| Target | Iterations | Design $R(T)$ | Replay $R(T)$ | Design extinction | Replay extinction |
| --- | --- | --- | --- | --- | --- |
| 0.5 | 38 | 0.500000 | 0.506030 | 3.68e-06 | 6.78e-03 |
| free | 14 | 1.002471 | 1.002475 | 2.85e-06 | 9.27e-04 |
| 1.5 | 25 | 1.500000 | 1.497339 | 3.22e-05 | 6.37e-04 |

Extinction is the terminal spherical temperature norm divided by its peak
value during the run. The animation shows the design-grid trajectories, whose
prescribed terminal radii are reached to rounding error. The doubled-grid
replay errors are about $0.00603$ for contraction and $0.00266$ for expansion;
these larger changes in radius have greater discretization error than the
paper's smaller target changes.

## Model and discretization

With $s=r/R(t)$, $W(s,t)=R(t)\,s\,z(R(t)s,t)$ and $q(r,t)=r\,v(r,t)$, the radial
problem reads

$$
W_t-R^{-2}W_{ss}-s\frac{R'}{R}W_s+Rs\,f\!\left(\frac{W}{Rs}\right)=q(Rs,t),
\qquad W(0,t)=W(1,t)=0,
$$

and the Stefan law in the volume variable $p=R^3/3$ is $p'(t)=-W_s(1,t)$. The
reaction is $f(y)=-4\tanh y$, written as the secant potential
$g(y)=f(y)/y$ with $g(0)=f'(0)$. The interior nodes are $s_i=i/(N+1)$,
centered differences approximate the spatial derivatives, backward Euler
advances the state, and the flux $W_s(1,t_n)$ is the second-order one-sided
difference $(W_{N-1}^n-4W_N^n)/(2h)$. The discrete Stefan law is
$p^n=p^{n-1}-\Delta t\,W_s(1,t_n)$ and $R^n=(3p^n)^{1/3}$; the transformed
control is piecewise affine in the physical radius on $[0,b]$ and piecewise
constant in time, so the actuator stays the fixed ball $B_b$. The control mass
matrix includes the factor $4\pi$, so the discrete control norm is the norm in
$L^2(B_b\times(0,T))$, and terminal norms include the factor $4\pi R(T)$.

The numerical control basis covers the physical ball over the whole time
interval. The proof uses a smaller moving annulus and starts the control after
an initial uncontrolled evolution. Both supports are admissible for the same
physical controllability problem; the computations do not reproduce that
particular support restriction from the proof.

## Constrained HUM for the frozen problem

For frozen state and radius, the terminal map is $x\mapsto Hx+d$ in the
weighted control variables $x^n=\sqrt{\Delta t}\,U_qC^n$, and the terminal
volume functional is affine, $x\mapsto g^\top x+\ell_0$. Both are assembled with
the exact discrete adjoints. The vector $g$ is the gradient of the integrated
boundary flux in the weighted physical-control coordinates. Its Euclidean norm
$|g|$ is the norm of this discrete Riesz representative.

The state recursion uses the normalized operator $\mathcal N=R^{-2}\mathcal M$,
where $\mathcal M$ is the operator in the proof. Hence
$\mathcal M^*\psi=\mathcal N^*(R^2\psi)$: the raw discrete adjoint $\Lambda$
satisfies $-\Lambda/h\approx R(t)^2\psi$, and the comparison with the
unit-boundary auxiliary adjoint is

$$
\psi(s_i,t_n)\approx-\frac{\Lambda_i^n}{hR(t_n)^2}.
$$

The adjoint figure uses this normalization and omits the final spatial node,
where the one-sided flux stencil creates a discrete source artifact. All
remaining displayed values determine its color range. Neither the raw adjoint
nor the flux gradient used by the control solver is altered for plotting.
The physical control metric and the proof's transformed $L^2$ metric differ,
so $g$ should not be identified directly with sampled values of $\psi$.
The control solves

$$
\min_x\ \tfrac12|x|^2+\tfrac{1}{2\eta}|Hx+d|^2
\quad\text{subject to}\quad g^\top x=m,
\qquad m=\frac{R_0^3-R_f^3}{3}-\ell_0 .
$$

With the thin SVD $H=U\Sigma V^\top$ the unconstrained minimizer is
$x_0=-V\,\mathrm{diag}\!\big(\sigma_j/(\sigma_j^2+\eta)\big)U^\top d$, the
constraint direction is $z=(H^\top H+\eta I)^{-1}g$, and
$x=x_0+\lambda z$ with $\lambda=(m-g^\top x_0)/(g^\top z)$. Writing
$x_\parallel=(g^\top x/|g|^2)\,g$ and $x_\perp=x-x_\parallel$ reproduces the
orthogonal decomposition in the physical control metric: $x_\parallel$ is
the multiple of $g$ that adjusts the volume and $x_\perp$ leaves it unchanged.
The quantities
`energy_parallel`, `energy_orthogonal`, `psinorm` $=|g|$ and `lambda` in
`summary.csv` are these objects. In the free-radius mode the constraint is
dropped and $\lambda=0$.

## Fixed-point iteration and replay

The outer iteration freezes the state and the radius, solves the constrained
problem, integrates the discrete Stefan law from the computed flux, and
relaxes with factor 0.8. The initial frozen radius interpolates linearly from
$R_0$ to $R_f$. The iteration stops when the maximum of the state, radius and
radius-speed defects divided by $\delta_0$ is below $10^{-9}$; every radius
iterate must stay in its configured admissible interval $(R_*,E)$, and lack
of convergence is reported as failure. At convergence the computed radius satisfies the discrete Stefan
law exactly and $R(T)=R_f$ to rounding, and the residuals of the nonlinear
discrete equations are checked by substitution.

Each control is then replayed by a fresh coupled nonlinear forward solve on the
design grid; the 27 runs outside the landscape sweep are also replayed on
$(2N+1,2M)$. Each replay solves state and Stefan law together at every step
to a tolerance of $2\times10^{-14}$ for the paper experiments. The replay
imposes no terminal condition: its terminal temperature and radius are measured. The same-grid
replay agrees with the optimized trajectory to $10^{-9}$ in every stored run.

## Terminal accuracy

For each trajectory, the reported temperature norm is the spherical $L^2$ norm
on its current ball. The extinction ratio divides the terminal norm by the
peak norm over the trajectory. Because moving the interface can create a large
temperature pulse, this ratio need not be small relative to the initial norm.

`tables/terminal_accuracy.tex` and `.csv` report the absolute initial,
design-terminal and doubled-grid replay-terminal norms for the seven targets,
the free-radius reference and the two zero-initial-temperature cases. They
also report the replay terminal norm divided by its replay-grid initial norm;
this ratio is marked undefined for zero initial temperature. The design and
replay initial norms use their respective grids. The same diagnostics for all
126 experiments are in `tables/terminal_accuracy_all.csv`, whose `replay_grid`
column distinguishes the landscape's same-grid replay from a doubled-grid
replay. These norms complement the peak-normalized extinction ratios in the
other tables.

## Experiments

The 126 paper runs use $R_0=1$, $R_*=0.8$, $E=1.2$ and, unless stated otherwise,
$T=0.25$, $b=0.55$, $\eta=10^{-10}$ and the positive profile
$z_0(r)=a(1-s^2)(1+\tfrac14\cos 2\pi s)$ with $a=0.02$. The oscillatory profile
is $a(1-s^2)\cos 2\pi s$.

- **fan**: $R_f\in\{0.85,0.90,\dots,1.15\}$ on $(N,M,K)=(127,512,32)$, plus the
  temperature-only control with free terminal radius.
- **pure**: $a=0$ and $R_f\in\{0.90,1.10\}$ on the fine mesh.
- **mesh**: $(31,128,8)$, $(63,256,16)$, $(127,512,32)$, $(255,1024,64)$ with
  $R_f=0.95$.
- **eta**: $\eta\in\{10^{-4},10^{-6},10^{-8},10^{-12}\}$ with $R_f=0.95$ on
  $(63,256,16)$; the row $\eta=10^{-10}$ is `mesh_63`.
- **horizon**: $T\in\{0.15,0.20,0.35\}$ and $b\in\{0.40,0.475,0.65\}$ with
  $R_f=0.95$ on $(63,256,16)$.
- **oscillatory**: $R_f\in\{0.97,1.00,1.03\}$ on $(63,256,16)$.
- **landscape**: $a\in\{0,0.005,\dots,0.04\}$ and
  $R_f-R_0\in\{-0.05,-0.04,\dots,0.05\}$ on $(63,256,16)$, 99 runs without
  stored trajectories.

Every parameter is set through the `Config` keyword constructor in
`src/StefanDomainControl.jl`. The recorded parameters are in `data/summary.csv`
and, for runs with stored trajectories, `data/runs/NAME_config.toml`.
