# Methods and verification

The public repository includes the solver, tests, experiment records and
generated figures, tables and animation. The scripts that automate the
experiments and generate the visualizations are maintained privately and are
not included. The methods and parameters below describe the published results.

## Environment and included records

Use the Julia environment pinned in `Project.toml` and `Manifest.toml`. The
recorded environment is Julia 1.13.0 with CairoMakie 0.15.15 and one BLAS
thread. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
shasum -a 256 -c checksums.sha256
```

The tests cover the spherical weighting of the control, the physical support,
the sign of the Stefan law, manufactured consistency of the ALE scheme, the
discrete forward–adjoint identity, the affine structure of the volume
functional and its gradient, the KKT conditions and the orthogonal
decomposition of the control, the coupled solve with prescribed radius, the
independent replay, pure domain control from zero temperature, the zero-data
invariance and the plain-text record round trip.

The checksum command verifies the included files in `data/`, `figures/` and
`tables/` against `checksums.sha256`. It checks file integrity; the Julia tests
independently exercise the solver and its numerical identities.

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

The animation compares prescribed targets $R_f=0.5$ and $1.5$ with a
recomputed temperature-only control on the same mesh $(N,M,K)=(127,512,32)$.
All three runs use $R_0=1$, $T=0.25$, $b=0.4$, $(R_*,E)=(0.45,1.6)$,
$\eta=10^{-10}$ and the same positive initial profile with amplitude $a=0.02$.
The smaller actuator stays inside the contracting ball. The three trajectories
share one spatial scale and the cividis palette. Each evolution has its own
colorbar, whose limits are the minimum and maximum temperatures over that
run’s full trajectory. These three ranges stay fixed across frames. Playback
lasts 13.6 seconds, including the final-state hold; this is half the previous
playback speed and does not change the simulation horizon.

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

## Constrained HUM for the frozen problem

For frozen state and radius, the terminal map is $x\mapsto Hx+d$ in the
weighted control variables $x^n=\sqrt{\Delta t}\,U_qC^n$, and the terminal
volume functional is affine, $x\mapsto g^\top x+\ell_0$. Both are assembled with
the exact discrete adjoints; the vector $g$ is the discrete adjoint of the flux
functional, that is, the auxiliary adjoint $\psi$ restricted to the actuator in
the control metric, and equals $-h$ times the adjoint with unit boundary value.
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
decomposition of the paper: $x_\parallel$ is the multiple of $\psi$ that
adjusts the volume and $x_\perp$ leaves it unchanged. The quantities
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

Each control is then replayed by a fresh coupled nonlinear forward solve, on the
design grid and on $(2N+1,2M)$, solving state and Stefan law together at every
step to a tolerance of $2\times10^{-14}$ for the paper experiments. The
replay imposes no terminal
condition: its terminal temperature and radius are measured. The same-grid
replay agrees with the optimized trajectory to $10^{-9}$ in every stored run.

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
