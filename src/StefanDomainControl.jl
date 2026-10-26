"""
    StefanDomainControl

Numerical control of the three-dimensional radial one-phase Stefan problem with a
prescribed terminal radius. The moving ball is mapped onto the fixed interval
(0,1) through `W(s,t) = R(t) s z(R(t) s, t)` and the Stefan law is written in the
volume variable `p = R^3/3`, for which `p' = -W_s(1,t)`. For frozen coefficients the
prescribed terminal volume is one linear constraint on the control, and the
discrete control is the minimizer of a Tikhonov-regularized HUM functional under
this constraint. The outer fixed-point iteration updates the frozen interface and
the frozen semilinearity until the nonlinear discrete equations hold.
"""
module StefanDomainControl

using LinearAlgebra, Printf, DelimitedFiles, TOML

export Config, solve_control, frozen_control, replay, diagnostics, initial_profile,
       initial_state, physical_state, physical_control, control_value, control_matrix,
       operators, mass_control, flux, flux_stencil, grid, volume, radius, h1_seminorm,
       delta0, save_run, load_run, config_dict, config_from_dict

"""
    Config

Parameters of one experiment. `N` counts interior state nodes on (0,1), `M` time
intervals and `K` interior control nodes on the physical interval (0,b). When
`constrained` is true the terminal radius `Rf` is prescribed exactly; otherwise only
the terminal temperature is penalized and the terminal radius is free.
"""
Base.@kwdef struct Config
    N::Int = 63
    M::Int = 256
    K::Int = 16
    T::Float64 = 0.25
    b::Float64 = 0.55
    R0::Float64 = 1.0
    Rf::Float64 = 0.95
    Rmin::Float64 = 0.8
    Rmax::Float64 = 1.2
    amplitude::Float64 = 0.02
    profile::Symbol = :positive
    reaction::Float64 = -4.0
    eta::Float64 = 1e-10
    constrained::Bool = true
end

volume(R) = R^3/3
radius(v) = cbrt(3v)
grid(p::Config) = collect(1:p.N) ./ (p.N+1)

"Radial initial temperature z_0 as a function of s = r/R_0."
function initial_profile(s, p::Config)
    shape = p.profile == :positive ? 1 + 0.25*cospi(2s) :
            p.profile == :oscillatory ? cospi(2s) :
            error("Unknown profile $(p.profile)")
    p.amplitude*(1-s^2)*shape
end

"Derivative of the initial profile with respect to s."
function initial_slope(s, p::Config)
    if p.profile == :positive
        p.amplitude*(-2s*(1+0.25cospi(2s)) - 0.5pi*(1-s^2)*sinpi(2s))
    else
        p.amplitude*(-2s*cospi(2s) - 2pi*(1-s^2)*sinpi(2s))
    end
end

"Transformed initial state W_0(s) = R_0 s z_0(R_0 s) at the interior nodes."
initial_state(p::Config) = p.R0 .* grid(p) .* initial_profile.(grid(p), Ref(p))

"Gradient seminorm of y_0 on the ball B_{R_0}, by fine midpoint quadrature."
function h1_seminorm(p::Config; n=20000)
    ss = ((1:n) .- 0.5) ./ n
    sqrt(4pi*p.R0*sum(ss.^2 .* initial_slope.(ss, Ref(p)).^2)/n)
end

"Smallness parameter of the main theorem, delta_0 = |y_0|_{H^1_0} + |R_0 - R_f|."
delta0(p::Config) = h1_seminorm(p) + abs(p.R0-p.Rf)

"Secant potential g(y) = f(y)/y for f(y) = reaction*tanh(y), with g(0) = f'(0)."
secant(y, p::Config) = abs(y) < 1e-8 ? p.reaction*(1-y^2/3) : p.reaction*tanh(y)/y

"Second-order one-sided approximation of W_s(1,t), using W(1,t) = 0."
flux(w, p::Config) = (w[end-1]-4w[end])*(p.N+1)/2

"Vector e with flux(w,p) == dot(e,w)."
function flux_stencil(p::Config)
    e = zeros(p.N); e[end-1] = (p.N+1)/2; e[end] = -2(p.N+1); e
end

"Exact mass matrix of the physical hat functions, including the factor 4pi."
function mass_control(p::Config)
    d = p.b/(p.K+1)
    SymTridiagonal(fill(4pi*2d/3, p.K), fill(4pi*d/6, p.K-1))
end

"q(r) = r*v(r), piecewise affine on [0,b] with q(0) = q(b) = 0, and zero outside."
function control_value(r, c, p::Config)
    (r <= 0 || r >= p.b) && return 0.0
    ξ = r*(p.K+1)/p.b
    j = floor(Int, ξ)
    θ = ξ-j
    left = j == 0 ? 0.0 : c[j]
    right = j == p.K ? 0.0 : c[j+1]
    (1-θ)*left + θ*right
end

"Evaluation of the physical control basis at the moving nodes r = R s_i."
function control_matrix(R, p::Config)
    B = zeros(p.N, p.K)
    for (i,s) in enumerate(grid(p))
        r = R*s
        r >= p.b && continue
        ξ = r*(p.K+1)/p.b
        j = floor(Int, ξ)
        θ = ξ-j
        j > 0 && (B[i,j] = 1-θ)
        j < p.K && (B[i,j+1] = θ)
    end
    B
end

"Backward Euler matrix of the ALE equation for frozen radius, mesh velocity and potential."
function operators(wbar, R, speed, p::Config)
    h = 1/(p.N+1); dt = p.T/p.M; s = grid(p)
    diffusion = 1/(R*h)^2
    a = secant.(wbar ./ (R.*s), Ref(p))
    adv = s .* (speed/R)/(2h)
    A = Tridiagonal(-dt .* (diffusion .- adv[2:end]),
                    1 .+ dt .* (2diffusion .+ a),
                    -dt .* (diffusion .+ adv[1:end-1]))
    A, control_matrix(R, p)
end

"""
    frozen_control(Wbar, Rbar, p)

Constrained regularized HUM for the linear problem with frozen state `Wbar` and
frozen radius `Rbar`. The terminal map `x -> H x + d` and the volume functional
`x -> g'x + l0` are assembled with the exact discrete adjoints. The control
minimizes `|x|^2/2 + |Hx+d|^2/(2 eta)` subject to `g'x = m`, where `m` is the
flux integral that steers the volume to `R_f^3/3`. The multiplier `lambda` and the
discrete auxiliary adjoint `psi` are returned together with the trajectory.
"""
function frozen_control(Wbar, Rbar, p::Config)
    dt = p.T/p.M; h = 1/(p.N+1)
    Uc = cholesky(mass_control(p)).U
    factors = Vector{Any}(undef, p.M)
    Bs = Vector{Matrix{Float64}}(undef, p.M)
    W0 = initial_state(p)
    d = copy(W0); l0 = 0.0
    for n in 1:p.M
        A, B = operators(Wbar[:,n+1], Rbar[n+1], (Rbar[n+1]-Rbar[n])/dt, p)
        factors[n] = lu(A); Bs[n] = B
        d = factors[n] \ d
        l0 += dt*flux(d, p)
    end
    P = Matrix{Float64}(I, p.N, p.N)
    H = zeros(p.N, p.K*p.M)
    gvec = zeros(p.K*p.M)
    psi = zeros(p.N, p.M)
    e = flux_stencil(p)
    lam = zeros(p.N)
    for n in p.M:-1:1
        # Right division implements P*A_n^{-1} without forming an inverse.
        P = Matrix(transpose(transpose(factors[n]) \ Matrix(transpose(P))))
        H[:,(n-1)*p.K+1:n*p.K] = sqrt(dt) .* (P*Bs[n]/Uc)
        # Discrete auxiliary adjoint: the flux functional enters as a source.
        lam = transpose(factors[n]) \ (dt .* e .+ lam)
        psi[:,n] = lam
        gvec[(n-1)*p.K+1:n*p.K] = sqrt(dt) .* (Uc' \ (Bs[n]'*lam))
    end
    terminal_weight = sqrt(4pi*Rbar[end]*h)
    H .*= terminal_weight; d .*= terminal_weight
    S = svd(H; full=false)
    shrink = S.S ./ (S.S.^2 .+ p.eta)
    x0 = -S.V * (shrink .* (S.U' * d))
    if p.constrained
        m = (volume(p.R0) - volume(p.Rf)) - l0
        Vg = S.V' * gvec
        z = S.V * (Vg ./ (S.S.^2 .+ p.eta)) .+ (gvec .- S.V*Vg) ./ p.eta
        lambda = (m - dot(gvec, x0)) / dot(gvec, z)
        x = x0 .+ lambda .* z
    else
        m = NaN; lambda = 0.0; x = x0
    end
    C = (Uc \ reshape(x, p.K, p.M)) ./ sqrt(dt)
    W = zeros(p.N, p.M+1); W[:,1] = W0
    vol = zeros(p.M+1); vol[1] = volume(p.R0)
    for n in 1:p.M
        W[:,n+1] = factors[n] \ (W[:,n] + dt*Bs[n]*C[:,n])
        vol[n+1] = vol[n] - dt*flux(W[:,n+1], p)
    end
    gradient = p.eta .* x .+ H'*(H*x .+ d) .- lambda .* gvec
    optimality = norm(gradient) / max(norm(H'*d), abs(lambda)*norm(gvec), eps())
    constraint = p.constrained ? abs(dot(gvec, x) - m) : 0.0
    gnorm = norm(gvec)
    parallel = (dot(gvec, x)/gnorm^2) .* gvec
    info = (lambda=lambda, m=m, l0=l0, optimality=optimality, constraint=constraint,
            psinorm=gnorm, energy=norm(x), energy_parallel=norm(parallel),
            energy_orthogonal=norm(x .- parallel), gvec=gvec, x=x)
    W, vol, C, psi, info
end

"""
    solve_control(p; maxiter=100, tol=1e-9, damping=0.8, verbose=true)

Relaxed fixed-point iteration on the frozen state and radius. The iteration stops
when the state, radius and radius-speed updates, divided by `delta0(p)`, are below
`tol`. Every radius iterate must stay in `(Rmin, Rmax)`; failure is an error.
"""
function solve_control(p::Config; maxiter=100, tol=1e-9, damping=0.8, verbose=true)
    @assert 0 < p.b < p.Rmin < p.R0 < p.Rmax
    @assert p.Rmin < p.Rf < p.Rmax
    dt = p.T/p.M
    W = repeat(initial_state(p), 1, p.M+1)
    R = p.constrained ? collect(range(p.R0, p.Rf, length=p.M+1)) : fill(p.R0, p.M+1)
    scale = max(delta0(p), 1e-12)
    history = NamedTuple[]
    for k in 1:maxiter
        V, vol, C, psi, info = frozen_control(W, R, p)
        Q = radius.(vol)
        ew = maximum(abs, V-W)/scale
        er = maximum(abs, Q-R)/scale
        ev = maximum(abs, diff(Q-R))/dt/scale
        defect = max(ew, er, ev)
        push!(history, (iteration=k, defect=defect, state=ew, radius=er, speed=ev,
                        optimality=info.optimality, constraint=info.constraint,
                        lambda=info.lambda, energy=info.energy))
        verbose && @printf("N=%d Rf=%.3f a=%.3g k=%02d defect=%.3e terminal=%.3e R(T)-Rf=%.2e\n",
            p.N, p.Rf, p.amplitude, k, defect,
            sqrt(4pi*Q[end]/(p.N+1))*norm(V[:,end]), Q[end]-p.Rf)
        all(p.Rmin .< Q .< p.Rmax) || error("Radius left the admissible interval")
        if defect < tol
            result = (p=p, W=V, R=Q, vol=vol, C=C, psi=psi, info=info, history=history)
            check = diagnostics(result)
            check.pde < 2e-7 || error("Nonlinear PDE residual too large: $(check.pde)")
            check.stefan < 2e-7 || error("Stefan residual too large: $(check.stefan)")
            return result
        end
        W = (1-damping).*W + damping.*V
        R = (1-damping).*R + damping.*Q
    end
    error("Fixed-point iteration did not converge; no result is reported")
end

"""
    replay(result; N, M, uncontrolled=false, tol=2e-14)

Independent nonlinear forward integration of the coupled discrete system with the
stored control, on the design grid or on a finer one. No terminal condition is
imposed, overwritten or projected: the terminal temperature and radius are measured.
`tol` is the absolute nonlinear iteration tolerance. Larger-amplitude states may
need a slightly larger value to stay above floating-point roundoff.
"""
function replay(result; N=result.p.N, M=result.p.M, uncontrolled=false, tol::Real=2e-14)
    isfinite(tol) && tol > 0 || throw(ArgumentError("Replay tolerance must be finite and positive"))
    old = result.p
    @assert M % old.M == 0 "Replay time grid must resolve every control discontinuity"
    p = Config(; (name=>getfield(old,name) for name in fieldnames(Config) if name ∉ (:N,:M))..., N=N, M=M)
    W = zeros(N, M+1); R = fill(p.R0, M+1); vol = fill(volume(p.R0), M+1); dt = p.T/M
    W[:,1] = initial_state(p)
    C = zeros(p.K, M)
    for n in 1:M
        c = uncontrolled ? zeros(p.K) : result.C[:,cld(n, M÷old.M)]
        C[:,n] = c
        w = copy(W[:,n]); rad = R[n]
        converged = false
        for k in 1:200
            A, B = operators(w, rad, (rad-R[n])/dt, p)
            v = A \ (W[:,n] + dt*B*c)
            qv = vol[n] - dt*flux(v, p)
            q = radius(qv)
            defect = max(maximum(abs, v-w), abs(q-rad))
            if defect < tol
                W[:,n+1] = v; R[n+1] = q; vol[n+1] = qv
                converged = true; break
            end
            w = 0.2w + 0.8v; rad = 0.2rad + 0.8q
        end
        converged || error("Forward nonlinear solve failed at step $n")
        p.Rmin < R[n+1] < p.Rmax || error("Replay radius out of bounds")
    end
    (p=p, W=W, R=R, vol=vol, C=C, history=NamedTuple[])
end

"Discrete H^1(0,T) norm of R - R_f with trapezoidal weights."
function radius_h1(R, p::Config)
    dt = p.T/p.M
    w = fill(dt, length(R)); w[1] /= 2; w[end] /= 2
    sqrt(sum(w .* (R .- p.Rf).^2) + sum(diff(R).^2)/dt)
end

"""
    diagnostics(result)

Spherical norms of the temperature, control cost, boundary excursion, terminal
radius error, and absolute residuals of the discrete equations evaluated with the
actual radius trajectory.
"""
function diagnostics(result)
    (;p, W, R, C) = result
    h = 1/(p.N+1); dt = p.T/p.M
    norms = [sqrt(4pi*R[n]*h)*norm(W[:,n]) for n in 1:p.M+1]
    energy = sqrt(dt*sum(C.*(mass_control(p)*C)))
    vol = volume.(R)
    stefan = 0.0; pde = 0.0
    for n in 1:p.M
        speed = (R[n+1]-R[n])/dt
        A, B = operators(W[:,n+1], R[n+1], speed, p)
        res = (A*W[:,n+1]-W[:,n])/dt - B*C[:,n]
        pde = max(pde, maximum(abs, res))
        stefan = max(stefan, abs((vol[n+1]-vol[n])/dt + flux(W[:,n+1], p)))
    end
    initial = norms[1] > 0 ? norms[1] : NaN
    (norms=norms, relative=norms[end]/initial, terminal=norms[end], energy=energy,
     excursion=maximum(abs, R.-p.R0), Rfinal=R[end], radius_error=abs(R[end]-p.Rf),
     Rlo=minimum(R), Rhi=maximum(R), radius_h1=radius_h1(R, p),
     pde=pde, stefan=stefan, h1=h1_seminorm(p), delta0=delta0(p))
end

"Temperature y = W/(R s) at all s-nodes including the endpoints."
function physical_state(result)
    (;p, W, R) = result
    Y = zeros(p.N+2, p.M+1)
    for n in 1:p.M+1
        Y[2:end-1,n] = W[:,n] ./ (R[n].*grid(p))
        # Odd quadratic-consistent extension of W at the origin: y(0) = W_s(0)/R.
        Y[1,n] = (8W[1,n]-W[2,n])*(p.N+1)/(6R[n])
    end
    Y
end

"Physical control u = q/r on a uniform grid of the actuator [0,b], one column per time interval."
function physical_control(result; nr=140)
    p = result.p; rr = collect(range(0, p.b, length=nr))
    V = [x == 0 ? result.C[1,n]*(p.K+1)/p.b : control_value(x, result.C[:,n], p)/x
         for x in rr, n in 1:size(result.C,2)]
    rr, V
end

config_dict(p::Config) = Dict{String,Any}(string(f) =>
    (getfield(p,f) isa Symbol ? string(getfield(p,f)) : getfield(p,f)) for f in fieldnames(Config))

function config_from_dict(d)
    kw = Dict{Symbol,Any}()
    for f in fieldnames(Config)
        v = d[string(f)]
        kw[f] = fieldtype(Config, f) == Symbol ? Symbol(v) : convert(fieldtype(Config, f), v)
    end
    Config(; kw...)
end

"""
    save_run(dir, name, result; psi=false, passive=nothing)

Portable plain-text records: `NAME_config.toml`, `NAME_state.csv` with columns
`t, R, W_1, ..., W_N`, `NAME_control.csv` with the interval right endpoint followed by
the hat-function values of `q`, `NAME_history.csv`, optionally `NAME_psi.csv` with the
discrete auxiliary adjoint, and `NAME_passive_state.csv` for the uncontrolled replay.
"""
function save_run(dir, name, result; psi=false, passive=nothing)
    mkpath(dir)
    p = result.p
    open(io -> TOML.print(io, config_dict(p)), joinpath(dir, name*"_config.toml"), "w")
    t = collect(range(0, p.T, length=p.M+1))
    writedlm(joinpath(dir, name*"_state.csv"), [t result.R permutedims(result.W)], ',')
    writedlm(joinpath(dir, name*"_control.csv"), [collect(1:p.M).*(p.T/p.M) permutedims(result.C)], ',')
    if !isempty(result.history)
        open(joinpath(dir, name*"_history.csv"), "w") do io
            println(io, join(string.(keys(first(result.history))), ","))
            for hrow in result.history; println(io, join(values(hrow), ",")); end
        end
    end
    psi && haskey(result, :psi) && writedlm(joinpath(dir, name*"_psi.csv"),
        [collect(1:p.M).*(p.T/p.M) permutedims(result.psi)], ',')
    passive === nothing || writedlm(joinpath(dir, name*"_passive_state.csv"),
        [t passive.R permutedims(passive.W)], ',')
    nothing
end

"Load a record written by `save_run`; the uncontrolled replay and psi are included when present."
function load_run(dir, name)
    p = config_from_dict(TOML.parsefile(joinpath(dir, name*"_config.toml")))
    state = readdlm(joinpath(dir, name*"_state.csv"), ',', Float64)
    ctrl = readdlm(joinpath(dir, name*"_control.csv"), ',', Float64)
    R = state[:,2]; W = permutedims(state[:,3:end]); C = permutedims(ctrl[:,2:end])
    out = Dict{Symbol,Any}(:p=>p, :W=>W, :R=>R, :vol=>volume.(R), :C=>C, :history=>NamedTuple[])
    psipath = joinpath(dir, name*"_psi.csv")
    isfile(psipath) && (out[:psi] = permutedims(readdlm(psipath, ',', Float64)[:,2:end]))
    passivepath = joinpath(dir, name*"_passive_state.csv")
    if isfile(passivepath)
        ps = readdlm(passivepath, ',', Float64)
        out[:passive] = (p=p, W=permutedims(ps[:,3:end]), R=ps[:,2], vol=volume.(ps[:,2]),
                         C=zeros(p.K, p.M), history=NamedTuple[])
    end
    (; out...)
end

end
