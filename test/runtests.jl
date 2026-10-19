using StefanDomainControl, Test, LinearAlgebra, Random
BLAS.set_num_threads(1)

@testset "Geometry, weighting and physical support" begin
    p = Config(N=31, M=64, K=8)
    c = collect(1.0:p.K)
    rs = ((1:100000) .- 0.5) .* p.b/100000
    quadrature = 4pi*p.b/length(rs)*sum(control_value(r, c, p)^2 for r in rs)
    @test quadrature ≈ dot(c, mass_control(p)*c) rtol=1e-8
    @test control_value(p.b+0.01, c, p) == 0
    @test control_value(0.0, c, p) == 0
    @test control_value(1e-6, c, p)/1e-6 ≈ control_value(2e-6, c, p)/2e-6
    s = grid(p)
    w = s .* (1 .- s)
    @test flux(w, p) ≈ -1
    @test dot(flux_stencil(p), w) ≈ flux(w, p)
    @test -flux(w, p) > 0            # positive temperature expands the ball
    @test radius(volume(0.93)) ≈ 0.93
    @test delta0(Config(amplitude=0.0, Rf=0.97)) ≈ 0.03
end

@testset "Manufactured ALE consistency" begin
    errors = Float64[]
    for N in (31, 63, 127)
        p = Config(N=N, M=(N+1)^2, reaction=0.0)
        h = 1/(N+1); dt = p.T/p.M; t = 0.1; R = 1+0.2t
        s = collect(1:N) .* h
        w = exp(-t) .* s .* (1 .- s.^2)
        prev = exp(-(t-dt)) .* s .* (1 .- s.^2)
        A, _ = operators(w, R, 0.2, p)
        exact = -w + 6exp(-t) .* s/R^2 - (0.2/R) .* s .* exp(-t) .* (1 .- 3s.^2)
        push!(errors, maximum(abs, (A*w-prev)/dt-exact))
    end
    @test errors[1]/errors[2] > 3.8
    @test errors[2]/errors[3] > 3.8
end

@testset "Discrete forward-adjoint identity" begin
    Random.seed!(47)
    p = Config(N=15, M=8, K=5)
    F = []; B = []
    for n in 1:p.M
        A, b = operators(0.01randn(p.N), 1+0.001n, 0.02, p)
        push!(F, lu(A)); push!(B, b)
    end
    v = randn(p.K, p.M); terminal = randn(p.N); w = zeros(p.N)
    for n in 1:p.M
        w = F[n] \ (w + (p.T/p.M)*B[n]*v[:,n])
    end
    λ = copy(terminal); pairing = 0.0
    for n in p.M:-1:1
        λ = transpose(F[n]) \ λ
        pairing += (p.T/p.M)*dot(B[n]'*λ, v[:,n])
    end
    @test dot(w, terminal) ≈ pairing rtol=2e-13
end

function volume_functional(x, Wbar, Rbar, p)
    dt = p.T/p.M; Uc = cholesky(mass_control(p)).U
    C = (Uc \ reshape(x, p.K, p.M)) ./ sqrt(dt)
    w = initial_state(p); l = 0.0
    for n in 1:p.M
        A, B = operators(Wbar[:,n+1], Rbar[n+1], (Rbar[n+1]-Rbar[n])/dt, p)
        w = A \ (w + dt*B*C[:,n]); l += dt*flux(w, p)
    end
    l
end

@testset "Volume constraint, auxiliary adjoint and KKT conditions" begin
    Random.seed!(3)
    p = Config(N=15, M=12, K=4, Rf=0.97)
    Wbar = 0.01 .* randn(p.N, p.M+1); Rbar = collect(range(1.0, 0.98, length=p.M+1))
    W, vol, C, psi, info = frozen_control(Wbar, Rbar, p)
    x = randn(p.K*p.M)
    # The volume functional is affine in the control with gradient given by psi.
    @test volume_functional(x, Wbar, Rbar, p) ≈ dot(info.gvec, x) + info.l0 rtol=1e-11
    @test volume_functional(zeros(p.K*p.M), Wbar, Rbar, p) ≈ info.l0 rtol=1e-12
    @test size(psi) == (p.N, p.M)
    # The prescribed terminal volume holds exactly for the frozen problem.
    @test vol[end] ≈ volume(p.Rf) atol=1e-14
    @test info.constraint < 1e-13
    @test info.optimality < 1e-10
    # Decomposition into the psi-multiple and its orthogonal complement.
    parallel = (dot(info.gvec, info.x)/info.psinorm^2) .* info.gvec
    @test norm(parallel) ≈ abs(info.m)/info.psinorm
    @test info.energy^2 ≈ info.energy_parallel^2 + info.energy_orthogonal^2
    # Without the constraint the multiplier vanishes and the volume is not prescribed.
    Wf, volf, Cf, psif, infof = frozen_control(Wbar, Rbar, Config(N=15, M=12, K=4, Rf=0.97, constrained=false))
    @test infof.lambda == 0
    @test abs(volf[end] - volume(p.Rf)) > 1e-6
end

@testset "Coupled control with prescribed radius and independent replay" begin
    r = solve_control(Config(N=31, M=128, K=8, Rf=0.95); verbose=false)
    d = diagnostics(r)
    @test d.radius_error < 1e-12
    @test d.terminal/maximum(d.norms) < 1e-3
    @test d.pde < 2e-7
    @test d.stefan < 2e-7
    @test r.info.optimality < 1e-8
    f = replay(r)
    @test maximum(abs, f.W - r.W) < 1e-9
    @test maximum(abs, f.R - r.R) < 1e-9
    @test abs(f.R[end] - 0.95) < 1e-9
    passive = replay(r; uncontrolled=true)
    @test abs(passive.R[end] - 0.95) > 1e-2
    fine = replay(r; N=63, M=256)
    @test abs(fine.R[end] - 0.95) < 1e-3
    # Pure domain control from zero temperature.
    z = solve_control(Config(N=31, M=128, K=8, Rf=1.03, amplitude=0.0); verbose=false)
    dz = diagnostics(z)
    @test dz.radius_error < 1e-12
    @test maximum(dz.norms) > 0
    @test dz.terminal/maximum(dz.norms) < 1e-3
    # Zero data and matched radius give the zero control and a resting ball.
    p0 = Config(N=15, M=16, K=4, amplitude=0.0, Rf=1.0)
    W0, vol0, C0, _, info0 = frozen_control(zeros(p0.N, p0.M+1), fill(1.0, p0.M+1), p0)
    @test all(iszero, C0)
    @test all(==(volume(1.0)), vol0)
    zr = replay((p=p0, C=zeros(p0.K, p0.M)))
    @test all(iszero, zr.W)
    @test all(==(p0.R0), zr.R)
end

@testset "Plain-text records round trip" begin
    dir = mktempdir()
    r = solve_control(Config(N=15, M=16, K=4, Rf=0.98); verbose=false)
    save_run(dir, "tiny", r; psi=true, passive=replay(r; uncontrolled=true))
    back = load_run(dir, "tiny")
    @test back.p == r.p
    @test back.W ≈ r.W
    @test back.R ≈ r.R
    @test back.C ≈ r.C
    @test back.psi ≈ r.psi
    @test haskey(back, :passive)
end
