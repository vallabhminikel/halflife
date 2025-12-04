"""
Mixture Decay Model - Using DifferentialEquations.jl

This implementation uses the DifferentialEquations.jl ecosystem for more
robust ODE solving with adaptive timesteps and better numerical stability.
"""

using DifferentialEquations, Optim, CairoMakie, LaTeXStrings, Interpolations

# ============================================================================
# ODE System Definitions
# ============================================================================

"""
    protein_dynamics_single!(du, u, p, t)

Single-rate protein dynamics ODE system:
    dP/dt = λ*R(t) - λ*P(t)

Arguments:
- u: state vector [P]
- p: parameters (λ, R_interp)
- t: time
"""
function protein_dynamics_single!(du, u, p, t)
    P = u[1]
    λ, R_interp = p

    R_t = R_interp(t)
    du[1] = λ * R_t - λ * P

    return nothing
end

"""
    protein_dynamics_mixture!(du, u, p, t)

Mixture model protein dynamics ODE system:
    dPₖ/dt = λₖ*R(t) - λₖ*Pₖ(t),  k=1,...,K

Arguments:
- u: state vector [P₁, P₂, ..., Pₖ]
- p: parameters (λ_vec, R_interp)
- t: time
"""
function protein_dynamics_mixture!(du, u, p, t)
    λ_vec, R_interp = p

    R_t = R_interp(t)

    for k in eachindex(λ_vec)
        du[k] = λ_vec[k] * R_t - λ_vec[k] * u[k]
    end

    return nothing
end

# ============================================================================
# ODE Solvers
# ============================================================================

"""
    solve_single_rate_ode(R_interp, tspan, λ)

Solve single-rate protein dynamics using DifferentialEquations.jl
"""
function solve_single_rate_ode(R_interp, tspan, λ; saveat=nothing)
    u0 = [1.0]  # Initial condition: P(0) = 1
    p = (λ, R_interp)

    prob = ODEProblem(protein_dynamics_single!, u0, tspan, p)

    if isnothing(saveat)
        sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    else
        sol = solve(prob, Tsit5(), saveat=saveat, reltol=1e-8, abstol=1e-8)
    end

    return sol
end

"""
    solve_mixture_ode(R_interp, tspan, λ_vec, π_vec)

Solve mixture model using DifferentialEquations.jl
"""
function solve_mixture_ode(R_interp, tspan, λ_vec, π_vec; saveat=nothing)
    u0 = π_vec  # Initial conditions: Pₖ(0) = πₖ
    p = (λ_vec, R_interp)

    prob = ODEProblem(protein_dynamics_mixture!, u0, tspan, p)

    if isnothing(saveat)
        sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    else
        sol = solve(prob, Tsit5(), saveat=saveat, reltol=1e-8, abstol=1e-8)
    end

    return sol
end

# ============================================================================
# Data Generation
# ============================================================================

function generate_rna_knockdown(t_days)
    k_rna = log(2) / 2.0
    return @. exp(-k_rna * t_days)
end

function generate_mixture_protein_data(t_data, λ_true, π_true, R_data; noise_σ=0.03)
    R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
    tspan = (0.0, maximum(t_data))

    sol = solve_mixture_ode(R_interp, tspan, λ_true, π_true, saveat=t_data)

    # Sum across subpopulations and add noise
    P_total = [sum(sol.u[i]) for i in eachindex(sol.u)]
    P_measured = P_total .+ noise_σ * randn(length(P_total))

    return P_measured
end

# ============================================================================
# Objective Functions for Fitting
# ============================================================================

"""
    objective_single_diffeq(params, t_data, R_data, P_data)

Objective function using DifferentialEquations.jl solver
"""
function objective_single_diffeq(params, t_data, R_data, P_data)
    λ = params[1]

    if λ <= 0
        return Inf
    end

    try
        R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
        tspan = (0.0, maximum(t_data))

        sol = solve_single_rate_ode(R_interp, tspan, λ, saveat=t_data)

        P_pred = [sol.u[i][1] for i in eachindex(sol.u)]

        return sum((P_data .- P_pred).^2)
    catch e
        return Inf
    end
end

"""
    objective_mixture_diffeq(params, t_data, R_data, P_data, K)

Objective function for mixture model using DifferentialEquations.jl
"""
function objective_mixture_diffeq(params, t_data, R_data, P_data, K)
    λ_vec = params[1:K]

    if K > 1
        π_vec = [params[K+1:end]..., 1 - sum(params[K+1:end])]
    else
        π_vec = [1.0]
    end

    # Parameter validation
    if any(λ_vec .<= 0) || any(π_vec .<= 0) || any(π_vec .>= 1) || sum(π_vec) > 1.001
        return Inf
    end

    try
        R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
        tspan = (0.0, maximum(t_data))

        sol = solve_mixture_ode(R_interp, tspan, λ_vec, π_vec, saveat=t_data)

        # Sum subpopulations
        P_pred = [sum(sol.u[i]) for i in eachindex(sol.u)]

        return sum((P_data .- P_pred).^2)
    catch e
        return Inf
    end
end

# ============================================================================
# Fitting Functions
# ============================================================================

function fit_single_rate_diffeq(t_data, R_data, P_data; λ_init=log(2)/5)
    params_init = [λ_init]

    result = optimize(
        p -> objective_single_diffeq(p, t_data, R_data, P_data),
        params_init,
        BFGS(),
        Optim.Options(show_trace=false)
    )

    λ_fit = result.minimizer[1]
    t_half = log(2) / λ_fit

    println("Single-rate fit (DiffEq):")
    println("  λ = $(round(λ_fit, digits=4))")
    println("  t₁/₂ = $(round(t_half, digits=2)) days")
    println("  RSS = $(round(result.minimum, digits=6))")

    return λ_fit, result
end

function fit_mixture_diffeq(t_data, R_data, P_data, K=2;
                           λ_init=nothing, π_init=nothing)
    if isnothing(λ_init)
        λ_init = [log(2)/2, log(2)/10]
    end
    if isnothing(π_init)
        π_init = fill(1/K, K-1)
    end

    params_init = [λ_init..., π_init...]

    lower = [zeros(K)..., zeros(K-1)...]
    upper = [fill(Inf, K)..., fill(0.99, K-1)...]

    result = optimize(
        p -> objective_mixture_diffeq(p, t_data, R_data, P_data, K),
        lower,
        upper,
        params_init,
        Fminbox(BFGS()),
        Optim.Options(show_trace=false, iterations=10000)
    )

    λ_vec = result.minimizer[1:K]
    π_vec = [result.minimizer[K+1:end]..., 1 - sum(result.minimizer[K+1:end])]
    t_half_vec = log(2) ./ λ_vec

    println("\n$(K)-component mixture fit (DiffEq):")
    for k in 1:K
        println("  Component $k: λ = $(round(λ_vec[k], digits=4)), " *
                "t₁/₂ = $(round(t_half_vec[k], digits=2)) days, " *
                "π = $(round(π_vec[k], digits=3))")
    end
    println("  RSS = $(round(result.minimum, digits=6))")

    return λ_vec, π_vec, result
end

# ============================================================================
# Visualization
# ============================================================================

function plot_diffeq_comparison(t_data, R_data, P_data, λ_single, λ_vec, π_vec)
    # Generate predictions
    R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
    tspan = (0.0, maximum(t_data))
    t_fine = range(tspan[1], tspan[2], length=1000)

    # Single rate
    sol_single = solve_single_rate_ode(R_interp, tspan, λ_single, saveat=t_fine)
    P_single = [sol_single.u[i][1] for i in eachindex(sol_single.u)]

    # Mixture
    sol_mixture = solve_mixture_ode(R_interp, tspan, λ_vec, π_vec, saveat=t_fine)
    P_sub = hcat([sol_mixture.u[i] for i in eachindex(sol_mixture.u)]...)'
    P_mixture = sum(P_sub, dims=2)[:]

    # Create figure
    fig = Figure(size=(1000, 600), fontsize=14)

    # Panel A: Linear scale comparison
    ax1 = Axis(fig[1, 1],
        xlabel="Time (days)",
        ylabel="Protein Fraction Remaining",
        title="DifferentialEquations.jl Implementation")

    scatter!(ax1, t_data, P_data, label="Data", color=:black, markersize=10)
    lines!(ax1, t_fine, P_single, label="Single Rate",
           color=:red, linestyle=:dash, linewidth=2)
    lines!(ax1, t_fine, P_mixture, label="Mixture (K=2)",
           color=:blue, linewidth=3)

    axislegend(ax1, position=:rt)

    # Panel B: Subpopulation dynamics
    ax2 = Axis(fig[1, 2],
        xlabel="Time (days)",
        ylabel="Subpopulation Fraction",
        title="Subpopulation Trajectories")

    for k in 1:size(P_sub, 2)
        t_half_k = log(2) / λ_vec[k]
        lines!(ax2, t_fine, P_sub[:, k],
               label="P$(k) (t₁/₂=$(round(t_half_k, digits=1))d, π=$(round(π_vec[k], digits=2)))",
               linewidth=2)
    end

    axislegend(ax2, position=:rt)

    return fig
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("="^70)
    println("Mixture Decay Model with DifferentialEquations.jl")
    println("="^70)

    # Generate synthetic data
    println("\n1. Generating synthetic data...")
    t_data = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0]
    R_data = generate_rna_knockdown(t_data)

    λ_true = [log(2)/2, log(2)/8]
    π_true = [0.7, 0.3]

    P_data = generate_mixture_protein_data(t_data, λ_true, π_true, R_data, noise_σ=0.02)

    println("True parameters:")
    println("  Component 1: t₁/₂ = 2.0 days, π = 0.7")
    println("  Component 2: t₁/₂ = 8.0 days, π = 0.3")

    # Fit models
    println("\n2. Fitting single-rate model...")
    λ_single, res_single = fit_single_rate_diffeq(t_data, R_data, P_data)

    println("\n3. Fitting mixture model...")
    λ_vec, π_vec, res_mixture = fit_mixture_diffeq(t_data, R_data, P_data, 2)

    # Compare
    println("\n4. Model comparison:")
    println("  RSS improvement: $(round((1 - res_mixture.minimum/res_single.minimum)*100, digits=1))%")

    # Visualize
    println("\n5. Creating plot...")
    fig = plot_diffeq_comparison(t_data, R_data, P_data, λ_single, λ_vec, π_vec)

    save("mixture_decay_diffeq.png", fig)
    println("Saved: mixture_decay_diffeq.png")

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
