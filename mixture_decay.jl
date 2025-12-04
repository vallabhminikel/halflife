"""
Mixture Decay Model - Balanced Equilibrium (Fixed)

Implements:
    dPₖ/dt = (ψₖ * πₖ) * R(t)  -  ψₖ * Pₖ(t)
"""

using Optim, CairoMakie, LaTeXStrings, LinearAlgebra, Interpolations, Statistics
using DataFrames, CSV, Colors, HypothesisTests, Printf

# ============================================================================
# Core ODE Solver
# ============================================================================

function solve_protein_dynamics(R_interp, t, λ)
    P = zeros(length(t))
    P[1] = 1.0
    for i in 2:length(t)
        dt = t[i] - t[i-1]
        R_prev = R_interp(t[i-1])
        dP = λ * dt * R_prev - λ * dt * P[i-1]
        P[i] = P[i-1] + dP
    end
    return P
end

"""
    solve_mixture_dynamics(R_interp, t, ψ_vec, π_vec)

Solve the balanced mixture model.
Production of k is scaled by π[k]*ψ[k] to ensure initial equilibrium.
"""
function solve_mixture_dynamics(R_interp, t, ψ_vec, π_vec)
    K = length(ψ_vec)
    n = length(t)
    P_sub = zeros(n, K)

    # Initial conditions
    P_sub[1, :] = π_vec


    prod_rates = ψ_vec .* π_vec

    for k in 1:K
        for i in 2:n
            dt = t[i] - t[i-1]
            R_prev = R_interp(t[i-1])
            
            term_prod = prod_rates[k] * R_prev
            term_deg  = ψ_vec[k] * P_sub[i-1, k]
            
            dP = (term_prod - term_deg) * dt
            
            P_sub[i, k] = P_sub[i-1, k] + dP
        end
    end

    P_total = sum(P_sub, dims=2)[:]
    return P_total, P_sub
end

# ============================================================================
# Data Generation & Objective
# ============================================================================

function generate_rna_knockdown(t_days)
    k_rna = log(2) / 0.5
    return @. exp(-k_rna * t_days)
end

function generate_mixture_protein_data(t_days, ψ_true, π_true, R_values; noise_σ=0.03)
    t_fine = range(0, maximum(t_days), length=1000)
    R_interp = LinearInterpolation(t_days, R_values, extrapolation_bc=Line())
    P_true, _ = solve_mixture_dynamics(R_interp, t_fine, ψ_true, π_true)
    P_interp = LinearInterpolation(t_fine, P_true, extrapolation_bc=Line())
    return [P_interp(t) for t in t_days] .+ noise_σ * randn(length(t_days))
end

# Objective for Single Rate Model
function objective_single_rate(params, t_data, R_data, P_data)
    λ = params[1]
    if λ <= 0 return Inf end

    t_fine = range(0, maximum(t_data), length=1000)
    R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
    
    P_pred = solve_protein_dynamics(R_interp, t_fine, λ)
    P_interp = LinearInterpolation(t_fine, P_pred, extrapolation_bc=Line())
    
    return sum((P_data .- [P_interp(t) for t in t_data]).^2)
end

# Objective for Mixture Model
function objective_mixture(params, t_data, R_data, P_data, K)
    ψ_vec = params[1:K]
    if K > 1
        π_partial = params[K+1:end]
        π_vec = [π_partial..., 1 - sum(π_partial)]
    else
        π_vec = [1.0]
    end

    if any(ψ_vec .<= 0) || any(π_vec .< 0) || sum(π_vec) > 1.001
        return Inf
    end

    t_fine = range(0, maximum(t_data), length=1000)
    R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
    P_pred, _ = solve_mixture_dynamics(R_interp, t_fine, ψ_vec, π_vec)
    P_interp = LinearInterpolation(t_fine, P_pred, extrapolation_bc=Line())
    
    return sum((P_data .- [P_interp(t) for t in t_data]).^2)
end

# ============================================================================
# Fitting
# ============================================================================

function fit_single_rate_model(t_data, R_data, P_data)
    # FIX: Cleaned up optimize call
    result = optimize(
        p -> objective_single_rate(p, t_data, R_data, P_data),
        [log(2)/5], 
        BFGS(), 
        Optim.Options(show_trace=false)
    )
    return result.minimizer[1]
end

function fit_mixture_model(t_data, R_data, P_data, K=2)
    ψ_init = [log(2)/2, log(2)/10] 
    π_init = fill(1/K, K-1)
    params_init = [ψ_init..., π_init...]

    lower = [fill(1e-5, K)..., fill(0.0, K-1)...]
    upper = [fill(Inf, K)..., fill(0.99, K-1)...]

    result = optimize(
        p -> objective_mixture(p, t_data, R_data, P_data, K),
        lower, upper, params_init, Fminbox(BFGS()), Optim.Options(show_trace=false)
    )

    ψ_vec = result.minimizer[1:K]
    π_partial = result.minimizer[K+1:end]
    π_vec = [π_partial..., 1 - sum(π_partial)]
    
    # Calculate Total Effective Translation Rate
    λ_effective = sum(ψ_vec .* π_vec)

    println("\n$(K)-component Balanced Fit:")
    println("  Effective Total λ = $(round(λ_effective, digits=4))")
    for k in 1:K
        println("  Component $k: ψ=$(round(ψ_vec[k], digits=4)), π=$(round(π_vec[k], digits=3))")
    end

    return λ_effective, ψ_vec, π_vec
end

# ============================================================================
# Plotting & Main
# ============================================================================

function main()
    println("Running Balanced Mixture Model (No Overshoot)...")

    # 1. Generate Data
    t_data = collect(0.0:2.0:8.0)
    R_data = generate_rna_knockdown(t_data)
    
    ψ_true = [0.5, 0.005]
    π_true = [0.5, 0.5]
    
    # Generate data using the NEW balanced solver
    P_data = generate_mixture_protein_data(t_data, ψ_true, π_true, R_data, noise_σ=0.01)

    # 2. Fit
    λ_single = fit_single_rate_model(t_data, R_data, P_data)
    λ_mix, ψ_fit, π_fit = fit_mixture_model(t_data, R_data, P_data, 2)

    # 3. Plot
    t_fine = range(0, 8, length=200)
    R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
    
    P_single = solve_protein_dynamics(R_interp, t_fine, λ_single)
    P_mix, P_sub = solve_mixture_dynamics(R_interp, t_fine, ψ_fit, π_fit)

    fig = Figure(size=(500, 250),fontsize=10)
    ax1 = Axis(fig[1,1], title="Fits", xlabel="Time", ylabel="Unlabeled Protein")
    scatter!(ax1, t_data, P_data, color=:black, label="Data")
    lines!(ax1, t_fine, P_single, color=:red, linestyle=:dash, label="Single")
    lines!(ax1, t_fine, P_mix, color=:blue, linewidth=2, label="Mixture")
    axislegend(ax1)

    ax2 = Axis(fig[1,2], title="Subpopulations", xlabel="Time")
    lines!(ax2, t_fine, P_sub[:,1], color=:purple, label="Fast")
    lines!(ax2, t_fine, P_sub[:,2], color=:orange, label="Slow")
    axislegend(ax2,position=:lb)

    save("balanced_decay.png", fig)
    println("Saved plot: balanced_decay.png")
end

main()