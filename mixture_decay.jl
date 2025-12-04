"""
Mixture Decay Model - Balanced Equilibrium (Fixed)

Implements:
    dPₖ/dt = (ψₖ * πₖ) * R(t)  -  ψₖ * Pₖ(t)
"""

using Optim, CairoMakie, LaTeXStrings, LinearAlgebra, Interpolations, Statistics
using DataFrames, CSV, Colors, HypothesisTests, Printf, Distributions

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
    k_rna = log(2) / 0.1
    return @. 0.5 + (0.5 * exp(-k_rna * t_days))
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
    t_data = collect(0.0:10.0:50.0)
    R_data = generate_rna_knockdown(t_data)
    
    ψ_true = [1.0, 0.05]
    π_true = [0.5, 0.5]
    
    # Generate data using the NEW balanced solver
    P_data = generate_mixture_protein_data(t_data, ψ_true, π_true, R_data, noise_σ=0.01)

    # 2. Fit
    λ_single = fit_single_rate_model(t_data, R_data, P_data)
    λ_mix, ψ_fit, π_fit = fit_mixture_model(t_data, R_data, P_data, 2)

    # 3. Plot
    t_fine = range(0, 50, length=200)
    R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())
    
    P_single = solve_protein_dynamics(R_interp, t_fine, λ_single)
    P_mix, P_sub = solve_mixture_dynamics(R_interp, t_fine, ψ_fit, π_fit)

    fig = Figure(size=(500, 250),fontsize=10)
    ax1 = Axis(fig[1,1], title="Simulation", xlabel="Time Steps", ylabel="Measured Quantity")
    scatter!(ax1, t_data, P_data, color=:black, label="Data")
    lines!(ax1, t_fine, P_single, color=:red, linestyle=:dash, label=L"dP = λR_{t-1}dt-λP_{t-1}dt")
    lines!(ax1, t_fine, P_mix, color=:blue, linewidth=2, label=L"dP = \sum_k{\left(\sum_k{w_k λ_k}\right)R_{t-1}dt - λ_k P_{k,t-1}dt}")
    Legend(fig[1,2],ax1)

    # ax2 = Axis(fig[1,2], title="Subpopulations", xlabel="Time")
    # lines!(ax2, t_fine, P_sub[:,1], color=:purple, label="Fast")
    # lines!(ax2, t_fine, P_sub[:,2], color=:orange, label="Slow")
    # axislegend(ax2,position=:rt)

    save("balanced_decay.png", fig)
    println("Saved plot: balanced_decay.png")
end

main()

# ============================================================================
# FFI Analysis - Helper Functions
# ============================================================================

# Time step for simulation (same as R code)
const dt = 0.01

"""
Free lysine availability function from Fornasiero 2018 Figure S3H
Models the gradual incorporation of labeled lysine into the free amino acid pool
"""
function free_lysine(t)
    return max(0.0, 1.0 - 0.503 * exp(-t * log(2) * 0.799) - 0.503 * exp(-t * log(2) / 39.423))
end

"""
Always 100% labeled lysine available (alternative model)
"""
function always_all(t)
    return 1.0
end

"""
Simulate proportion labeled over time for SINGLE RATE model
Based on: dP/dt = λ * (R_avail(t) - P(t))
where R_avail is the availability of labeled lysine
"""
function proportion_labeled_single(thalf, t_points, avails_func=free_lysine)
    lambda = log(2) / thalf
    t = collect(0:dt:maximum(t_points))
    proportion_heavy = zeros(length(t))

    for i in 2:length(t)
        protein_turned_over = lambda * dt
        original = (1 - protein_turned_over) * proportion_heavy[i-1]
        nascent = protein_turned_over * avails_func(t[i])
        proportion_heavy[i] = original + nascent
    end

    # Interpolate to requested time points
    interp = LinearInterpolation(t, proportion_heavy, extrapolation_bc=Line())
    return [interp(tp) for tp in t_points]
end

"""
Simulate proportion labeled over time for MIXTURE model
Two subpopulations with different decay rates
"""
function proportion_labeled_mixture(thalves, proportions, t_points, avails_func=free_lysine)
    t = collect(0:dt:maximum(t_points))
    K = length(thalves)

    # Track each subpopulation separately
    prop_heavy_sub = zeros(length(t), K)

    for k in 1:K
        lambda_k = log(2) / thalves[k]
        for i in 2:length(t)
            protein_turned_over = lambda_k * dt
            original = (1 - protein_turned_over) * prop_heavy_sub[i-1, k]
            nascent = protein_turned_over * avails_func(t[i])
            prop_heavy_sub[i, k] = original + nascent
        end
    end

    # Weight by proportions
    total_prop_heavy = prop_heavy_sub * proportions

    # Interpolate to requested time points
    interp = LinearInterpolation(t, total_prop_heavy, extrapolation_bc=Line())
    return [interp(tp) for tp in t_points]
end

"""
Calculate 95% confidence interval lower bound using t-distribution
"""
function lower(x)
    n = length(x)
    if n < 2
        return mean(x)
    end
    m = mean(x)
    s = std(x)
    t_val = quantile(TDist(n-1), 0.975)
    return m - t_val * s / sqrt(n)
end

"""
Calculate 95% confidence interval upper bound using t-distribution
"""
function upper(x)
    n = length(x)
    if n < 2
        return mean(x)
    end
    m = mean(x)
    s = std(x)
    t_val = quantile(TDist(n-1), 0.975)
    return m + t_val * s / sqrt(n)
end

"""
Fit single-rate isotopic half-life from labeling time course data using free lysine model
Uses optimization to find half-life that best matches observed labeling
"""
function fit_isotopic_thalf_single(chow_days, prop_labeled, avails_func=free_lysine)
    if length(chow_days) < 2
        return NaN
    end

    # Objective: minimize sum of squared residuals
    function objective(thalf)
        if thalf <= 0.01 || thalf > 100
            return Inf
        end
        pred = proportion_labeled_single(thalf, chow_days, avails_func)
        return sum((prop_labeled .- pred).^2)
    end

    # Optimize
    try
        result = optimize(objective, 0.1, 100.0)
        return result.minimizer
    catch
        return NaN
    end
end

"""
Fit mixture model isotopic half-lives from labeling time course data
Returns two half-lives and their proportions
"""
function fit_isotopic_thalf_mixture(chow_days, prop_labeled, avails_func=free_lysine)
    if length(chow_days) < 4  # Need more data for mixture model
        return (thalves=[NaN, NaN], proportions=[0.5, 0.5], effective_thalf=NaN)
    end

    # Objective: minimize sum of squared residuals
    function objective(params)
        thalf1, thalf2, pi1 = params

        # Constraints
        if thalf1 <= 0.01 || thalf1 > 100 ||
           thalf2 <= 0.01 || thalf2 > 100 ||
           pi1 <= 0 || pi1 >= 1
            return Inf
        end

        # Ensure thalf1 < thalf2 (fast first)
        if thalf1 > thalf2
            thalf1, thalf2 = thalf2, thalf1
        end

        pi2 = 1 - pi1
        pred = proportion_labeled_mixture([thalf1, thalf2], [pi1, pi2], chow_days, avails_func)
        return sum((prop_labeled .- pred).^2)
    end

    # Try optimization with multiple initial conditions
    best_result = nothing
    best_obj = Inf

    for thalf1_init in [1.0, 3.0, 5.0]
        for thalf2_init in [8.0, 15.0, 25.0]
            for pi1_init in [0.3, 0.5, 0.7]
                try
                    result = optimize(
                        objective,
                        [thalf1_init, thalf2_init, pi1_init],
                        NelderMead()
                    )
                    if result.minimum < best_obj
                        best_obj = result.minimum
                        best_result = result
                    end
                catch
                    continue
                end
            end
        end
    end

    if best_result === nothing
        return (thalves=[NaN, NaN], proportions=[0.5, 0.5], effective_thalf=NaN)
    end

    thalf1, thalf2, pi1 = best_result.minimizer
    if thalf1 > thalf2
        thalf1, thalf2 = thalf2, thalf1
    end
    pi2 = 1 - pi1

    # Calculate effective half-life (weighted average)
    effective_thalf = pi1 * thalf1 + pi2 * thalf2

    return (thalves=[thalf1, thalf2], proportions=[pi1, pi2], effective_thalf=effective_thalf)
end

"""
Calculate labeling differences between FFI and control groups using BOTH models
Returns results for both single-rate and mixture models
"""
function label_difference(prop_labeled, chow_days, genotype, avails_func=free_lysine)
    # Identify FFI vs control groups
    is_ffi = genotype .== "129(TT-3F4-FFI)HOZ"

    # SINGLE RATE MODEL
    thalf_single_ffi = fit_isotopic_thalf_single(chow_days[is_ffi], prop_labeled[is_ffi], avails_func)
    thalf_single_ctrl = fit_isotopic_thalf_single(chow_days[.!is_ffi], prop_labeled[.!is_ffi], avails_func)
    thalf_single_ratio = thalf_single_ffi / thalf_single_ctrl

    # MIXTURE MODEL
    mixture_ffi = fit_isotopic_thalf_mixture(chow_days[is_ffi], prop_labeled[is_ffi], avails_func)
    mixture_ctrl = fit_isotopic_thalf_mixture(chow_days[.!is_ffi], prop_labeled[.!is_ffi], avails_func)
    thalf_mixture_ratio = mixture_ffi.effective_thalf / mixture_ctrl.effective_thalf

    # Simple t-test on labeling rates
    br_ffi_pval = 1.0  # Default value
    try
        test = UnequalVarianceTTest(prop_labeled[is_ffi], prop_labeled[.!is_ffi])
        br_ffi_pval = pvalue(test)
    catch
        br_ffi_pval = 1.0
    end

    # Interaction p-value (simplified - would need proper linear model)
    br_interaction_pval = 1.0  # Placeholder

    return (
        # Single rate results
        thalf_single_ratio = thalf_single_ratio,
        thalf_single_ffi = thalf_single_ffi,
        thalf_single_ctrl = thalf_single_ctrl,
        # Mixture results
        thalf_mixture_ratio = thalf_mixture_ratio,
        thalf_mixture_ffi = mixture_ffi.effective_thalf,
        thalf_mixture_ctrl = mixture_ctrl.effective_thalf,
        mixture_ffi_thalves = mixture_ffi.thalves,
        mixture_ffi_proportions = mixture_ffi.proportions,
        mixture_ctrl_thalves = mixture_ctrl.thalves,
        mixture_ctrl_proportions = mixture_ctrl.proportions,
        # P-values
        br_ffi_pval = br_ffi_pval,
        br_interaction_pval = br_interaction_pval
    )
end


"""
Generate comparison figure showing label accumulation curves for both models
Side-by-side comparison: Single (left) vs Mixture (right)
"""
function generate_figure_model_comparison()
    println("Creating model comparison figure...")

    # Load data
    ffi_all = CSV.read("data/iqp/ffi_turnover.tsv", DataFrame)
    ffi_all[!, :total] = ffi_all.light .+ ffi_all.heavy
    ffi_all[!, :prop_labeled] = ffi_all.heavy ./ ffi_all.total

    # Genotype legend
    leg = DataFrame(
        genotype = ["129(TT-3F4-FFI)HOZ", "129(TT-3F4WT)", "B6/N"],
        disp = ["ki-3F4-FFI", "ki-3F4-WT", "C57BL/6N WT"],
        color = [colorant"#D95F02", colorant"#22127A", colorant"#77127A"],
        xgeno = [1, 2, 3],
        ttest_grouping = ["test", "control", "control"]
    )

    # Peptides to plot
    use_peptides = ["GENFTETDVK", "VVEQMCVTQYQK"]
    ages = ["young", "aged"]

    # Create figure with two GridLayouts (wider to accommodate two 2x2 grids)
    fig = Figure(size=(1600, 800), fontsize=12)

    # Add main title labels at the top
    Label(fig[1, 1], "Single-Rate Model", fontsize=16, font=:bold, halign=:center)
    Label(fig[1, 2], "Mixture Model", fontsize=16, font=:bold, halign=:center)

    # Create GridLayouts for the plots
    gl_single = fig[2, 1] = GridLayout()
    gl_mixture = fig[2, 2] = GridLayout()

    # Panel labels
    panel_letters = ['C', 'D', 'G', 'H']
    avails_func = free_lysine

    # Loop through ages and peptides
    for (age_idx, this_age) in enumerate(ages)
        for (pep_idx, this_peptide) in enumerate(use_peptides)
            panel_idx = (age_idx - 1) * 2 + pep_idx

            # Create axes for both models with square aspect
            ax_single = Axis(gl_single[age_idx, pep_idx],
                xlabel=age_idx == 2 ? "day" : "",
                ylabel=pep_idx == 1 ? "proportion labeled" : "",
                title="$(this_age) $(this_peptide[1:4])",
                width = 200,
                height = 200
                )

            ax_mixture = Axis(gl_mixture[age_idx, pep_idx],
                xlabel=age_idx == 2 ? "day" : "",
                ylabel=pep_idx == 1 ? "proportion labeled" : "",
                title="$(this_age) $(this_peptide[1:4])",
                width = 200,
                height = 200
                )

            xlims!(ax_single, -0.5, 8.5)
            xlims!(ax_mixture, -0.5, 8.5)
            ylims!(ax_single, 0, 0.5)
            ylims!(ax_mixture, 0, 0.5)

            # Fine time grid for model predictions
            t_fine = collect(0:0.1:8)

            # Plot each genotype
            for geno_row in eachrow(leg)
                # Filter data
                curve_data = filter(row -> row.genotype == geno_row.genotype &&
                                          row.peptide == this_peptide &&
                                          row.age == this_age, ffi_all)

                if nrow(curve_data) == 0
                    continue
                end

                # Calculate summary statistics by day
                summary_df = combine(groupby(curve_data, :chow_days)) do df
                    DataFrame(
                        mean_prop = mean(df.prop_labeled),
                        l95 = lower(df.prop_labeled),
                        u95 = upper(df.prop_labeled)
                    )
                end
                sort!(summary_df, :chow_days)

                # SINGLE MODEL FIT
                if nrow(curve_data) > 0
                    thalf_single = fit_isotopic_thalf_single(
                        curve_data.chow_days,
                        curve_data.prop_labeled,
                        avails_func
                    )

                    if !isnan(thalf_single)
                        pred_single = proportion_labeled_single(thalf_single, t_fine, avails_func)
                        lines!(ax_single, t_fine, pred_single,
                              color=geno_row.color, linewidth=2, linestyle=:dash)
                    end
                end

                # MIXTURE MODEL FIT
                if nrow(curve_data) > 0
                    mixture_fit = fit_isotopic_thalf_mixture(
                        curve_data.chow_days,
                        curve_data.prop_labeled,
                        avails_func
                    )

                    if !isnan(mixture_fit.effective_thalf)
                        pred_mixture = proportion_labeled_mixture(
                            mixture_fit.thalves,
                            mixture_fit.proportions,
                            t_fine,
                            avails_func
                        )
                        lines!(ax_mixture, t_fine, pred_mixture,
                              color=geno_row.color, linewidth=2, linestyle=:dash)
                    end
                end

                # Plot observed data on both panels
                for ax in [ax_single, ax_mixture]
                    band!(ax, summary_df.chow_days, summary_df.l95, summary_df.u95,
                         color=RGBA(geno_row.color, 0.3))
                    scatter!(ax, curve_data.chow_days, curve_data.prop_labeled,
                            color=geno_row.color, markersize=6,
                            marker=:circle, strokewidth=1, strokecolor=geno_row.color)
                end
            end

            # # Add panel labels
            # text!(ax_single, -0.3, 0.47, text=string(panel_letters[panel_idx]),
            #      fontsize=14, font=:bold)
            # text!(ax_mixture, -0.3, 0.47, text=string(panel_letters[panel_idx]),
            #      fontsize=14, font=:bold)
        end
    end

    # Add single legend for the entire figure at the bottom
    legend_elements = [
        [MarkerElement(marker=:circle, color=c, strokecolor=c, markersize=10) for c in leg.color]...,
        LineElement(color=:black, linewidth=2, linestyle=:dash)
    ]
    legend_labels = [leg.disp..., "model fit"]

    Legend(fig[3, :],
          legend_elements, legend_labels,
          orientation=:horizontal,
          framevisible=false,
          halign=:center)

    filename = "display_items/figure-model-comparison.png"
    resize_to_layout!(fig)

    save(filename, fig, px_per_unit=2)
    println("done.\nModel comparison figure saved to $(filename)")
#    rowsize!(fig.layout,1,Aspect(1,1.0))
    return fig
end

# Run Figure 4 generation with both models
#generate_all_figure4()

# Generate model comparison figure
generate_figure_model_comparison()