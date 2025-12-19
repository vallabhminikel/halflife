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
            dt = t[i] - t[i-1]a
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

    function objective(params)
        thalf1, thalf2, pi1 = params

        if thalf1 <= 0.01 || thalf1 > 100 ||
           thalf2 <= 0.01 || thalf2 > 100 ||
           pi1 <= 0 || pi1 >= 1
            return Inf
        end

        # Keep parameters as-is, don't swap in objective
        # The optimizer will find the right assignment
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
        pi1 = 1 - pi1  # Also swap the proportions!
    end
    pi2 = 1 - pi1

    # Calculate effective half-life: time when mixture population decays to 50%
    # π₁ × (1/2)^(t/t½₁) + π₂ × (1/2)^(t/t½₂) = 0.5
    function mixture_remaining(t)
        return pi1 * (0.5^(t/thalf1)) + pi2 * (0.5^(t/thalf2))
    end
    
    # Find t where mixture_remaining(t) = 0.5
    function objective_halflife(t)
        return abs(mixture_remaining(t) - 0.5)
    end
    
    effective_thalf = pi1 * thalf1 + pi2 * thalf2  
    try
        result = optimize(objective_halflife, 0.0, max(thalf1, thalf2) * 2.0)
        effective_thalf = result.minimizer
    catch
        @error "Error"
    end

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
Generate text report of mixture model components and half-lives
"""
function generate_mixture_report()
    println("=" ^ 80)
    println("MIXTURE MODEL ANALYSIS REPORT")
    println("=" ^ 80)
    println()

    # Load data
    ffi_all = CSV.read("data/iqp/ffi_turnover.tsv", DataFrame)
    ffi_all[!, :total] = ffi_all.light .+ ffi_all.heavy
    ffi_all[!, :prop_labeled] = ffi_all.heavy ./ ffi_all.total

    # Genotype info
    genotypes = ["129(TT-3F4-FFI)HOZ", "129(TT-3F4WT)", "B6/N"]
    genotype_labels = ["ki-3F4-FFI", "ki-3F4-WT", "C57BL/6N WT"]

    # Peptides and ages
    peptides = ["GENF", "VVEQ"]
    ages = ["young", "aged"]

    avails_func = free_lysine

    # Generate report for each combination
    for this_age in ages
        println("\n" * "=" ^ 80)
        println("AGE: $(uppercase(this_age))")
        println("=" ^ 80)

        for this_peptide in peptides
            println("\n" * "-" ^ 80)
            println("Peptide: $(this_peptide)")
            println("-" ^ 80)

            for (geno_idx, this_genotype) in enumerate(genotypes)
                # Filter data
                data = filter(row -> row.genotype == this_genotype &&
                                    row.peptide == this_peptide &&
                                    row.age == this_age, ffi_all)

                if nrow(data) == 0
                    println("\n  $(genotype_labels[geno_idx]): No data")
                    continue
                end

                # Fit single model
                thalf_single = fit_isotopic_thalf_single(
                    data.chow_days,
                    data.prop_labeled,
                    avails_func
                )

                # Fit mixture model
                mixture_fit = fit_isotopic_thalf_mixture(
                    data.chow_days,
                    data.prop_labeled,
                    avails_func
                )

                # Print results
                println("\n  $(genotype_labels[geno_idx]):")
                println("    n = $(nrow(data)) observations")

                if !isnan(thalf_single)
                    @printf("    Single-rate model: t½ = %.2f days\n", thalf_single)
                else
                    println("    Single-rate model: Failed to converge")
                end

                if !isnan(mixture_fit.effective_thalf)
                    println("    Mixture model:")
                    @printf("      Component 1 (fast):  t½ = %.2f days, π = %.1f%%\n",
                           mixture_fit.thalves[1], mixture_fit.proportions[1] * 100)
                    @printf("      Component 2 (slow):  t½ = %.2f days, π = %.1f%%\n",
                           mixture_fit.thalves[2], mixture_fit.proportions[2] * 100)
                    @printf("      Effective t½ = %.2f days\n", mixture_fit.effective_thalf)

                    # Calculate improvement over single model
                    if !isnan(thalf_single)
                        ratio = mixture_fit.effective_thalf / thalf_single
                        @printf("      Ratio (mixture/single) = %.2f\n", ratio)
                    end
                else
                    println("    Mixture model: Failed to converge")
                end
            end
        end
    end

    println("\n" * "=" ^ 80)
    println("END OF REPORT")
    println("=" ^ 80)
end

"""
Prepare data for Figure 4: Load data, compute genotypic differences with mixture model
Returns: (ffi_all, genotypic_diffs, leg)
"""
function prepare_figure_4_data()
    println("Loading and processing data for Figure 4...")

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

    # Compute genotypic differences
    println("Computing genotypic differences with mixture model...")
    genotypic_diffs = combine(groupby(
        innerjoin(ffi_all, leg, on=:genotype),
        [:protein, :peptide, :age]
    )) do df
        # Filter to peptides detected in both ages
        n_ffi = sum(df.genotype .== "129(TT-3F4-FFI)HOZ")
        n_ctrl = sum(df.genotype .!= "129(TT-3F4-FFI)HOZ")
        
        if n_ffi == 0 || n_ctrl == 0
            return DataFrame()
        end

        # Abundance test
        total_ffi = df.total[df.genotype .== "129(TT-3F4-FFI)HOZ"]
        total_ctrl = df.total[df.genotype .!= "129(TT-3F4-FFI)HOZ"]
        
        pval_total = 1.0
        ratio_total = 1.0
        try
            test = UnequalVarianceTTest(total_ffi, total_ctrl)
            pval_total = pvalue(test)
            ratio_total = mean(total_ffi) / mean(total_ctrl)
        catch
        end

        # Labeling/half-life test using mixture model
        label_obj = label_difference(df.prop_labeled, df.chow_days, df.genotype, free_lysine)
        
        DataFrame(
            n = nrow(df),
            pval_total = pval_total,
            ratio_total = ratio_total,
            ratio_thalf = label_obj.thalf_mixture_ratio,  # Use mixture model ratio
            pval_label_ffi = label_obj.br_ffi_pval,
            pval_label_interaction = label_obj.br_interaction_pval,
            log2_ratio_total = log2(ratio_total),
            log2_ratio_thalf = log2(label_obj.thalf_mixture_ratio),
            pepnick = first(df.peptide, 4)
        )
    end

    # Filter to peptides in both ages
    peptide_counts = combine(groupby(genotypic_diffs, [:protein, :peptide])) do df
        DataFrame(n_ages = length(unique(df.age)))
    end
    genotypic_diffs = innerjoin(genotypic_diffs, peptide_counts, on=[:protein, :peptide])
    filter!(row -> row.n_ages == 2, genotypic_diffs)
    select!(genotypic_diffs, Not(:n_ages))

    # Assign colors
    genotypic_diffs[!, :total_color] = map(genotypic_diffs.log2_ratio_total) do x
        x < 0 ? colorant"#FF0000" : x > 0 ? colorant"#0000FF" : colorant"#000000"
    end
    genotypic_diffs[!, :label_ffi_color] = map(genotypic_diffs.log2_ratio_thalf) do x
        x < 0 ? colorant"#FF0000" : x > 0 ? colorant"#0000FF" : colorant"#000000"
    end

    println("Data preparation complete.")
    return ffi_all, genotypic_diffs, leg
end

"""
Generate supplementary table with mixture model results
Includes: peptide, genotype, age, mixture components, runs test on residuals
"""
function generate_mixture_model_table()
    println("Generating mixture model supplementary table...")
    
    # Load data
    ffi_all = CSV.read("data/iqp/ffi_turnover.tsv", DataFrame)
    ffi_all[!, :total] = ffi_all.light .+ ffi_all.heavy
    ffi_all[!, :prop_labeled] = ffi_all.heavy ./ ffi_all.total
    
    # Genotypes and peptides
    genotypes = ["129(TT-3F4-FFI)HOZ", "129(TT-3F4WT)", "B6/N"]
    genotype_labels = ["ki-3F4-FFI", "ki-3F4-WT", "C57BL/6N WT"]
    peptides = ["GENFTETDVK", "VVEQMCVTQYQK"]
    ages = ["young", "aged"]
    
    results = DataFrame()
    
    for this_age in ages
        for this_peptide in peptides
            for (geno_idx, this_genotype) in enumerate(genotypes)
                # Filter data
                data = filter(row -> row.genotype == this_genotype &&
                                    row.peptide == this_peptide &&
                                    row.age == this_age, ffi_all)
                
                if nrow(data) == 0
                    continue
                end
                
                # Fit mixture model
                mixture_fit = fit_isotopic_thalf_mixture(
                    data.chow_days,
                    data.prop_labeled,
                    free_lysine
                )
                
                # Compute residuals for runs test using SINGLE exponential fit
                thalf_single = fit_isotopic_thalf_single(
                    data.chow_days,
                    data.prop_labeled,
                    free_lysine
                )
                
                if !isnan(thalf_single)
                    # Predict using single exponential model
                    pred_single = proportion_labeled_single(
                        thalf_single,
                        data.chow_days,
                        free_lysine
                    )
                    residuals = data.prop_labeled .- pred_single
                    log_residuals = log.(abs.(residuals) .+ 1e-10)  # Add small constant to avoid log(0)
                    runs_result = runs_test(log_residuals)
                else
                    runs_result = (statistic=NaN, p_value=NaN, n_runs=0, n_pos=0, n_neg=0)
                end
                
                # Create row for table
                row = DataFrame(
                    protein = data.protein[1],
                    peptide = this_peptide,
                    age = this_age,
                    genotype = this_genotype,
                    genotype_label = genotype_labels[geno_idx],
                    n_obs = nrow(data),
                    # Mixture component 1 (fast)
                    thalf_fast = !isnan(mixture_fit.effective_thalf) ? round(mixture_fit.thalves[1], digits=4) : NaN,
                    proportion_fast = !isnan(mixture_fit.effective_thalf) ? round(mixture_fit.proportions[1], digits=4) : NaN,
                    # Mixture component 2 (slow)
                    thalf_slow = !isnan(mixture_fit.effective_thalf) ? round(mixture_fit.thalves[2], digits=4) : NaN,
                    proportion_slow = !isnan(mixture_fit.effective_thalf) ? round(mixture_fit.proportions[2], digits=4) : NaN,
                    # Effective (weighted average)
                    thalf_effective = isnan(mixture_fit.effective_thalf) ? NaN : round(mixture_fit.effective_thalf, digits=4),
                    # Runs test on log residuals
                    runs_test_statistic = isnan(runs_result.statistic) ? NaN : round(runs_result.statistic, digits=4),
                    runs_test_pvalue = isnan(runs_result.p_value) ? NaN : round(runs_result.p_value, digits=4),
                    runs_test_n_runs = runs_result.n_runs,
                    runs_test_n_pos = runs_result.n_pos,
                    runs_test_n_neg = runs_result.n_neg
                )
                
                results = vcat(results, row)
            end
        end
    end
    
    # Save table
    output_file = "display_items/table-mixture-model-results.tsv"
    CSV.write(output_file, results, delim='\t')
    println("Mixture model table saved to $(output_file)")
    
    return results
end

ffi_all, genotypic_diffs, leg = prepare_figure_4_data()

fig = Figure(size=(800,400))
set_theme!(Axis = (xgridvisible = false, ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        rightticksvisible = false,
        topticksvisible = false,))

use_peptides = ["GENF", "VVEQ"]
full_peptides = ["GENFTETDVK", "VVEQMCVTQYQK"]
volcano_axes = []
for (ind,letter) ∈ zip([(1,1),(1,2),(3,1),(3,2)],["A","B","E","F"])
    ax = Axis(fig[ind[1],ind[2]],
    width = 200,
    height = 200)

    Label(fig[ind[1],ind[2],TopLeft()],letter,
        font = :bold,
        fontsize=18,
        padding = (0,0,13,0),
        halign= :left,
        )

    if ind[2] == 1
        ax.xlabel = "L2FC abundance"
    else
        ax.xlabel = "L2FC half-life"
    end

    if ind[2] == 1
        ax.ylabel = "-log10(p)"
    end
    push!(volcano_axes,ax)
end

curve_axes = []
for (ind,letter) ∈ zip([(1,3),(1,4),(3,3),(3,4)],["C","D","G","H"])
    ax = Axis(fig[ind[1],ind[2]],
    width = 200,
    height = 200)
    Label(fig[ind[1],ind[2],TopLeft()],letter,
        font = :bold,
        fontsize=18,
        padding = (0,0,13,0),
        halign= :left
        )
    ax.xlabel = "day"
    if ind[2] == 1
        ax.ylabel = "proportion labeled"
    end
    push!(curve_axes,ax)
end

Label(fig[0,:],"young (9 weeks)", font = :bold,
    fontsize=18
    )
Label(fig[2,:],"old (63 weeks)", font = :bold,
    fontsize=18
    )

ages = ["young", "aged"]
    
for (age_idx, this_age) in enumerate(ages)
    age_data = filter(row -> row.age == this_age, genotypic_diffs)
    
    n_tests = nrow(age_data)
    age_data[!, :total_color_alpha] = map(age_data.pval_total) do p
        p > 0.05 / n_tests ? 0.1 : 1.0
    end
    age_data[!, :label_ffi_color_alpha] = map(age_data.pval_label_ffi) do p
        p > 0.05 / n_tests ? 0.1 : 1.0
    end

    # Abundance volcano
    volcano_idx = (age_idx - 1) * 2 + 1
    ax_abun = volcano_axes[volcano_idx]
    ax_abun.yticks = (0:10, [i < 10 ? "$i" : "\u226510" for i in 0:10])
    ax_abun.xticks = (-3:3)
    xlims!(ax_abun, -3, 3)
    ylims!(ax_abun, 0, 10.5)
    
    hlines!(ax_abun, [-log10(0.05 / n_tests)], color=:black, linewidth=1,linestyle = :dash )
    
    for row in eachrow(age_data)
        col = RGBA(row.total_color, row.total_color_alpha)
        marker_shape = row.peptide == "GENFTETDVK" ? :diamond : row.peptide == "VVEQMCVTQYQK" ? :rect : :circle
        scatter!(ax_abun, [row.log2_ratio_total], [min(-log10(row.pval_total), 10)],
                color=col, markersize=6, marker=marker_shape)
    end
    
    prnp_data = filter(row -> row.protein == "PRNP", age_data)
    
    # Plot all PRNP data points
    for row in eachrow(prnp_data)
        marker_shape = row.peptide == "GENFTETDVK" ? :diamond : row.peptide == "VVEQMCVTQYQK" ? :rect : :circle
        scatter!(ax_abun, [row.log2_ratio_total], [min(-log10(row.pval_total), 10)],
                color=:transparent, strokecolor=:black, strokewidth=1.5, markersize=8, marker=marker_shape)
    end

    volcano_idx = (age_idx - 1) * 2 + 2
    ax_hl = volcano_axes[volcano_idx]
    ax_hl.yticks = (0:10, [i < 10 ? "$i" : "\u226510" for i in 0:10])
    ax_hl.xticks = (-3:3)
    xlims!(ax_hl, -3, 3)
    ylims!(ax_hl, 0, 10.5)
    
    hlines!(ax_hl, [-log10(0.05 / n_tests)], color=:black, linewidth=1,linestyle=:dash)
    
    for row in eachrow(age_data)
        col = RGBA(row.label_ffi_color, row.label_ffi_color_alpha)
        marker_shape = row.peptide == "GENFTETDVK" ? :diamond : row.peptide == "VVEQMCVTQYQK" ? :rect : :circle
        scatter!(ax_hl, [row.log2_ratio_thalf], [min(-log10(row.pval_label_ffi), 10)],
                color=col, markersize=6, marker=marker_shape)
    end
    
    for row in eachrow(prnp_data)
        marker_shape = row.peptide == "GENFTETDVK" ? :diamond : row.peptide == "VVEQMCVTQYQK" ? :rect : :circle
        scatter!(ax_hl, [row.log2_ratio_thalf], [min(-log10(row.pval_label_ffi), 10)],
                color=:transparent, strokecolor=:black, strokewidth=1.5, markersize=8, marker=marker_shape)
    end
    

    for pep_idx in 1:2
        this_peptide_full = full_peptides[pep_idx]
        this_peptide_short = use_peptides[pep_idx]
        curve_idx = (age_idx - 1) * 2 + pep_idx
        ax = curve_axes[curve_idx]
        ax.xlabel = "day"
        ax.ylabel = pep_idx == 1 ? "proportion labeled" : ""
        ax.title = this_peptide_short
        
        
        xlims!(ax, -0.5, 8.5)
        ylims!(ax, 0, 0.50)

        t_fine = collect(0:0.1:8)
        
        for geno_row in eachrow(leg)
            curve_data = filter(row -> row.genotype == geno_row.genotype &&
                                        row.peptide == this_peptide_full &&
                                        row.age == this_age, ffi_all)

            if nrow(curve_data) == 0
                continue
            end

            summary_df = combine(groupby(curve_data, :chow_days)) do df
                DataFrame(
                    mean_prop = mean(df.prop_labeled),
                    l95 = lower(df.prop_labeled),
                    u95 = upper(df.prop_labeled)
                )
            end
            sort!(summary_df, :chow_days)

            mixture_fit = fit_isotopic_thalf_mixture(
                curve_data.chow_days,
                curve_data.prop_labeled,
                free_lysine
            )

            if !isnan(mixture_fit.effective_thalf)
                pred = proportion_labeled_mixture(
                    mixture_fit.thalves,
                    mixture_fit.proportions,
                    t_fine,
                    free_lysine
                )
                lines!(ax, t_fine, pred, color=geno_row.color, linewidth=2, linestyle=:dash)
            end

            band!(ax, summary_df.chow_days, summary_df.l95, summary_df.u95,
                    color=RGBA(geno_row.color, 0.35))
            p = scatter!(ax, curve_data.chow_days, curve_data.prop_labeled,
                    color=geno_row.color, markersize=6, strokewidth=0.5,
                    strokecolor=geno_row.color)
            
        end
    end
end

linkaxes!(volcano_axes...)
linkaxes!(curve_axes...)

leg_colors = [colorant"#D95F02", colorant"#22127A", colorant"#77127A"];
leg_labels = ["ki-3F4-FFI", "ki-3F4-WT", "C57BL/6N WT"];

# Create GridLayout for legend column
legend_grid = fig[:,5] = GridLayout()

# Peptide legend (above genotype legend)
peptide_elements = [
    MarkerElement(color = :black, marker = :diamond, markersize = 10, strokecolor = :black),
    MarkerElement(color = :gray, marker = :rect, markersize = 10, strokecolor = :black)
]
peptide_labels = ["GENF", "VVEQ"]

Legend(legend_grid[1,1], peptide_elements, peptide_labels, "Peptides", framevisible = false, valign = :bottom, halign = :left, titleposition = :left, titleorientation = :vertical)

# Genotype legend
legend_elements = [
    [MarkerElement(color = c, marker = :circle, markersize = 10, strokecolor = c, strokewidth = 0.5),
     LineElement(color = c, linestyle = :dash, linewidth = 2)]
    for c in leg_colors
];

Legend(legend_grid[2,1], legend_elements, leg_labels, "Curves", framevisible = false, valign = :top, halign = :left, titleposition = :left, titleorientation = :vertical)
colsize!(fig.layout,5,Aspect(1,1.0))
resize_to_layout!(fig)
save("test.png",fig,px_per_unit=3)


generate_mixture_model_table()