"""
Analyze CRL Proteomics Data - Direct Translation from R Script

This replicates the CRL data reading and processing from src/halflife_figures.R
Fits mixture models with lysine labeling to brain and colon data separately.
"""

using CSV, DataFrames, Statistics, Distributions
include("mixture_decay_optim.jl")

# ============================================================================
# Helper Functions (from R script)
# ============================================================================

"""Calculate confidence intervals"""
function calculate_ci(x; ci=0.95)
    n = length(x)
    if n == 0
        return (mean=NaN, lower=NaN, upper=NaN)
    end

    m = mean(x)
    if n == 1
        return (mean=m, lower=m, upper=m)
    end

    α = 1 - ci
    sds = quantile(Normal(), 1 - α/2)
    se = std(x) / sqrt(n)

    return (mean=m, lower=m - sds*se, upper=m + sds*se)
end

# ============================================================================
# Read and Process CRL Data (matching R script lines 621-631)
# ============================================================================

println("="^70)
println("CRL Proteomics Data Analysis")
println("="^70)

# Read peptide standards (line 621)
println("\n1. Reading data...")
pep_stds = CSV.read("data/crl_pepstds.tsv", DataFrame, types=Dict(:back_calc_ng_ml=>Float64))
bql = pep_stds[1, :back_calc_ng_ml] / 10  # Below quantification limit
println("  BQL threshold: $(round(bql, digits=4)) ng/mL")

# Tissue metadata (lines 623-625)
tissue_meta = DataFrame(
    tissue = ["brain", "colon"],
    color = ["#CCCF00", "#EEBB77"],
    x = [1, 2]
)

# Read CRL proteomics data (lines 627-631)
# Explicitly specify column types to ensure numeric parsing
crl_prot = CSV.read("data/crl_proteomics.tsv", DataFrame,
                    types=Dict(:animal=>Int, :day=>Float64, :light=>Float64,
                              :heavy=>Union{Float64,Missing}, :tissue=>String))

# Process exactly like R script:
# - Replace NA heavy with BQL
# - Calculate total = light + heavy
# - Calculate proportion_heavy = heavy / total
# - Join with tissue metadata

crl_prot[!, :heavy] = coalesce.(crl_prot.heavy, bql)
crl_prot[!, :total] = crl_prot.light .+ crl_prot.heavy
crl_prot[!, :proportion_heavy] = crl_prot.heavy ./ crl_prot.total
crl_prot = innerjoin(crl_prot, tissue_meta, on=:tissue)

println("  Total observations: $(nrow(crl_prot))")
println("  Tissues: $(join(unique(crl_prot.tissue), ", "))")
println("  Days: $(join(sort(unique(crl_prot.day)), ", "))")
println("  Proportion heavy range: $(round(minimum(crl_prot.proportion_heavy), digits=3)) to $(round(maximum(crl_prot.proportion_heavy), digits=3))")

# ============================================================================
# Summarize by Tissue and Day (matching R script lines 904-909)
# ============================================================================

println("\n2. Calculating summary statistics by tissue...")

crl_summary = combine(groupby(crl_prot, [:tissue, :day, :color])) do subdf
    ci = calculate_ci(subdf.proportion_heavy)
    (mean = ci.mean, l95 = ci.lower, u95 = ci.upper, n = nrow(subdf))
end

crl_summary = sort(crl_summary, [:tissue, :day])

println("\nSummary by tissue:")
for tissue in ["brain", "colon"]
    tissue_data = crl_summary[crl_summary.tissue .== tissue, :]
    println("\n$tissue:")
    for row in eachrow(tissue_data)
        println("  Day $(row.day): $(round(row.mean, digits=3)) [$(round(row.l95, digits=3)), $(round(row.u95, digits=3))], n=$(row.n)")
    end
end

# Save summary table
CSV.write("crl_summary.csv", crl_summary)
println("\nSaved summary: crl_summary.csv")

# ============================================================================
# Calculate Effective BQL by Tissue (matching R script lines 898-903)
# ============================================================================

println("\n3. Calculating effective BQL (below quantification limit) for each tissue...")

effective_bql = Dict()
for tissue in ["brain", "colon"]
    tissue_day0 = crl_prot[(crl_prot.tissue .== tissue) .& (crl_prot.day .== 0) .& (crl_prot.heavy .== bql), :]
    if nrow(tissue_day0) > 0
        effective_bql[tissue] = maximum(tissue_day0.proportion_heavy)
        println("  $tissue LLQ: $(round(effective_bql[tissue], digits=4))")
    else
        effective_bql[tissue] = bql / (bql + mean(crl_prot[crl_prot.tissue .== tissue, :light]))
        println("  $tissue LLQ (estimated): $(round(effective_bql[tissue], digits=4))")
    end
end

# ============================================================================
# Fit Mixture Models to Each Tissue
# ============================================================================

println("\n" * "="^70)
println("FITTING MIXTURE MODELS")
println("="^70)

results_by_tissue = Dict()

for tissue in ["brain", "colon"]
    println("\n" * "="^70)
    println("Tissue: $tissue")
    println("="^70)

    # Extract data for this tissue
    tissue_data = crl_prot[crl_prot.tissue .== tissue, :]

    # Get unique days and calculate mean proportion for each
    days = sort(unique(tissue_data.day))
    P_data = [mean(tissue_data[tissue_data.day .== d, :proportion_heavy]) for d in days]

    # RNA is constant (no knockdown in labeling experiment)
    R_data = ones(length(days))

    println("\nData summary:")
    println("  Time points: $(length(days))")
    println("  Days: $(join(days, ", "))")
    println("  Proportion labeled range: $(round(minimum(P_data), digits=3)) to $(round(maximum(P_data), digits=3))")
    println("  Effective BQL: $(round(effective_bql[tissue], digits=4))")

    # Fit single-rate model
    println("\n1. Fitting single-rate model...")
    λ_single, res_single = fit_single_rate_model(days, R_data, P_data)

    # Fit 2-component mixture (unconstrained) with labeling
    println("\n2. Fitting 2-component mixture WITH LABELING (unconstrained)...")
    λ_mixture, ψ_vec, π_vec, res_mixture = fit_mixture_model_with_labeling(
        days, R_data, P_data, 2, constrained=false
    )

    # Model comparison
    println("\n3. Model comparison:")
    println("  Single-rate RSS: $(round(res_single.minimum, digits=6))")
    println("  Mixture RSS:     $(round(res_mixture.minimum, digits=6))")
    improvement = (1 - res_mixture.minimum/res_single.minimum) * 100
    println("  Improvement:     $(round(improvement, digits=1))%")

    # Calculate AIC
    n = length(P_data)
    k_single = 1
    k_mixture = 4  # λ + 2 ψs + 1 π

    AIC_single = n * log(res_single.minimum / n) + 2 * k_single
    AIC_mixture = n * log(res_mixture.minimum / n) + 2 * k_mixture

    ΔAIC = AIC_single - AIC_mixture

    println("\n4. Information criteria:")
    println("  Single-rate AIC: $(round(AIC_single, digits=2))")
    println("  Mixture AIC:     $(round(AIC_mixture, digits=2))")
    println("  ΔAIC:            $(round(ΔAIC, digits=2))")

    if ΔAIC > 4
        println("  → STRONG evidence for mixture model (ΔAIC > 4)")
    elseif ΔAIC > 2
        println("  → Moderate evidence for mixture model (ΔAIC > 2)")
    else
        println("  → Weak evidence for mixture model")
    end

    # Store results
    results_by_tissue[tissue] = (
        days = days,
        P_data = P_data,
        R_data = R_data,
        λ_single = λ_single,
        λ_mixture = λ_mixture,
        ψ_vec = ψ_vec,
        π_vec = π_vec,
        AIC_single = AIC_single,
        AIC_mixture = AIC_mixture,
        ΔAIC = ΔAIC,
        improvement = improvement
    )

    # Generate predictions for plotting
    println("\n5. Generating predictions...")
    t_fine = range(0, maximum(days), length=1000)
    R_interp = LinearInterpolation(days, R_data, extrapolation_bc=Line())

    # Lysine availability
    avails_fine = free_lysine.(t_fine)
    avails_interp = LinearInterpolation(t_fine, avails_fine, extrapolation_bc=Line())

    # Single rate prediction with labeling
    P_single_total = solve_protein_dynamics(R_interp, t_fine, λ_single)
    P_single_obs = P_single_total .* avails_fine

    # Mixture prediction with labeling
    P_obs_mixture, P_obs_sub, P_total_mixture, L_sub =
        solve_mixture_dynamics_with_labeling(R_interp, avails_interp, t_fine, λ_mixture, ψ_vec, π_vec)

    # Create figure
    fig = Figure(size=(1400, 1000), fontsize=14)

    # Panel A: Data and model fits
    ax1 = Axis(fig[1, 1],
        xlabel="Day",
        ylabel="Proportion Labeled",
        title="A. $tissue: Model Fits to Data")

    scatter!(ax1, days, P_data, label="Data", color=:black, markersize=12)
    lines!(ax1, t_fine, P_single_obs, label="Single Rate",
           color=:red, linestyle=:dash, linewidth=3)
    lines!(ax1, t_fine, P_obs_mixture, label="Mixture (K=2)",
           color=:blue, linewidth=3)
    hlines!(ax1, [effective_bql[tissue]], color=:gray, linestyle=:dot, label="LLQ")

    axislegend(ax1, position=:rb)

    # Panel B: Free lysine availability
    ax2 = Axis(fig[1, 2],
        xlabel="Day",
        ylabel="Lysine Availability",
        title="B. Free Lysine Dynamics")

    lines!(ax2, t_fine, avails_fine, color=:green, linewidth=3, label="Free lysine")
    hlines!(ax2, [avails_fine[end]], color=:gray, linestyle=:dash, label="Asymptote")

    axislegend(ax2, position=:rb)

    # Panel C: Total protein (both labeled and unlabeled)
    ax3 = Axis(fig[2, 1],
        xlabel="Day",
        ylabel="Total Protein",
        title="C. Total Protein (Labeled + Unlabeled)")

    lines!(ax3, t_fine, P_total_mixture, color=:blue, linewidth=3, label="Total")

    for k in 1:2
        t_half_k = log(2) / ψ_vec[k]
        amp_k = λ_mixture / ψ_vec[k]
        P_k = [sum(P_total_mixture[i] * π_vec[k] for i in 1:length(P_total_mixture))]
        lines!(ax3, t_fine, P_total_mixture .* π_vec[k],
               label="Comp $k (t₁/₂=$(round(t_half_k, digits=1))d, amp=$(round(amp_k, digits=2)))",
               linestyle=:dash, linewidth=2)
    end

    axislegend(ax3, position=:rt)

    # Panel D: Labeled fractions
    ax4 = Axis(fig[2, 2],
        xlabel="Day",
        ylabel="Labeled Fraction",
        title="D. Labeling Dynamics by Component")

    for k in 1:2
        t_half_k = log(2) / ψ_vec[k]
        lines!(ax4, t_fine, L_sub[:, k],
               label="Component $k (t₁/₂=$(round(t_half_k, digits=1))d)",
               linewidth=3)
    end

    lines!(ax4, t_fine, avails_fine, color=:green, linestyle=:dot, linewidth=2, label="Lysine avail")

    axislegend(ax4, position=:rb)

    # Panel E: Observed protein (only labeled)
    ax5 = Axis(fig[3, 1:2],
        xlabel="Day",
        ylabel="Observed Protein (Proportion Labeled)",
        title="E. Observed = Total × Labeled Fraction")

    scatter!(ax5, days, P_data, label="Data", color=:black, markersize=12)
    lines!(ax5, t_fine, P_obs_mixture, label="Mixture total", color=:blue, linewidth=3)

    for k in 1:2
        lines!(ax5, t_fine, P_obs_sub[:, k],
               label="Component $k (observed)",
               linestyle=:dash, linewidth=2)
    end

    axislegend!(ax5, position=:rb)

    # Save figure
    save("crl_$(tissue)_mixture_fit.png", fig)
    println("\nSaved: crl_$(tissue)_mixture_fit.png")
end

# ============================================================================
# Final Summary
# ============================================================================

println("\n" * "="^70)
println("FINAL SUMMARY: All Tissues")
println("="^70)

for tissue in ["brain", "colon"]
    results = results_by_tissue[tissue]
    println("\n$tissue:")
    println("  Single-rate:")
    println("    t₁/₂: $(round(log(2)/results.λ_single, digits=2)) days")
    println("    RSS: $(round(results.AIC_single, digits=2))")

    println("  Mixture model:")
    for k in 1:2
        t_half_k = log(2) / results.ψ_vec[k]
        amp_k = results.λ_mixture / results.ψ_vec[k]
        println("    Component $k: t₁/₂=$(round(t_half_k, digits=1))d, amp=$(round(amp_k, digits=2)), π=$(round(results.π_vec[k], digits=2))")
    end
    println("    RSS: $(round(results.AIC_mixture, digits=2))")

    println("  Model selection:")
    println("    ΔAIC: $(round(results.ΔAIC, digits=2))")
    println("    Improvement: $(round(results.improvement, digits=1))%")

    if results.ΔAIC > 4
        println("    → HETEROGENEOUS population detected!")
    else
        println("    → Homogeneous population (single rate sufficient)")
    end
end

println("\n" * "="^70)
println("KEY FINDINGS")
println("="^70)
println("\nWith lysine labeling dynamics:")
println("1. Observed protein = Total protein × Labeled fraction")
println("2. Labeled fraction approaches lysine availability (~60%)")
println("3. Fast-degrading proteins label quickly (high turnover)")
println("4. Slow-degrading proteins label slowly (low turnover)")
println("\nIf ΔAIC > 4: Strong evidence for heterogeneous decay rates!")
println("This indicates mixture of fast/slow degrading subpopulations.")

println("\n" * "="^70)
println("Analysis complete!")
println("="^70)
