"""
Fit Mixture Model to CRL Isotope Labeling Data

This fits the mixture decay model with lysine labeling dynamics to
actual proteomics data measuring proportion of heavy-labeled protein over time.
"""

using CSV, DataFrames, Statistics, Distributions
include("mixture_decay_optim.jl")
include("demo_labeling.jl")

"""
    calculate_confidence_intervals(x; ci=0.95)

Calculate mean and 95% CI for a vector
"""
function calculate_confidence_intervals(x; ci=0.95)
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

"""
    read_crl_data(filepath)

Read CRL proteomics labeling data.
Expected columns: day, proportion_heavy, tissue
"""
function read_crl_data(filepath)
    if !isfile(filepath)
        error("Data file not found: $filepath\nPlease provide path to CRL proteomics data CSV")
    end

    df = CSV.read(filepath, DataFrame)

    # Check required columns
    required_cols = [:day, :proportion_heavy, :tissue]
    for col in required_cols
        if !(col in names(df))
            error("Missing required column: $col")
        end
    end

    return df
end

"""
    summarize_by_tissue(df)

Group by tissue and day, calculate means and confidence intervals
"""
function summarize_by_tissue(df)
    grouped = combine(groupby(df, [:tissue, :day])) do subdf
        ci = calculate_confidence_intervals(subdf.proportion_heavy)
        (mean = ci.mean, l95 = ci.lower, u95 = ci.upper, n = nrow(subdf))
    end

    return sort(grouped, [:tissue, :day])
end

"""
    fit_tissue_labeling_data(df_tissue, tissue_name; constrained=false)

Fit mixture model with labeling to a single tissue's data
"""
function fit_tissue_labeling_data(df_tissue, tissue_name; constrained=false)
    println("\n" * "="^70)
    println("Fitting tissue: $tissue_name")
    println("="^70)

    # Extract data
    days = sort(unique(df_tissue.day))

    # For each day, get mean proportion labeled
    P_data = [mean(df_tissue[df_tissue.day .== d, :proportion_heavy]) for d in days]

    # RNA is constant (no knockdown in labeling experiment)
    R_data = ones(length(days))

    println("\nData points: $(length(days))")
    println("Day range: $(minimum(days)) to $(maximum(days))")
    println("Proportion labeled range: $(round(minimum(P_data), digits=3)) to $(round(maximum(P_data), digits=3))")

    # Fit single-rate model
    println("\n1. Fitting single-rate model...")
    λ_single, res_single = fit_single_rate_model(days, R_data, P_data)

    # Fit 2-component mixture (unconstrained)
    println("\n2. Fitting 2-component mixture...")
    λ_mixture, ψ_vec, π_vec, res_mixture = fit_mixture_model_with_labeling(
        days, R_data, P_data, 2, constrained=false
    )

    # Fit 2-component mixture (constrained)
    if constrained
        println("\n3. Fitting 2-component mixture (CONSTRAINED: ∑ψₖ = λ)...")
        λ_constrained, ψ_vec_constrained, π_vec_constrained, res_constrained =
            fit_mixture_model_with_labeling(days, R_data, P_data, 2, constrained=true)
    end

    # Model comparison
    println("\n4. Model comparison:")
    println("  Single-rate RSS: $(round(res_single.minimum, digits=6))")
    println("  Mixture RSS:     $(round(res_mixture.minimum, digits=6))")
    println("  Improvement:     $(round((1 - res_mixture.minimum/res_single.minimum)*100, digits=1))%")

    # Calculate AIC
    n = length(P_data)
    k_single = 1
    k_mixture = 4  # λ + 2 ψs + 1 π

    AIC_single = n * log(res_single.minimum / n) + 2 * k_single
    AIC_mixture = n * log(res_mixture.minimum / n) + 2 * k_mixture

    println("\n5. Information criteria:")
    println("  Single-rate AIC: $(round(AIC_single, digits=2))")
    println("  Mixture AIC:     $(round(AIC_mixture, digits=2))")
    println("  ΔAIC:            $(round(AIC_single - AIC_mixture, digits=2))")

    if AIC_single - AIC_mixture > 4
        println("  → Strong evidence for mixture model (ΔAIC > 4)")
    elseif AIC_single - AIC_mixture > 2
        println("  → Moderate evidence for mixture model (ΔAIC > 2)")
    else
        println("  → Weak evidence for mixture model")
    end

    return (
        days = days,
        P_data = P_data,
        R_data = R_data,
        λ_single = λ_single,
        λ_mixture = λ_mixture,
        ψ_vec = ψ_vec,
        π_vec = π_vec,
        AIC_single = AIC_single,
        AIC_mixture = AIC_mixture
    )
end

"""
    plot_tissue_fits(results, tissue_name)

Create comparison plot for a tissue
"""
function plot_tissue_fits(results, tissue_name)
    days = results.days
    P_data = results.P_data
    R_data = results.R_data
    λ_single = results.λ_single
    λ_mixture = results.λ_mixture
    ψ_vec = results.ψ_vec
    π_vec = results.π_vec

    # Generate fine time grid for predictions
    t_fine = range(0, maximum(days), length=1000)
    R_interp = LinearInterpolation(days, R_data, extrapolation_bc=Line())

    # Lysine availability
    avails_fine = free_lysine.(t_fine)
    avails_interp = LinearInterpolation(t_fine, avails_fine, extrapolation_bc=Line())

    # Single rate prediction
    P_single = solve_protein_dynamics(R_interp, t_fine, λ_single)

    # Mixture prediction with labeling
    P_obs_mixture, P_obs_sub, P_total_mixture, L_sub =
        solve_mixture_dynamics_with_labeling(R_interp, avails_interp, t_fine, λ_mixture, ψ_vec, π_vec)

    # Create figure
    fig = Figure(size=(1200, 900), fontsize=14)

    # Panel A: Data and fits
    ax1 = Axis(fig[1, 1],
        xlabel="Day",
        ylabel="Proportion Labeled",
        title="A. $tissue_name: Model Fits")

    scatter!(ax1, days, P_data, label="Data", color=:black, markersize=10)
    lines!(ax1, t_fine, P_single .* avails_fine, label="Single Rate",
           color=:red, linestyle=:dash, linewidth=2)
    lines!(ax1, t_fine, P_obs_mixture, label="Mixture (K=2)",
           color=:blue, linewidth=3)

    axislegend(ax1, position=:rb)

    # Panel B: Lysine availability
    ax2 = Axis(fig[1, 2],
        xlabel="Day",
        ylabel="Availability",
        title="B. Free Lysine Availability")

    lines!(ax2, t_fine, avails_fine, color=:green, linewidth=3)
    hlines!(ax2, [avails_fine[end]], color=:gray, linestyle=:dash, label="Asymptote")

    axislegend(ax2, position=:rb)

    # Panel C: Total protein (unlabeled + labeled)
    ax3 = Axis(fig[2, 1],
        xlabel="Day",
        ylabel="Total Protein",
        title="C. Total Protein (Both Labeled & Unlabeled)")

    lines!(ax3, t_fine, P_total_mixture, color=:blue, linewidth=3, label="Total")

    for k in 1:2
        t_half_k = log(2) / ψ_vec[k]
        lines!(ax3, t_fine, P_total_mixture .* π_vec[k],
               label="Component $k (t₁/₂=$(round(t_half_k, digits=1))d)",
               linestyle=:dash, linewidth=2)
    end

    axislegend(ax3, position=:rb)

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

    return fig
end

"""
    main_analysis(data_filepath)

Main analysis pipeline
"""
function main_analysis(data_filepath)
    println("="^70)
    println("CRL Isotope Labeling Data Analysis")
    println("="^70)

    # Read data
    println("\nReading data from: $data_filepath")
    df = read_crl_data(data_filepath)

    println("Data summary:")
    println("  Total observations: $(nrow(df))")
    println("  Tissues: $(join(unique(df.tissue), ", "))")
    println("  Day range: $(minimum(df.day)) to $(maximum(df.day))")

    # Summarize by tissue
    println("\nCalculating summary statistics...")
    df_summary = summarize_by_tissue(df)

    # Fit each tissue separately
    tissues = unique(df.tissue)
    results_by_tissue = Dict()

    for tissue in tissues
        df_tissue = df[df.tissue .== tissue, :]
        results = fit_tissue_labeling_data(df_tissue, tissue, constrained=false)
        results_by_tissue[tissue] = results

        # Plot
        fig = plot_tissue_fits(results, tissue)
        save("$(tissue)_labeling_fit.png", fig)
        println("\nSaved: $(tissue)_labeling_fit.png")
    end

    # Summary comparison
    println("\n" * "="^70)
    println("SUMMARY: All Tissues")
    println("="^70)

    for tissue in tissues
        results = results_by_tissue[tissue]
        println("\n$tissue:")
        println("  Single-rate t₁/₂: $(round(log(2)/results.λ_single, digits=2)) days")
        println("  Mixture component 1: t₁/₂=$(round(log(2)/results.ψ_vec[1], digits=2))d, π=$(round(results.π_vec[1], digits=2))")
        println("  Mixture component 2: t₁/₂=$(round(log(2)/results.ψ_vec[2], digits=2))d, π=$(round(results.π_vec[2], digits=2))")
        println("  ΔAIC: $(round(results.AIC_single - results.AIC_mixture, digits=2))")
    end

    return df_summary, results_by_tissue
end

# Usage instructions
println("""
Usage:
------
1. Prepare your data as CSV with columns: day, proportion_heavy, tissue
2. Run: df_summary, results = main_analysis("path/to/your/data.csv")

Example data format:
day,proportion_heavy,tissue
0,0.01,brain
0,0.02,brain
1,0.15,brain
...

If you don't have the data file, you can create synthetic data:
>>> include("demo_labeling.jl")
""")
