include("TurnoverModels.jl")
using .TurnoverModels
using DataFrames, CSV, Printf, Statistics

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

function main()
    println("="^80)
    println("FIG2 CRL PROTEOMICS MIXTURE MODEL ANALYSIS")
    println("="^80)
    println()

    # Load CRL proteomics data
    println("Loading CRL proteomics data...")
    raw_data = CSV.read("data/crl_proteomics.tsv", DataFrame, types=Dict(:light => Float64, :heavy => Union{Float64, Missing}))

    # Load peptide standards to determine BQL
    pep_stds = CSV.read("data/crl_pepstds.tsv", DataFrame)
    bql = pep_stds.back_calc_ng_ml[1] / 10
    println("BQL value: $bql")

    # Handle missing heavy values (set to BQL)
    raw_data.heavy = coalesce.(raw_data.heavy, bql)

    # Calculate total and proportion labeled
    raw_data.total = raw_data.light .+ raw_data.heavy
    raw_data.prop_labeled = raw_data.heavy ./ raw_data.total

    # Subtract LLQ from all values, floor at 0
    println("\nApplying baseline subtraction...")
    println("="^80)
    for tissue in unique(raw_data.tissue)
        tissue_mask = raw_data.tissue .== tissue

        # Find LLQ: max proportion at day 0 where heavy == bql
        day0_bql_mask = tissue_mask .& (raw_data.day .== 0) .& (raw_data.heavy .== bql)

        if any(day0_bql_mask)
            llq_value = maximum(raw_data.prop_labeled[day0_bql_mask])
        else
            llq_value = 0.0
        end

        # Calculate statistics before adjustment
        before_mean = mean(raw_data.prop_labeled[tissue_mask])
        before_min = minimum(raw_data.prop_labeled[tissue_mask])
        before_max = maximum(raw_data.prop_labeled[tissue_mask])

        println("\n$tissue:")
        println("  LLQ value (from day 0 BQL samples): $(round(llq_value, digits=6))")
        println("  Before adjustment:")
        println("    Range: $(round(before_min, digits=4)) to $(round(before_max, digits=4))")
        println("    Mean: $(round(before_mean, digits=4))")

        # Subtract LLQ and floor at 0
        raw_data.prop_labeled[tissue_mask] = max.(0.0, raw_data.prop_labeled[tissue_mask] .- llq_value)

        # Calculate statistics after adjustment
        after_mean = mean(raw_data.prop_labeled[tissue_mask])
        after_min = minimum(raw_data.prop_labeled[tissue_mask])
        after_max = maximum(raw_data.prop_labeled[tissue_mask])
        n_zeros = sum(raw_data.prop_labeled[tissue_mask] .== 0.0)

        println("  After adjustment:")
        println("    Range: $(round(after_min, digits=4)) to $(round(after_max, digits=4))")
        println("    Mean: $(round(after_mean, digits=4))")
        println("    Data points floored to 0: $n_zeros")
    end
    println("\n" * "="^80)

    println("  Loaded $(nrow(raw_data)) observations")
    println("  Tissues: $(unique(raw_data.tissue))")
    println("  Time points: $(sort(unique(raw_data.day)))")

    # Create output directory
    println("\nCreating output directory...")
    mkpath("data/crl")

    # Save baseline-subtracted data for verification/plotting
    println("Saving baseline-subtracted data...")
    adjusted_data_output = select(raw_data, :animal, :day, :tissue, :light, :heavy, :total, :prop_labeled)
    CSV.write("data/crl/baseline_subtracted_data.tsv", adjusted_data_output, delim='\t')
    println("  Saved to: data/crl/baseline_subtracted_data.tsv")

    results = DataFrame()

    # Analyze each tissue separately
    println("\n" * "="^80)
    println("FITTING MODELS ON BASELINE-SUBTRACTED DATA")
    println("="^80)
    for tissue in unique(raw_data.tissue)
        println("\nAnalyzing tissue: $tissue")
        tissue_data = filter(row -> row.tissue == tissue, raw_data)

        days = tissue_data.day
        props = tissue_data.prop_labeled  # Already baseline-subtracted and floored at 0

        println("  Using $(length(days)) observations with adjusted prop_labeled values")

        # Fit single-rate model
        thalf_single = fit_single_isotopic(days, props)
        println("  Single-rate t½: $thalf_single days")

        # Fit mixture model
        mix = fit_mixture_isotopic(days, props)
        println("  Mixture model:")
        println("    Fast: $(mix.thalves[1]) d ($(mix.proportions[1]*100)%)")
        println("    Slow: $(mix.thalves[2]) d ($(mix.proportions[2]*100)%)")
        println("    Effective: $(mix.effective_thalf) d")

        # Calculate model comparison metrics
        n = length(days)

        # Single model AIC/BIC
        if !isnan(thalf_single)
            pred_single = proportion_labeled_single(thalf_single, days)
            rss_single = sum((props .- pred_single).^2)
            k_single = 2  # thalf, sigma
            aic_single = calculate_aic(rss_single, n, k_single)
            bic_single = calculate_bic(rss_single, n, k_single)
        else
            aic_single = NaN
            bic_single = NaN
        end

        # Mixture model AIC/BIC
        if !isnan(mix.thalves[1])
            pred_mixture = proportion_labeled_mixture(mix.thalves, mix.proportions, days)
            rss_mixture = sum((props .- pred_mixture).^2)
            k_mixture = 4  # t1, t2, pi1, sigma
            aic_mixture = calculate_aic(rss_mixture, n, k_mixture)
            bic_mixture = calculate_bic(rss_mixture, n, k_mixture)
        else
            aic_mixture = NaN
            bic_mixture = NaN
        end

        # Calculate delta values
        delta_aic = isnan(aic_single) ? NaN : aic_mixture - aic_single
        delta_bic = isnan(bic_single) ? NaN : bic_mixture - bic_single

        println("  ΔAIC: $delta_aic")
        println("  ΔBIC: $delta_bic")

        push!(results, (
            tissue = tissue,
            n_obs = n,

            # Single model
            thalf_single = thalf_single,
            aic_single = aic_single,
            bic_single = bic_single,

            # Mixture model
            thalf_fast = mix.thalves[1],
            thalf_slow = mix.thalves[2],
            prop_fast = mix.proportions[1],
            prop_slow = mix.proportions[2],
            thalf_effective = mix.effective_thalf,
            aic_mixture = aic_mixture,
            bic_mixture = bic_mixture,

            # Comparison
            delta_aic = delta_aic,
            delta_bic = delta_bic,
            model_ratio = isnan(thalf_single) ? NaN : mix.effective_thalf / thalf_single
        ))
    end

    # Save results
    println("\nSaving model fitting results...")
    CSV.write("data/crl/mixture_params.tsv", results, delim='\t')
    println("  Saved to: data/crl/mixture_params.tsv")

    # Generate report
    generate_report(results)

    println("Analysis complete!")
end

function generate_report(results)
    println("\n" * "="^80)
    println("CRL PROTEOMICS MIXTURE MODEL ANALYSIS REPORT")
    println("="^80)

    for row in eachrow(results)
        println("\n$(uppercase(row.tissue)):")
        @printf("  n = %d observations\n", row.n_obs)

        if !isnan(row.thalf_single)
            @printf("  Single-rate model:\n")
            @printf("    t½ = %.2f days\n", row.thalf_single)
            @printf("    AIC = %.1f, BIC = %.1f\n", row.aic_single, row.bic_single)
        else
            println("  Single-rate model: FAILED")
        end

        if !isnan(row.thalf_fast)
            @printf("  Mixture model:\n")
            @printf("    Fast pool: t½ = %.2f days (%.1f%%)\n", row.thalf_fast, row.prop_fast*100)
            @printf("    Slow pool: t½ = %.2f days (%.1f%%)\n", row.thalf_slow, row.prop_slow*100)
            @printf("    Effective t½ = %.2f days\n", row.thalf_effective)
            @printf("    AIC = %.1f, BIC = %.1f\n", row.aic_mixture, row.bic_mixture)

            if !isnan(row.delta_aic)
                @printf("  Model comparison:\n")
                @printf("    ΔAIC = %.1f", row.delta_aic)
                if row.delta_aic < -3
                    print(" (mixture substantially better)")
                elseif row.delta_aic < 0
                    print(" (mixture slightly better)")
                else
                    print(" (single-rate preferred)")
                end
                println()
                @printf("    ΔBIC = %.1f", row.delta_bic)
                if row.delta_bic < -3
                    print(" (strong evidence for mixture)")
                elseif row.delta_bic < 0
                    print(" (weak evidence for mixture)")
                else
                    print(" (no evidence for mixture)")
                end
                println()
            end
        else
            println("  Mixture model: FAILED")
        end
    end

    println("\n" * "="^80 * "\n")
end

# Run main analysis
main()
