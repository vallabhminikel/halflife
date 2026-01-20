include("TurnoverModels.jl")
using .TurnoverModels
using DataFrames, CSV, Statistics, Printf, Interpolations

"""
Load and prepare Panel A data: ASO6 WT mice

Loads ELISA data from plates 103 and 104, RNA data from PRP220307,
and normalizes protein to control means.

Returns: DataFrame with columns (day, rna, protein)
"""
function load_panel_a()
    # Load ELISA data
    elisa_103 = CSV.read("data/103_summary.tsv", DataFrame)
    elisa_103[!, :plate] .= 103
    elisa_104 = CSV.read("data/104_summary.tsv", DataFrame)
    elisa_104[!, :plate] .= 104
    elisa = vcat(elisa_103, elisa_104)

    # Load RNA/meta data
    meta = CSV.read("data/PRP220307.tsv", DataFrame)

    # Clean animal IDs (remove hyphens) - note: column name has space and capitals
    meta[!, :animal] = replace.(meta[!, Symbol("Subject ID")], "-" => "")
    meta[!, :day] = meta[!, Symbol("Necropsy Study Day (computed)")]

    # Classify treatment - note: column name has capitals and spaces
    meta[!, :tx] = ifelse.(meta[!, Symbol("Test Article (from group)")] .== "1087171", "ASO6", "none")

    # Normalize RNA (already in percent of PBS) - note exact capitalization
    meta[!, :rna] = meta[!, Symbol("mPrnp %PBS from hemisphere")] ./ 100

    # Select needed columns
    meta_clean = select(meta, :animal, :day, :tx, :rna)

    # Join ELISA with meta
    elisa_clean = select(elisa, :sample => :animal, :ngml_av, :plate)
    joined = innerjoin(elisa_clean, meta_clean, on=:animal)

    # Normalize protein to control mean within each plate
    grouped = groupby(joined, :plate)
    combined = combine(grouped) do df
        ctl_mean = mean(df.ngml_av[df.tx .== "none"])
        df[!, :protein] = df.ngml_av ./ ctl_mean
        return df
    end

    # Filter for ASO6 treated animals only
    aso6_data = filter(row -> row.tx == "ASO6", combined)

    # Return only needed columns
    return select(aso6_data, :day, :rna, :protein)
end

"""
Load and prepare Panel B data: ASON ki817 mice

Loads ELISA data from plates 175 and 177, RNA data from PRP230911,
and normalizes both protein and RNA to control means.

Returns: DataFrame with columns (day, rna, protein)
"""
function load_panel_b()
    # Load ELISA data for plates 175 and 177 only
    elisa_175 = CSV.read("data/175_summary.tsv", DataFrame)
    elisa_175[!, :plate] .= 175
    elisa_177 = CSV.read("data/177_summary.tsv", DataFrame)
    elisa_177[!, :plate] .= 177
    ason_prot = vcat(elisa_175, elisa_177)

    # Load RNA/PK data
    ason_pk_rna = CSV.read("data/PRP230911.tsv", DataFrame)
    ason_pk_rna[!, :animal] = replace.(ason_pk_rna[!, :animal], "-" => "")

    # Convert day column from "D1" format to numeric
    ason_pk_rna[!, :day_numeric] = parse.(Int, replace.(string.(ason_pk_rna[!, :day]), "D" => ""))

    # Join protein and RNA data
    ason_prot_clean = rename(ason_prot, :sample => :animal)
    joined = innerjoin(ason_prot_clean, ason_pk_rna, on=:animal)

    # Normalize RNA to PBS mean
    pbs_rna_mean = mean(joined.rna[joined.tx .== "PBS"])
    joined[!, :rna_norm] = joined.rna ./ pbs_rna_mean

    # Normalize protein to PBS mean within each plate
    grouped = groupby(joined, :plate)
    combined = combine(grouped) do df
        pbs_mean = mean(df.ngml_av[df.tx .== "PBS"])
        df[!, :protein] = df.ngml_av ./ pbs_mean
        return df
    end

    # Filter for non-PBS (treated) animals
    ason_data = filter(row -> row.tx != "PBS", combined)

    # Return only needed columns (use numeric day)
    return select(ason_data, :day_numeric => :day, :rna_norm => :rna, :protein)
end

"""
Analyze one ASO panel: fit both single and mixture models

Parameters:
- panel_name: Panel identifier (e.g., "ASO6_WT")
- days: Vector of time points
- rna: Vector of RNA measurements (fraction of control)
- protein: Vector of protein measurements (fraction of control)

Returns: NamedTuple with all fit results
"""
function analyze_aso_panel(panel_name, days, rna, protein)
    println("Analyzing panel: $panel_name")

    # Fit single-rate residual decay model
    thalf_single = fit_single_residual(days, rna, protein)

    # Calculate single model diagnostics
    aic_single = NaN
    bic_single = NaN

    if !isnan(thalf_single)
        # Create RNA interpolation (same as in fit function)
        rna_df = DataFrame(day=days, rna=rna)
        rna_agg = combine(groupby(rna_df, :day), :rna => mean => :rna)
        sort!(rna_agg, :day)

        rna_days = vcat([0.0], rna_agg.day)
        rna_vals = vcat([1.0], rna_agg.rna)
        R_interp = LinearInterpolation(rna_days, rna_vals, extrapolation_bc=Line())

        # Predict and calculate residuals
        lambda_single = log(2) / thalf_single
        pred_single = residual_decay_single(R_interp, lambda_single, days)
        resid_single = protein .- pred_single
        rss_single = sum(resid_single.^2)

        n = length(days)
        k_single = 2  # lambda, sigma
        aic_single = calculate_aic(rss_single, n, k_single)
        bic_single = calculate_bic(rss_single, n, k_single)
    end

    # Fit mixture model
    mix = fit_mixture_residual(days, rna, protein)

    # Calculate delta AIC/BIC
    delta_aic = mix.aic - aic_single
    delta_bic = mix.bic - bic_single

    # Calculate model ratio
    model_ratio = isnan(thalf_single) ? NaN : mix.effective_thalf / thalf_single

    return (
        panel = panel_name,
        n_obs = length(days),

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
        aic_mixture = mix.aic,
        bic_mixture = mix.bic,

        # Comparison
        delta_aic = delta_aic,
        delta_bic = delta_bic,
        model_ratio = model_ratio
    )
end

"""
Generate text report summarizing analysis results
"""
function generate_report(results)
    println("\n" * "="^80)
    println("ASO KNOCKDOWN MIXTURE MODEL ANALYSIS")
    println("="^80)

    for row in eachrow(results)
        println("\n$(row.panel):")
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
            if !isnan(row.thalf_effective)
                @printf("    Effective t½ = %.2f days\n", row.thalf_effective)
            else
                @printf("    Effective t½ = N/A (protein not decaying to steady state)\n")
            end
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

"""
Main analysis pipeline
"""
function main()
    println("="^80)
    println("FIG3 ASO KNOCKDOWN MIXTURE MODEL ANALYSIS")
    println("="^80)
    println()

    results = DataFrame()

    # Panel A: ASO6 WT mice
    println("Loading Panel A data (ASO6 WT mice)...")
    try
        panel_a_data = load_panel_a()
        @printf("  Loaded %d observations\n", nrow(panel_a_data))

        # Analyze panel
        result_a = analyze_aso_panel(
            "ASO6_WT",
            panel_a_data.day,
            panel_a_data.rna,
            panel_a_data.protein
        )

        push!(results, result_a)
    catch e
        println("  ERROR loading/analyzing Panel A: $e")
    end

    # Panel B: ASON ki817 mice
    println("\nLoading Panel B data (ASON ki817 mice)...")
    try
        panel_b_data = load_panel_b()
        @printf("  Loaded %d observations\n", nrow(panel_b_data))

        # Analyze panel
        result_b = analyze_aso_panel(
            "ASON_ki817",
            panel_b_data.day,
            panel_b_data.rna,
            panel_b_data.protein
        )

        push!(results, result_b)
    catch e
        println("  ERROR loading/analyzing Panel B: $e")
    end

    # Panel C: ASO6 RML mice (use only plate 226, plate 237 missing)
    println("\nLoading Panel C data (ASO6 RML prion-infected mice)...")
    try
        # Load ELISA data (only plate 226)
        elisa_226 = CSV.read("data/226_summary.tsv", DataFrame)
        elisa_226[!, :plate] .= 226

        # Load metadata
        meta = CSV.read("data/PRP240407.tsv", DataFrame)
        meta[!, :animal] = string.(meta[!, :animal])

        # Load RNA data
        rna_data = CSV.read("data/PRP240407_rna.tsv", DataFrame)
        rna_data[!, :animal] = string.(rna_data[!, :mouse])
        select!(rna_data, Not(:mouse))

        # Join and process
        elisa_clean = select(elisa_226, :sample => :animal, :ngml_av, :plate)
        joined = innerjoin(elisa_clean, meta, on=:animal)

        # Normalize protein to control mean
        grouped = groupby(joined, :plate)
        combined = combine(grouped) do df
            ctl_mean = mean(df.ngml_av[df.tx .== "none"])
            df[!, :protein] = df.ngml_av ./ ctl_mean
            return df
        end

        # Join with RNA
        final = innerjoin(combined, rna_data, on=:animal)

        # Filter for ASO6 treated animals
        aso6_data = filter(row -> row.tx == "ASO6", final)

        @printf("  Loaded %d observations\n", nrow(aso6_data))

        # Analyze panel
        result_c = analyze_aso_panel(
            "ASO6_RML",
            aso6_data.day,
            aso6_data.rna,
            aso6_data.protein
        )

        push!(results, result_c)
    catch e
        println("  ERROR loading/analyzing Panel C: $e")
    end

    # Panel D: Rats - NOTE: No RNA data available
    println("\nPanel D (ASO6 rats): SKIPPED (no RNA measurements available for residual decay model)")
    println("  Note: Residual decay mixture model requires RNA time series")
    println("  Rat data only has protein measurements, cannot fit Pt_mixture model")

    # Generate output directory
    println("\nCreating output directory...")
    mkpath("data/aso")

    # Save results
    println("Saving results...")
    CSV.write("data/aso/mixture_params.tsv", results, delim='\t')
    println("  Saved to: data/aso/mixture_params.tsv")

    # Generate report
    generate_report(results)

    println("Analysis complete!")
end

# Run main analysis
main()
