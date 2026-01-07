include("TurnoverModels.jl")
using .TurnoverModels
using DataFrames, CSV, HypothesisTests, Printf, Distributions, Colors

# ============================================================================
# CORE ANALYSIS PIPELINE
# ============================================================================

"""
Main Analysis Function:
Iterates through the data once, fits both models, computes statistics/residuals.
Returns a master DataFrame with all results.
"""
function analyze_turnover_models(df)
    println("Running core model analysis (Single & Mixture fits)...")

    results = DataFrame()

    grouped = groupby(df, [:protein, :peptide, :age, :genotype])

    for subdf in grouped
        prot, pep, age, geno = subdf.protein[1], subdf.peptide[1], subdf.age[1], subdf.genotype[1]
        days = subdf.chow_days
        props = subdf.prop_labeled

        t_single = fit_single_isotopic(days, props)

        mix = fit_mixture_isotopic(days, props)

        runs_stat = (stat=NaN, p=NaN, n_runs=0, n_pos=0, n_neg=0)
        if !isnan(t_single)
            preds = proportion_labeled_single(t_single, days)
            resids = props .- preds
            log_resids = log.(abs.(resids) .+ 1e-10)

            try
                rt = RunsTest(log_resids)
                runs_stat = (stat=rt.z, p=pvalue(rt), n_runs=rt.n_runs, n_pos=rt.n1, n_neg=rt.n0)
            catch
            end
        end

        push!(results, (
            protein = prot, peptide = pep, age = age, genotype = geno,
            n_obs = length(days),

            # Single Results
            thalf_single = t_single,

            # Mixture Results
            thalf_fast      = mix.thalves[1],
            thalf_slow      = mix.thalves[2],
            prop_fast       = mix.proportions[1],
            prop_slow       = mix.proportions[2],
            thalf_effective = mix.effective_thalf,

            # Comparison
            model_ratio     = isnan(t_single) ? NaN : mix.effective_thalf / t_single,

            # Residual Stats
            runs_pvalue     = runs_stat.p,
            runs_n          = runs_stat.n_runs
        ))
    end
    return results
end

"""
Compute Genotypic Differences (FFI vs Control) based on the Master Results
"""
function compute_genotypic_differences(raw_df, results_df, leg_df)
    println("Computing genotypic differences...")

    full_res = innerjoin(results_df, leg_df, on=:genotype)

    diffs = combine(groupby(full_res, [:protein, :peptide, :age])) do df
        is_ffi = df.genotype .== "129(TT-3F4-FFI)HOZ"

        if sum(is_ffi) == 0 || sum(.!is_ffi) == 0
            return DataFrame()
        end

        thalf_ffi = df.thalf_effective[is_ffi][1]
        thalf_ctrl = mean(df.thalf_effective[.!is_ffi]) # Average if multiple controls

        raw_sub = filter(r -> r.peptide == df.peptide[1] && r.age == df.age[1], raw_df)
        raw_ffi = raw_sub.total[raw_sub.genotype .== "129(TT-3F4-FFI)HOZ"]
        raw_ctrl = raw_sub.total[raw_sub.genotype .!= "129(TT-3F4-FFI)HOZ"]

        pval_total = 1.0
        try
            pval_total = pvalue(UnequalVarianceTTest(raw_ffi, raw_ctrl))
        catch; end

        prop_ffi = raw_sub.prop_labeled[raw_sub.genotype .== "129(TT-3F4-FFI)HOZ"]
        prop_ctrl = raw_sub.prop_labeled[raw_sub.genotype .!= "129(TT-3F4-FFI)HOZ"]
        pval_label = 1.0
        try
            pval_label = pvalue(UnequalVarianceTTest(prop_ffi, prop_ctrl))
        catch; end

        ratio_total = mean(raw_ffi) / mean(raw_ctrl)
        ratio_thalf = thalf_ffi / thalf_ctrl

        DataFrame(
            pval_total = pval_total,
            ratio_total = ratio_total,
            ratio_thalf = ratio_thalf,
            pval_label_ffi = pval_label,
            log2_ratio_total = log2(ratio_total),
            log2_ratio_thalf = log2(ratio_thalf),
            pepnick = first(df.peptide[1], 4)
        )
    end

    pep_counts = combine(groupby(diffs, [:protein, :peptide]), :age => (x -> length(unique(x))) => :n_ages)
    diffs = innerjoin(diffs, pep_counts, on=[:protein, :peptide])
    filter!(r -> r.n_ages == 2, diffs)
    select!(diffs, Not(:n_ages))

    return diffs
end

# ============================================================================
# REPORTING & OUTPUTS
# ============================================================================

function generate_text_report(results_df, leg_df)
    println("\n" * "="^80)
    println("MIXTURE MODEL ANALYSIS REPORT")
    println("="^80)

    # Add friendly names
    df = innerjoin(results_df, leg_df, on=:genotype)

    for subdf in groupby(df, [:age, :peptide])
        println("\nAGE: $(uppercase(subdf.age[1])) | PEPTIDE: $(subdf.peptide[1])")
        println("-"^60)

        for row in eachrow(subdf)
            println("  $(row.disp):")
            @printf("    n = %d obs\n", row.n_obs)
            if !isnan(row.thalf_single)
                @printf("    Single Model: t½ = %.2f days\n", row.thalf_single)
            else
                println("    Single Model: Failed")
            end

            if !isnan(row.thalf_effective)
                println("    Mixture Model:")
                @printf("      Fast: %.2f d (%.1f%%) | Slow: %.2f d (%.1f%%)\n",
                        row.thalf_fast, row.prop_fast*100, row.thalf_slow, row.prop_slow*100)
                @printf("      Effective t½ = %.2f days (Ratio Mix/Single = %.2f)\n",
                        row.thalf_effective, row.model_ratio)
            else
                println("    Mixture Model: Failed")
            end
        end
    end
    println("\n" * "="^80 * "\n")
end

# ============================================================================
# MAIN EXECUTION
# ============================================================================

function main()
    println("Starting Analysis...")

    raw_data = CSV.read("data/iqp/ffi_turnover.tsv", DataFrame)
    raw_data.total = raw_data.light .+ raw_data.heavy
    raw_data.prop_labeled = raw_data.heavy ./ raw_data.total

    leg = DataFrame(
        genotype = ["129(TT-3F4-FFI)HOZ", "129(TT-3F4WT)", "B6/N"],
        disp = ["ki-3F4-FFI", "ki-3F4-WT", "C57BL/6N WT"],
        color = ["#D95F02", "#22127A", "#77127A"]
    )

    master_results = analyze_turnover_models(raw_data)

    genotypic_diffs = compute_genotypic_differences(raw_data, master_results, leg)

    generate_text_report(master_results, leg)

    out_table = select(innerjoin(master_results, leg, on=:genotype),
                      :protein, :peptide, :age, :genotype, :disp => :genotype_label,
                      :n_obs, :thalf_fast, :prop_fast, :thalf_slow, :prop_slow,
                      :thalf_effective, :runs_pvalue, :runs_n)
    CSV.write("display_items/table-mixture-model-results.tsv", out_table, delim='\t')

    CSV.write("data/iqp/mixture_params.tsv", master_results, delim='\t')

    CSV.write("data/iqp/genotypic_diffs.tsv", genotypic_diffs, delim='\t')

    println("All files saved successfully.")
end

# Run
main()
