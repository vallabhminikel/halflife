"""
Generate comparison figure showing label accumulation curves for both models
Side-by-side comparison: Single (left) vs Mixture (right)
"""

using CairoMakie, LaTeXStrings, DataFrames, CSV, Colors, Printf
using Statistics: mean,std
using Distributions: TDist, quantile
include("TurnoverModels.jl")
using .TurnoverModels

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

    use_peptides = ["GENFTETDVK", "VVEQMCVTQYQK"]
    ages = ["young", "aged"]

    fig = Figure(size=(1400, 1000), fontsize=11)
    
    # Create 2×2 GridLayouts for each genotype (2 ages × 2 peptides)
    gl = [GridLayout() for _ in 1:3]
    
    #Label(fig[1, :, Top()], "Young", fontsize=13, font=:bold, padding=(0, 0, 20, 0))
    
    fig[1, 1] = gl[1]  # geno 1
    fig[1, 2] = gl[2]  # geno 2
    fig[1, 3] = gl[3]  # geno 3
    
    #Label(fig[2, :, Top()], "Aged", fontsize=13, font=:bold, padding=(0, 0, 20, 0))
    
    gl_aged = [GridLayout() for _ in 1:3]
    fig[2, 1] = gl_aged[1]  # geno 1
    fig[2, 2] = gl_aged[2]  # geno 2
    fig[2, 3] = gl_aged[3]  # geno 3

    # Panel labels (A, B, C for the three genotypes)
    avails_func = free_lysine

    # Loop through genotypes (3 columns)
    for (geno_idx, geno_row) in enumerate(eachrow(leg))
        # Add panel label for this genotype
        panel_letter = string(Char('A' + geno_idx - 1))
        
        # Create 2×2 grid: rows=ages, cols=peptides
        for (age_idx, this_age) in enumerate(ages)
            # Select appropriate grid layout
            current_gl = age_idx == 1 ? gl[geno_idx] : gl_aged[geno_idx]
            row_in_gl = 1  # Always row 1 within each GridLayout
            
            for (pep_idx, this_peptide) in enumerate(use_peptides)
                
                # Create axis
                ax = Axis(current_gl[row_in_gl, pep_idx],
                    xlabel= "day",
                    ylabel=pep_idx == 1 ? "proportion labeled" : "",
                    title=age_idx == 1 ? "Young " * this_peptide[1:4] : "Aged " * this_peptide[1:4],
                    titlesize=11,
                    width = 125,
                    height = 125,
                    yticks=0:0.1:0.5
                )
                hidespines!(ax, :t, :r)
                limits!(ax, 0, 8.5, 0, 0.55)

                # Add panel label to top-left of first subplot in each genotype column (young only)
                if age_idx == 1 && pep_idx == 1
                    Label(current_gl[row_in_gl, pep_idx, TopLeft()], 
                          panel_letter * ")", 
                          fontsize=14, font=:bold, padding=(0, 5, 5, 0))
                end

                t_fine = collect(0:0.1:8)
                t_fine = collect(0:0.1:8)

                curve_data = filter(row -> row.genotype == geno_row.genotype &&
                                          row.peptide == this_peptide &&
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

                # SINGLE MODEL FIT (dashed)
                if nrow(curve_data) > 0
                    thalf_single = fit_single_isotopic(
                        curve_data.chow_days,
                        curve_data.prop_labeled,
                        avails_func
                    )

                    if !isnan(thalf_single)
                        pred_single = proportion_labeled_single(thalf_single, t_fine, avails_func)
                        lines!(ax, t_fine, pred_single,
                              color=geno_row.color, linewidth=2, linestyle=:dash)
                    end
                end

                # MIXTURE MODEL FIT (solid)
                if nrow(curve_data) > 0
                    mixture_fit = fit_mixture_isotopic(
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
                        lines!(ax, t_fine, pred_mixture,
                              color=geno_row.color, linewidth=2, linestyle=:solid)
                    end
                end

                # Plot observed data with error bars
                scatter!(ax, summary_df.chow_days, summary_df.mean_prop,
                         color=geno_row.color, markersize=8,
                         marker=:circle, strokewidth=2, strokecolor=geno_row.color)
                errorbars!(ax, summary_df.chow_days, summary_df.mean_prop, 
                    summary_df.mean_prop .- summary_df.l95, 
                    summary_df.u95 .- summary_df.mean_prop,
                    color=geno_row.color, linewidth=2,
                    whiskerwidth=10)
            end  # end peptide loop
        end  # end age loop
    end  # end genotype loop

    # Create legend with both genotype colors and model types
    legend_elements = [
        [LineElement(color=leg[i, :color], linewidth=2) for i in 1:3]...,
        LineElement(color=:black, linewidth=2, linestyle=:dash),
        LineElement(color=:black, linewidth=2, linestyle=:solid)
    ]
    legend_labels = [
        leg[1, :disp],
        leg[2, :disp],
        leg[3, :disp],
        "single-rate fit",
        "mixture fit"
    ]

    Legend(fig[3, :],
          legend_elements, legend_labels,
          orientation=:horizontal,
          framevisible=false,
          halign=:center,
          nbanks=1)

    filename = "display_items/figure-s8.png"
    resize_to_layout!(fig)

    save(filename, fig, px_per_unit=2)
    println("done.\nModel comparison figure saved to $(filename)")
    return fig
end

# Generate model comparison figure
generate_figure_model_comparison()
