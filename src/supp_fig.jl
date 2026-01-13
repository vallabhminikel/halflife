"""
Generate comparison figure showing label accumulation curves for both models
Side-by-side comparison: Single (left) vs Mixture (right)
"""

using CairoMakie, LaTeXStrings, DataFrames, CSV, Colors, Printf
using Statistics: mean
include("TurnoverModels.jl")
using .TurnoverModels


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

    # Panel labels (A–D)
    panel_letters = ['A', 'B', 'C', 'D']
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
                height = 200,
                )
            hidespines!(ax_single, :t, :r)
            ax_mixture = Axis(gl_mixture[age_idx, pep_idx],
                xlabel=age_idx == 2 ? "day" : "",
                ylabel=pep_idx == 1 ? "proportion labeled" : "",
                title="$(this_age) $(this_peptide[1:4])",
                width = 200,
                height = 200
                )
            hidespines!(ax_mixture, :t, :r)

            xlims!(ax_single, -0.5, 8.5)
            xlims!(ax_mixture, -0.5, 8.5)
            ylims!(ax_single, 0, 0.5)
            ylims!(ax_mixture, 0, 0.5)

            # Add top-left panel labels (A–D on left, E–H on right)
            n_panels = length(ages) * length(use_peptides)
            left_letter = string(Char('A' + (panel_idx - 1))) * ")"
            right_letter = string(Char('A' + (panel_idx - 1 + n_panels))) * ")"
            Label(gl_single[age_idx, pep_idx, TopLeft()], left_letter, fontsize=14, font=:bold, padding=(0,10,10,0))
            Label(gl_mixture[age_idx, pep_idx, TopLeft()], right_letter, fontsize=14, font=:bold, padding=(0,10,10,0))

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
                    )
                end
                sort!(summary_df, :chow_days)

                # SINGLE MODEL FIT
                if nrow(curve_data) > 0
                    thalf_single = fit_single_isotopic(
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
                        lines!(ax_mixture, t_fine, pred_mixture,
                              color=geno_row.color, linewidth=2, linestyle=:dash)
                    end
                end

                # Plot observed data on both panels
                for ax in [ax_single, ax_mixture]
                    scatter!(ax, curve_data.chow_days, curve_data.prop_labeled,
                            color=:transparent, markersize=8,
                            marker=:circle, strokewidth=2, strokecolor=geno_row.color)
                end
            end
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

    filename = "display_items/figure-s8.png"
    resize_to_layout!(fig)

    save(filename, fig, px_per_unit=2)
    println("done.\nModel comparison figure saved to $(filename)")
    return fig
end

# Generate model comparison figure
generate_figure_model_comparison()
