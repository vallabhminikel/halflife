"""
Demo: Lysine Labeling Dynamics in Mixture Decay Model

This demonstrates how lysine availability affects observed protein levels.
"""

include("mixture_decay_optim.jl")

function demo_labeling()
    println("="^70)
    println("Lysine Labeling Demo")
    println("="^70)

    # Setup
    t_data = collect(0.0:0.5:100.0)
    R_data = ones(length(t_data))  # Constant RNA (no knockdown)

    # Parameters
    λ = log(2) / 5  # Translation rate
    ψ_vec = [log(2)/1, log(2)/12]  # Fast and slow degradation
    π_vec = [0.5, 0.5]

    println("\n1. System parameters:")
    println("  λ (translation) = $(round(λ, digits=4))")
    println("  ψ₁ (fast deg)   = $(round(ψ_vec[1], digits=4)), t₁/₂ = 3 days")
    println("  ψ₂ (slow deg)   = $(round(ψ_vec[2], digits=4)), t₁/₂ = 12 days")

    # Create interpolation functions
    t_fine = range(0, maximum(t_data), length=1000)
    R_interp = LinearInterpolation(t_data, R_data, extrapolation_bc=Line())

    # Lysine availability over time
    avails_fine = free_lysine.(t_fine)
    avails_interp = LinearInterpolation(t_fine, avails_fine, extrapolation_bc=Line())

    # Solve WITHOUT labeling (standard model)
    P_total_no_label, P_sub_no_label = solve_mixture_dynamics(R_interp, t_fine, λ, ψ_vec, π_vec)

    # Solve WITH labeling
    P_obs_with_label, P_obs_sub, P_total_with_label, L_sub =
        solve_mixture_dynamics_with_labeling(R_interp, avails_interp, t_fine, λ, ψ_vec, π_vec)

    # Create visualization
    fig = Figure(size=(1200, 800), fontsize=14)

    # Panel A: Lysine availability
    ax1 = Axis(fig[1, 1],
        xlabel="Time (days)",
        ylabel="Lysine Availability",
        title="A. Free Lysine Dynamics")

    lines!(ax1, t_fine, avails_fine, color=:green, linewidth=3)
    hlines!(ax1, [1.0], color=:gray, linestyle=:dash, label="Maximum")

    println("\n2. Lysine dynamics:")
    println("  Initial availability: 0.0 (no labeled lysine)")
    println("  Final availability: $(round(avails_fine[end], digits=3))")
    println("  → Protein can only be labeled when lysine is available!")

    # Panel B: Total protein (without labeling consideration)
    ax2 = Axis(fig[1, 2],
        xlabel="Time (days)",
        ylabel="Total Protein",
        title="B. Total Protein (Unlabeled + Labeled)")

    lines!(ax2, t_fine, P_total_no_label, color=:blue, linewidth=3, label="Total")
    lines!(ax2, t_fine, P_sub_no_label[:, 1], color=:purple, linestyle=:dash, label="Fast component")
    lines!(ax2, t_fine, P_sub_no_label[:, 2], color=:orange, linestyle=:dash, label="Slow component")

    axislegend(ax2, position=:rb)

    # Panel C: Labeled fraction
    ax3 = Axis(fig[2, 1],
        xlabel="Time (days)",
        ylabel="Labeled Fraction",
        title="C. Fraction of Protein That Is Labeled")

    lines!(ax3, t_fine, L_sub[:, 1], color=:purple, linewidth=3, label="Fast component")
    lines!(ax3, t_fine, L_sub[:, 2], color=:orange, linewidth=3, label="Slow component")
    hlines!(ax3, [avails_fine[end]], color=:green, linestyle=:dot, label="Equilibrium")

    axislegend(ax3, position=:rb)

    println("\n3. Labeled fraction at equilibrium:")
    println("  Fast component: $(round(L_sub[end, 1], digits=3))")
    println("  Slow component: $(round(L_sub[end, 2], digits=3))")
    println("  Lysine availability: $(round(avails_fine[end], digits=3))")
    println("  → At equilibrium, labeled fraction = lysine availability!")

    # Panel D: Observed protein (only labeled is detectable)
    ax4 = Axis(fig[2, 2],
        xlabel="Time (days)",
        ylabel="Observed Protein",
        title="D. Observed Protein (Only Labeled Is Detectable)")

    lines!(ax4, t_fine, P_obs_with_label, color=:blue, linewidth=3, label="Total observed")
    lines!(ax4, t_fine, P_obs_sub[:, 1], color=:purple, linestyle=:dash, label="Fast (observed)")
    lines!(ax4, t_fine, P_obs_sub[:, 2], color=:orange, linestyle=:dash, label="Slow (observed)")

    # Show asymptotic level
    asymptote = P_total_no_label[end] * avails_fine[end]
    hlines!(ax4, [asymptote], color=:red, linestyle=:dot, linewidth=2, label="Asymptote")

    axislegend(ax4, position=:rb)

    println("\n4. Key insight:")
    println("  WITHOUT labeling: Protein → $(round(P_total_no_label[end], digits=3)) (total protein)")
    println("  WITH labeling:    Protein → $(round(P_obs_with_label[end], digits=3)) (observed)")
    println("  Ratio: $(round(P_obs_with_label[end] / P_total_no_label[end], digits=3)) = $(round(avails_fine[end], digits=3))")
    println("\n  Observed protein = Total protein × Lysine availability")
    println("  → Protein doesn't go to zero, it goes to a fraction!")

    save("labeling_demo.png", fig)
    println("\n5. Saved: labeling_demo.png")

    println("\n" * "="^70)
    println("SUMMARY")
    println("="^70)
    println("\nIn isotope labeling experiments:")
    println("1. Initially, no labeled lysine → no labeled protein (all invisible)")
    println("2. As labeled lysine becomes available, new protein incorporates it")
    println("3. Old protein slowly turns over and is replaced by labeled protein")
    println("4. At equilibrium: labeled fraction = lysine availability")
    println("5. Observed protein = Total × Labeled_fraction")
    println("\nFast-degrading proteins reach labeling equilibrium faster!")
    println("Slow-degrading proteins take longer to fully label.")

    return fig
end

demo_labeling()
