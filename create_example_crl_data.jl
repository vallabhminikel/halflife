"""
Create Example CRL-style Labeling Data

This generates synthetic data in the format expected by fit_crl_labeling_data.jl
to demonstrate the analysis pipeline.
"""

using CSV, DataFrames, Random
include("mixture_decay_optim.jl")

function generate_crl_example_data(; save_path="example_crl_data.csv")
    Random.seed!(123)

    # Experimental design
    days = [0, 1, 2, 3, 4, 6, 8]
    tissues = ["brain", "liver", "spleen"]
    n_replicates = 4

    # True parameters by tissue (different heterogeneity)
    tissue_params = Dict(
        "brain" => (λ=log(2)/8, ψ=[log(2)/4, log(2)/15], π=[0.3, 0.7]),   # Mostly slow
        "liver" => (λ=log(2)/5, ψ=[log(2)/2, log(2)/10], π=[0.7, 0.3]),   # Mostly fast
        "spleen" => (λ=log(2)/6, ψ=[log(2)/3, log(2)/12], π=[0.5, 0.5])   # Balanced
    )

    # Generate data
    data = DataFrame(
        day = Float64[],
        proportion_heavy = Float64[],
        tissue = String[],
        replicate = Int[]
    )

    for tissue in tissues
        params = tissue_params[tissue]

        # Fine time grid for simulation
        t_fine = range(0, 8, length=1000)

        # RNA constant (no knockdown)
        R_fine = ones(length(t_fine))
        R_interp = LinearInterpolation(t_fine, R_fine, extrapolation_bc=Line())

        # Lysine availability
        avails_fine = free_lysine.(t_fine)
        avails_interp = LinearInterpolation(t_fine, avails_fine, extrapolation_bc=Line())

        # Generate true trajectory
        P_obs_true, _, _, _ = solve_mixture_dynamics_with_labeling(
            R_interp, avails_interp, t_fine, params.λ, params.ψ, params.π
        )

        # Sample at experimental timepoints
        P_interp = LinearInterpolation(t_fine, P_obs_true, extrapolation_bc=Line())

        for day in days
            P_expected = P_interp(day)

            # Add biological + technical noise
            noise_level = day == 0 ? 0.005 : 0.02  # Less noise at baseline

            for rep in 1:n_replicates
                P_observed = P_expected + noise_level * randn()
                P_observed = max(0.0, min(1.0, P_observed))  # Bound [0,1]

                push!(data, (day, P_observed, tissue, rep))
            end
        end
    end

    # Save
    CSV.write(save_path, data)
    println("Generated example CRL data: $save_path")
    println("  $(nrow(data)) total observations")
    println("  $(length(tissues)) tissues: $(join(tissues, ", "))")
    println("  $(length(days)) timepoints: $(join(days, ", "))")
    println("  $n_replicates replicates per condition")

    # Show summary
    println("\nTrue parameters:")
    for tissue in tissues
        params = tissue_params[tissue]
        println("\n$tissue:")
        println("  λ = $(round(params.λ, digits=4)) (t₁/₂ = $(round(log(2)/params.λ, digits=1))d)")
        println("  ψ₁ = $(round(params.ψ[1], digits=4)) (t₁/₂ = $(round(log(2)/params.ψ[1], digits=1))d), π₁ = $(params.π[1])")
        println("  ψ₂ = $(round(params.ψ[2], digits=4)) (t₁/₂ = $(round(log(2)/params.ψ[2], digits=1))d), π₂ = $(params.π[2])")
    end

    return data
end

# Generate the example data
println("="^70)
println("Creating Example CRL Labeling Data")
println("="^70)

data = generate_crl_example_data(save_path="example_crl_data.csv")

println("\n" * "="^70)
println("Next steps:")
println("="^70)
println("\n1. View the data:")
println("   >>> using CSV, DataFrames")
println("   >>> df = CSV.read(\"example_crl_data.csv\", DataFrame)")
println("   >>> first(df, 10)")

println("\n2. Run the analysis:")
println("   >>> include(\"fit_crl_labeling_data.jl\")")
println("   >>> df_summary, results = main_analysis(\"example_crl_data.csv\")")

println("\n3. This will:")
println("   - Fit single-rate and mixture models to each tissue")
println("   - Compare models using AIC")
println("   - Generate plots for each tissue")
println("   - Show which tissues exhibit heterogeneity")
