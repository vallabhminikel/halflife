module TurnoverModels

using Optim, LinearAlgebra, Interpolations, Statistics
using DataFrames, CSV, HypothesisTests, Printf, Distributions

export dt, free_lysine
export proportion_labeled_single, proportion_labeled_mixture
export residual_decay_single, residual_decay_mixture
export fit_single_isotopic, fit_mixture_isotopic
export fit_single_residual, fit_mixture_residual
export calculate_aic, calculate_bic, runs_test_residuals
export effective_halflife

# ============================================================================
# CONSTANTS
# ============================================================================

const dt = 0.01

# ============================================================================
# AVAILABILITY FUNCTIONS
# ============================================================================

"""
Free lysine availability function (Fornasiero 2018).
Models the proportion of newly synthesized protein that contains heavy isotope label.
"""
function free_lysine(t)
    return max(0.0, 1.0 - 0.503 * exp(-t * log(2) * 0.799) - 0.503 * exp(-t * log(2) / 39.423))
end

# ============================================================================
# ISOTOPIC LABELING MODELS
# ============================================================================

"""
Simulate proportion labeled: SINGLE RATE Model (Algebraic/Euler)

Parameters:
- thalf: Protein half-life (days)
- t_points: Time points to evaluate
- avails_func: Isotope availability function (default: free_lysine)

Returns: Proportion labeled at each time point
"""
function proportion_labeled_single(thalf, t_points, avails_func=free_lysine)
    lambda = log(2) / thalf
    t = collect(0:dt:maximum(t_points))
    proportion_heavy = zeros(length(t))

    for i in 2:lastindex(t)
        protein_turned_over = lambda * dt
        proportion_heavy[i] = (1 - protein_turned_over) * proportion_heavy[i-1] +
                              (protein_turned_over * avails_func(t[i]))
    end

    interp = LinearInterpolation(t, proportion_heavy, extrapolation_bc=Line())
    return [interp(tp) for tp in t_points]
end

"""
Simulate proportion labeled: MIXTURE Model (Algebraic/Euler)

Parameters:
- thalves: Vector of half-lives for each component (days)
- proportions: Vector of proportions for each component (sum to 1)
- t_points: Time points to evaluate
- avails_func: Isotope availability function (default: free_lysine)

Returns: Total proportion labeled at each time point
"""
function proportion_labeled_mixture(thalves, proportions, t_points, avails_func=free_lysine)
    t = collect(0:dt:maximum(t_points))
    K = length(thalves)
    prop_heavy_sub = zeros(length(t), K)

    for k in 1:K
        lambda_k = log(2) / thalves[k]
        for i in 2:lastindex(t)
            protein_turned_over = lambda_k * dt
            prop_heavy_sub[i, k] = (1 - protein_turned_over) * prop_heavy_sub[i-1, k] +
                                   (protein_turned_over * avails_func(t[i]))
        end
    end

    total_prop_heavy = prop_heavy_sub * proportions
    interp = LinearInterpolation(t, total_prop_heavy, extrapolation_bc=Line())
    return [interp(tp) for tp in t_points]
end

# ============================================================================
# RESIDUAL DECAY MODELS
# ============================================================================

"""
Simulate residual protein decay: SINGLE RATE Model

Parameters:
- R_interp: RNA interpolation object
- lambda: Decay rate constant (log(2)/t_half)
- t_points: Time points to evaluate

Returns: Residual protein levels at each time point

Differential equation: dP/dt = λ*R(t) - λ*P(t)
where P(t) is protein level, R(t) is RNA level, λ is decay rate
"""
function residual_decay_single(R_interp, lambda, t_points)
    t = collect(0:dt:maximum(t_points))
    P = ones(length(t))  # Start at 100% (normalized)

    for i in 2:lastindex(t)
        dt_val = t[i] - t[i-1]
        R = R_interp(t[i-1])
        dP = lambda * dt_val * R - lambda * dt_val * P[i-1]
        P[i] = P[i-1] + dP
    end

    interp = LinearInterpolation(t, P, extrapolation_bc=Line())
    return [interp(tp) for tp in t_points]
end

"""
Simulate residual protein decay: MIXTURE Model

Parameters:
- R_interp: RNA interpolation object
- lambdas: Vector of decay rate constants for each component
- weights: Vector of initial proportions for each component (sum to 1)
- t_points: Time points to evaluate

Returns: (total, components) tuple
- total: Total residual protein at each time point
- components: Matrix of component levels [time × components]

Models heterogeneous protein pools with different decay rates.
"""
function residual_decay_mixture(R_interp, lambdas, weights, t_points)
    t = collect(0:dt:maximum(t_points))
    K = length(lambdas)
    P_components = zeros(length(t), K)

    # Initialize each component with its weight
    for k in 1:K
        P_components[1, k] = weights[k]
    end

    for i in 2:lastindex(t)
        dt_val = t[i] - t[i-1]
        R = R_interp(t[i-1])

        for k in 1:K
            # dP = production - degradation
            # production: λ * R * weight (new synthesis proportional to RNA and component weight)
            # degradation: λ * P (degradation of existing protein in this component)
            dP = lambdas[k] * dt_val * R * weights[k] - lambdas[k] * dt_val * P_components[i-1, k]
            P_components[i, k] = P_components[i-1, k] + dP
        end
    end

    P_total = vec(sum(P_components, dims=2))

    interp_total = LinearInterpolation(t, P_total, extrapolation_bc=Line())
    return (
        total = [interp_total(tp) for tp in t_points],
        components = P_components
    )
end

# ============================================================================
# FITTING FUNCTIONS: ISOTOPIC LABELING
# ============================================================================

"""
Fit single-rate isotopic labeling model

Parameters:
- days: Observed time points
- prop: Observed proportion labeled
- avails_func: Isotope availability function (default: free_lysine)

Returns: Fitted half-life (days), or NaN if fit fails
"""
function fit_single_isotopic(days, prop, avails_func=free_lysine)
    length(days) < 2 && return NaN

    obj(thalf) = (thalf <= 0.01 || thalf > 100) ? Inf : sum((prop .- proportion_labeled_single(thalf, days, avails_func)).^2)

    try
        res = optimize(obj, 0.1, 100.0)
        return res.minimizer
    catch
        return NaN
    end
end

"""
Fit mixture isotopic labeling model (2 components)

Parameters:
- days: Observed time points
- prop: Observed proportion labeled
- avails_func: Isotope availability function (default: free_lysine)

Returns: NamedTuple with (thalves, proportions, effective_thalf)
"""
function fit_mixture_isotopic(days, prop, avails_func=free_lysine)
    length(days) < 4 && return (thalves=[NaN, NaN], proportions=[NaN, NaN], effective_thalf=NaN)

    function obj(p)
        t1, t2, pi1 = p
        (t1 <= 0.01 || t1 > 100 || t2 <= 0.01 || t2 > 100 || pi1 <= 0 || pi1 >= 1) && return Inf
        return sum((prop .- proportion_labeled_mixture([t1, t2], [pi1, 1-pi1], days, avails_func)).^2)
    end

    best_res = nothing
    best_val = Inf

    # Multiple starting points for robust optimization
    for t1 in [1.0, 5.0], t2 in [8.0, 25.0], pi in [0.3, 0.7]
        try
            res = optimize(obj, [t1, t2, pi], NelderMead())
            if res.minimum < best_val
                best_val = res.minimum
                best_res = res
            end
        catch; continue; end
    end

    if best_res === nothing
        return (thalves=[NaN, NaN], proportions=[NaN, NaN], effective_thalf=NaN)
    end

    t1, t2, pi1 = best_res.minimizer
    # Sort by half-life (fast first)
    if t1 > t2
        t1, t2 = t2, t1
        pi1 = 1 - pi1
    end
    pi2 = 1 - pi1

    eff_thalf = effective_halflife([t1, t2], [pi1, pi2])

    return (thalves=[t1, t2], proportions=[pi1, pi2], effective_thalf=eff_thalf)
end

# ============================================================================
# FITTING FUNCTIONS: RESIDUAL DECAY
# ============================================================================

"""
Fit single-rate residual decay model

Parameters:
- days: Observed time points
- rna: Observed RNA levels (fraction of control)
- protein: Observed protein levels (fraction of control)

Returns: Fitted half-life (days), or NaN if fit fails
"""
function fit_single_residual(days, rna, protein)
    length(days) < 2 && return NaN

    # Create RNA interpolation (prepend day 0 with RNA=1.0)
    # First, aggregate duplicates and sort
    rna_df = DataFrame(day=days, rna=rna)
    rna_agg = combine(groupby(rna_df, :day), :rna => mean => :rna)
    sort!(rna_agg, :day)

    rna_days = vcat([0.0], rna_agg.day)
    rna_vals = vcat([1.0], rna_agg.rna)
    R_interp = LinearInterpolation(rna_days, rna_vals, extrapolation_bc=Line())

    obj(lambda) = begin
        (lambda <= 0.001 || lambda > 10) && return Inf
        pred = residual_decay_single(R_interp, lambda, days)
        return sum((protein .- pred).^2)
    end

    try
        res = optimize(obj, 0.01, 5.0)
        lambda = res.minimizer
        return log(2) / lambda
    catch
        return NaN
    end
end

"""
Fit mixture residual decay model (2 components)

Parameters:
- days: Observed time points
- rna: Observed RNA levels (fraction of control)
- protein: Observed protein levels (fraction of control)

Returns: NamedTuple with (thalves, proportions, effective_thalf, aic, bic)
"""
function fit_mixture_residual(days, rna, protein)
    length(days) < 4 && return (thalves=[NaN, NaN], proportions=[NaN, NaN],
                                 effective_thalf=NaN, aic=NaN, bic=NaN)

    # Create RNA interpolation
    # First, aggregate duplicates and sort
    rna_df = DataFrame(day=days, rna=rna)
    rna_agg = combine(groupby(rna_df, :day), :rna => mean => :rna)
    sort!(rna_agg, :day)

    rna_days = vcat([0.0], rna_agg.day)
    rna_vals = vcat([1.0], rna_agg.rna)
    R_interp = LinearInterpolation(rna_days, rna_vals, extrapolation_bc=Line())

    function obj(p)
        lambda1, lambda2, w1 = p
        (lambda1 <= 0.001 || lambda1 > 10 || lambda2 <= 0.001 || lambda2 > 10 ||
         w1 <= 0 || w1 >= 1) && return Inf

        pred = residual_decay_mixture(R_interp, [lambda1, lambda2], [w1, 1-w1], days)
        return sum((protein .- pred.total).^2)
    end

    best_res = nothing
    best_val = Inf

    # Multiple starting points for robust fitting
    for l1 in [0.1, 0.5], l2 in [0.05, 0.2], w in [0.3, 0.7]
        try
            res = optimize(obj, [l1, l2, w], NelderMead())
            if res.minimum < best_val
                best_val = res.minimum
                best_res = res
            end
        catch; continue; end
    end

    if best_res === nothing
        return (thalves=[NaN, NaN], proportions=[NaN, NaN],
                effective_thalf=NaN, aic=NaN, bic=NaN)
    end

    lambda1, lambda2, w1 = best_res.minimizer

    # Sort by half-life (fast first)
    t1 = log(2) / lambda1
    t2 = log(2) / lambda2
    if t1 > t2
        t1, t2 = t2, t1
        w1 = 1 - w1
    end
    w2 = 1 - w1

    # Use residual-specific effective halflife calculation
    eff_thalf = effective_halflife([t1, t2], [w1, w2])

    # Calculate AIC/BIC
    n = length(days)
    pred = residual_decay_mixture(R_interp, [log(2)/t1, log(2)/t2], [w1, w2], days)
    rss = sum((protein .- pred.total).^2)
    k_mixture = 4  # lambda1, lambda2, w1, sigma
    aic_val = calculate_aic(rss, n, k_mixture)
    bic_val = calculate_bic(rss, n, k_mixture)

    return (thalves=[t1, t2], proportions=[w1, w2], effective_thalf=eff_thalf,
            aic=aic_val, bic=bic_val)
end

# ============================================================================
# STATISTICAL UTILITIES
# ============================================================================

"""
Calculate Akaike Information Criterion

Parameters:
- rss: Residual sum of squares
- n: Number of observations
- k: Number of parameters

Returns: AIC value (lower is better)
"""
function calculate_aic(rss, n, k)
    return n * log(rss/n) + 2 * k
end

"""
Calculate Bayesian Information Criterion

Parameters:
- rss: Residual sum of squares
- n: Number of observations
- k: Number of parameters

Returns: BIC value (lower is better)
"""
function calculate_bic(rss, n, k)
    return n * log(rss/n) + k * log(n)
end


"""
Calculate effective half-life for simple exponential mixture model

For simple exponential decay (isotopic labeling), the effective half-life
is the time at which the mixture reaches 50%: sum(w_i * 0.5^(t/t½_i)) = 0.5

Parameters:
- thalves: Vector of half-lives for each component
- proportions: Vector of proportions for each component

Returns: Effective half-life (days)
"""
function effective_halflife(thalves, proportions)
    # Mixture remaining at time t: sum(prop[i] * 0.5^(t/thalf[i]))
    mix_rem(t) = sum(proportions[i] * (0.5^(t/thalves[i])) for i in eachindex(thalves))
    eff_obj(t) = abs(mix_rem(t) - 0.5)

    try
        return optimize(eff_obj, 0.0, maximum(thalves)*2).minimizer
    catch
        return NaN
    end
end

end # module
