# Mixture Decay Model with Lysine Labeling

This Julia implementation fits heterogeneous protein decay models to isotope labeling data, accounting for lysine availability dynamics.

## Overview

The model simultaneously accounts for:
1. **Heterogeneous decay rates**: Different protein subpopulations with distinct degradation rates
2. **Lysine labeling dynamics**: Only labeled protein is detectable
3. **RNA-protein coupling**: Protein production depends on RNA levels

## Model Equations

### Basic Mixture Model
```
dP_k/dt = λ * R(t) - ψ_k * P_k(t)
```

Where:
- λ = single translation rate (same for all components)
- ψ_k = degradation rate for component k (different for each)
- R(t) = RNA level over time

### With Lysine Labeling
```
dP_k/dt = λ * R(t) - ψ_k * P_k(t)              [Total protein]
dL_k/dt = avails(t) * ψ_k * (1-L_k) - (1-avails(t)) * ψ_k * L_k   [Labeled fraction]

Observed = P_k(t) * L_k(t)                      [Only labeled is detectable]
```

**Key insight**: At equilibrium, observed protein → avails(∞) × Total protein, **NOT zero**!

## Files

### Core Implementation
- **`mixture_decay_optim.jl`**: Main implementation with all solvers and fitting functions
- **`mixture_decay_diffeq.jl`**: Alternative implementation using DifferentialEquations.jl

### Demonstrations
- **`demo_labeling.jl`**: Visualize lysine labeling dynamics
- **`create_example_crl_data.jl`**: Generate synthetic CRL-style data

### Data Analysis
- **`fit_crl_labeling_data.jl`**: Fit mixture models to actual proteomics data

## Quick Start

### 1. Generate Example Data
```julia
julia> include("create_example_crl_data.jl")
```

This creates `example_crl_data.csv` with synthetic data for 3 tissues.

### 2. Fit Models
```julia
julia> include("fit_crl_labeling_data.jl")
julia> df_summary, results = main_analysis("example_crl_data.csv")
```

This will:
- Fit single-rate and mixture models to each tissue
- Compare models using AIC
- Generate plots showing fits and subpopulation dynamics
- Identify which tissues show heterogeneity

### 3. Visualize Labeling Dynamics
```julia
julia> include("demo_labeling.jl")
```

Shows how lysine availability affects observed protein levels.

## Model Features

### Standard Mixture Model
- Fits λ (translation), ψ₁, ψ₂, π₁ (4 parameters)
- Allows any relationship between rates

### Constrained Model (∑ψₖ = λ)
- Enforces: ψ₁ + ψ₂ + ... + ψ_K = λ
- Reduces to 3 parameters
- Ensures total degradation equals production
- Better parameter identifiability

## Model Selection

Models are compared using:
- **RSS**: Sum of squared residuals
- **AIC**: Akaike Information Criterion
- **Residual patterns**: Systematic bias indicates heterogeneity

**Rule of thumb**:
- ΔAIC > 4: Strong evidence for mixture
- ΔAIC > 2: Moderate evidence for mixture
- ΔAIC < 2: Weak evidence

## Equilibrium Properties

For component k:
```
P_k,ss = (λ / ψ_k) × R_ss           [Total protein at equilibrium]
L_k,ss = avails(∞)                   [Labeled fraction at equilibrium]
Observed = (λ / ψ_k) × R_ss × avails(∞)
```

**Amplification factor**: λ / ψ_k
- If > 1: Component accumulates above RNA level
- If < 1: Component depleted below RNA level
- If = 1: Component tracks RNA level

**Labeling equilibrium**:
- Fast-degrading components label quickly (high turnover)
- Slow-degrading components label slowly (low turnover)
- All components eventually reach L_k,ss = avails(∞)

## Data Format

CSV file with columns:
- `day`: Time point (days)
- `proportion_heavy`: Fraction of protein that's heavy-labeled [0,1]
- `tissue`: Tissue identifier (e.g., "brain", "liver")
- `replicate` (optional): Biological replicate ID

Example:
```csv
day,proportion_heavy,tissue,replicate
0,0.01,brain,1
0,0.02,brain,2
1,0.15,brain,1
1,0.14,brain,2
...
```

## Interpreting Results

### Heterogeneity Indicators
1. **ΔAIC > 4**: Strong statistical evidence
2. **Different half-lives**: ψ₁ ≠ ψ₂ (> 2-fold difference is meaningful)
3. **Significant proportions**: Both π₁ and π₂ > 0.1
4. **Residual pattern**: Single-rate shows systematic bias

### Biological Interpretation

**Fast component** (short t₁/₂):
- Cytosolic, soluble protein
- Normal cellular turnover
- Rapidly labeled (high turnover)

**Slow component** (long t₁/₂):
- Membrane-bound or aggregated
- Resistant to degradation
- Slowly labeled (low turnover)
- May accumulate to higher levels (if λ/ψ > 1)

## Advanced Usage

### Fitting with Constraint
```julia
λ, ψ_vec, π_vec, result = fit_mixture_model_with_labeling(
    t_data, R_data, P_data, 2,
    constrained=true  # Enforce ∑ψₖ = λ
)
```

### More Components
```julia
# Fit 3-component mixture
λ, ψ_vec, π_vec, result = fit_mixture_model_with_labeling(
    t_data, R_data, P_data, 3
)
```

### Custom Initial Guesses
```julia
λ, ψ_vec, π_vec, result = fit_mixture_model_with_labeling(
    t_data, R_data, P_data, 2,
    λ_init = 0.15,
    ψ_init = [0.3, 0.05],  # Fast and slow
    π_init = [0.8]         # Proportion in first component
)
```

## Comparison with R Implementation

The Julia code directly implements the R functions:
- `free_lysine(t)` → lysine availability (Fornasiero 2018)
- `proportion_labeled()` → labeled fraction dynamics
- `Pt()` → protein dynamics with RNA coupling
- Mixture extension with separate ψ_k rates

**Advantages of Julia version**:
- Faster optimization (BFGS with automatic differentiation)
- Easy to extend (3+ components)
- Built-in constraint handling
- Comprehensive model comparison

## References

- Fornasiero et al. (2018) for lysine availability dynamics
- Your manuscript on PrP turnover and mixture models
- Jensen's inequality and log-convexity proof

## Troubleshooting

**Optimization fails**:
- Try different initial guesses
- Check data quality (outliers, missing values)
- Reduce number of components

**Poor fits**:
- Check if data spans enough timepoints (need ≥ 6-8 points)
- Verify data covers 2-3 half-lives of slower component
- Consider 3-component model if K=2 insufficient

**Unidentifiable parameters**:
- Add constraint (constrained=true)
- Collect more data points
- Fix one parameter if known a priori
