# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a research project analyzing protein half-life dynamics using isotopic labeling data. The codebase combines R for data analysis and figure generation with Julia for mathematical modeling and Bayesian inference of protein decay kinetics.

## Code Architecture

### Dual Language Structure

- **R code** (primary): Data processing, statistical analysis, and figure generation
- **Julia code**: Mathematical modeling of protein decay, particularly mixture models with multiple decay rates

### R Code Organization

**Entry point**: `src/halflife_figures.R` - orchestrates all figure generation

**Core utilities** (loaded first):
- `src/utils.R`: Global constants, data loading, output streams, helper functions
- `src/funcs.R`: Mathematical models for protein decay and isotopic labeling

**Figure scripts**: Individual files (`fig1.R`, `fig2.R`, etc.) generate specific manuscript figures and supplementary figures (`figS1.R`, etc.)

### Key Mathematical Models

The project implements two types of differential equation models for protein turnover:

#### 1. Isotopic Labeling Models (Figure 4)
Used when measuring incorporation of heavy isotope-labeled amino acids into newly synthesized protein.

**Single-rate**: `proportion_labeled_single(thalf, t, avails_func)`
- Tracks what fraction of protein molecules contain heavy isotope label
- Depends on amino acid availability function (e.g., `free_lysine(t)`)

**Mixture**: `proportion_labeled_mixture(thalves, proportions, t, avails_func)`
- Models heterogeneous protein pools with different turnover rates
- Each pool independently accumulates label according to its half-life

**Key R function**: `Lt(avails, t, lambda)` in `src/utils.R` - implements this model using Euler integration with `dt=0.01`

**Availability function**: `free_lysine(t)` based on Fornasiero 2018 - models the bi-exponential decay of free lysine pool after pulse labeling

#### 2. Residual Decay Models (Figure 3)
Used when measuring protein abundance after RNA knockdown (e.g., ASO treatment).

**Single-rate**: `residual_decay_single(R_interp, lambda, t)`
- Models protein decay when mRNA production is reduced
- Differential equation: dP/dt = λ × R(t) - λ × P(t)
- R(t) is the RNA time course (measured), P(t) is protein level

**Mixture**: `residual_decay_mixture(R_interp, lambdas, weights, t)`
- Models heterogeneous protein pools with different decay rates
- Each pool decays independently but shares the same RNA input

**Key R function**: `Pt(R, t, lambda)` and `Pt_mixture(R, t, lambdas, weights)` in `src/funcs.R`

**Critical difference**: Residual decay requires RNA measurements; isotopic labeling requires availability function

### Julia Mathematical Framework

The project uses a modular architecture:

**`src/TurnoverModels.jl`** - Reusable module providing:
- Optimized numerical integration for both model types
- Multi-start optimization for robust parameter estimation
- Model comparison tools (AIC, BIC calculations)
- Runs test for detecting systematic residual patterns
- Effective half-life calculation for mixture models

**`src/fig3_mixture_analysis.jl`** - ASO knockdown analysis:
- Loads ELISA and RNA data from ASO experiments
- Fits residual decay mixture models
- Outputs parameters to `data/aso/mixture_params.tsv`

**`src/fig4_mixture_analysis.jl`** - Isotopic labeling analysis:
- Loads heavy/light peptide measurements from pulse-chase experiments
- Fits isotopic labeling mixture models
- Computes genotypic differences (FFI vs controls)
- Outputs to `data/iqp/mixture_params.tsv` and `genotypic_diffs.tsv`

## Running the Code

### R Analysis

Generate all figures:
```r
# From R or RStudio
source("src/halflife_figures.R")
```

This will:
- Load dependencies (tidyverse, janitor, openxlsx, smoother, plotrix, magick, minpack.lm)
- Process data from `data/` directory
- Generate figures in `display_items/`
- Create Excel supplement (`display_items/supplement.xlsx`) with numbered tables
- Write statistics to `display_items/stats_for_text.txt`

### Julia Modeling

The Julia environment is defined in `Project.toml` with key dependencies:
- Optim.jl for maximum likelihood fitting
- Turing.jl for Bayesian inference
- DataFrames.jl, CSV.jl for data handling
- Interpolations.jl for RNA time series interpolation
- HypothesisTests.jl for statistical tests

#### Module Architecture

The project uses a shared `TurnoverModels.jl` module that provides:
- **Isotopic labeling models**: `proportion_labeled_single`, `proportion_labeled_mixture`
- **Residual decay models**: `residual_decay_single`, `residual_decay_mixture`
- **Fitting functions**: for both model types with multi-start optimization
- **Statistical utilities**: AIC, BIC calculations, runs tests, effective half-life

#### Running Figure-Specific Analyses

**Figure 3 (ASO knockdown, residual decay)**:
```bash
julia --project=. src/fig3_mixture_analysis.jl
```
Outputs: `data/aso/mixture_params.tsv`

**Figure 4 (isotopic labeling)**:
```bash
julia --project=. src/fig4_mixture_analysis.jl
```
Outputs: `data/iqp/mixture_params.tsv`, `data/iqp/genotypic_diffs.tsv`

The R scripts (`fig3.R`, `fig4.R`) automatically load these pre-computed parameters for plotting.

## Data Structure

- `data/`: TSV files with experimental measurements
- `data/iqp/`: IQ Proteomics processed data (key: `combined_pivot_tab.tsv`, `meta.tsv`)
- `data/gtex/`: GTEx expression data for tissue-specific analysis
- `display_items/`: Generated figures (PNG) and supplementary tables (XLSX, TSV)

## Output Conventions

**Figures**: PNG format at 300 DPI in `display_items/`, named as `figure-1.png`, `figure-s1.png`, etc.

**Supplementary tables**:
- Excel workbook with sequential sheets (s01, s02, etc.)
- Tab-separated copies (`table-s01.tsv`, etc.) for version control
- Table of contents auto-generated in first worksheet

**Statistics**: Text file `stats_for_text.txt` stores computed statistics referenced in manuscript

## Development Notes

- The R code uses `tell_user()` to print progress messages to stderr
- Confidence intervals use 95% by default via `upper()` and `lower()` functions
- Plotting uses transparency constant `ci_alpha = 0.35` for confidence interval shading
- The `write_supp_table()` function handles dual Excel/TSV output automatically
