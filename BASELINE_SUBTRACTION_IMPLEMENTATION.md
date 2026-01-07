# Baseline Subtraction (LLQ Correction) Implementation

## Overview
This document describes the baseline subtraction implementation for Figure 2 CRL proteomics analysis. The goal is to subtract the lower limit of quantification (LLQ) from all data points before model fitting and plotting.

## Implementation Summary

### Data Flow
```
Raw Data (data/crl_proteomics.tsv)
    ↓
Julia Analysis (src/fig2_mixture_analysis.jl)
    ├─→ Baseline Subtraction
    ├─→ Model Fitting on Adjusted Data
    ├─→ Save: data/crl/baseline_subtracted_data.tsv
    └─→ Save: data/crl/mixture_params.tsv
    ↓
R Plotting (src/fig2.R)
    ├─→ Load: data/crl/baseline_subtracted_data.tsv
    ├─→ Load: data/crl/mixture_params.tsv
    └─→ Plot Adjusted Data + Model Fits
```

### Algorithm

1. **Determine BQL (Below Quantification Limit)**
   - Load peptide standards from `data/crl_pepstds.tsv`
   - BQL = first back_calc_ng_ml value / 10 = 0.1111...

2. **Calculate Proportion Labeled**
   - For missing heavy values, set heavy = BQL
   - total = light + heavy
   - prop_labeled = heavy / total

3. **Baseline Subtraction (per tissue group)**
   - For each tissue (brain, colon):
     - Identify day 0 samples where heavy == BQL
     - LLQ = max(prop_labeled) among those samples
     - Subtract LLQ from all prop_labeled values in that tissue
     - Floor at 0: prop_labeled = max(0, prop_labeled - LLQ)

4. **Model Fitting**
   - Fit isotopic labeling models on baseline-subtracted data
   - Single-rate and mixture models use adjusted prop_labeled values

5. **Plotting**
   - Plot baseline-subtracted data points
   - Plot model predictions (which are already in adjusted scale)
   - LLQ appears at y=0 (indicated by dashed line)

## Results from Latest Run

### Brain
- LLQ value: 0.024602
- Before adjustment: range 0.0145 to 0.2338, mean 0.1095
- After adjustment: range 0.0 to 0.2092, mean 0.0857
- Data points floored to 0: 7

### Colon
- LLQ value: 0.133111
- Before adjustment: range 0.0783 to 0.5947, mean 0.3027
- After adjustment: range 0.0 to 0.4616, mean 0.1769
- Data points floored to 0: 9

## Files Modified

### Julia
- **src/fig2_mixture_analysis.jl**
  - Fixed BQL calculation to use peptide standards (was hardcoded 0.001)
  - Added detailed diagnostic output for baseline subtraction
  - Saves baseline-subtracted data to `data/crl/baseline_subtracted_data.tsv`
  - Models are fit on adjusted data

### R
- **src/fig2.R**
  - Loads baseline-subtracted data from Julia output
  - Falls back to figS2.R data if Julia output not available
  - Plots use the adjusted prop_labeled values
  - Model predictions are plotted without offset

- **src/figS2.R**
  - Enhanced comments documenting baseline subtraction logic
  - Implements same algorithm as Julia (for fallback)

### Module
- **src/TurnoverModels.jl**
  - No changes needed (already supports fitting on any data)

## How to Use

### Step 1: Run Julia Analysis
```bash
julia --project=. src/fig2_mixture_analysis.jl
```

This will:
- Load raw data
- Apply baseline subtraction
- Fit models on adjusted data
- Save `data/crl/baseline_subtracted_data.tsv`
- Save `data/crl/mixture_params.tsv`
- Print detailed diagnostics

### Step 2: Generate Figures in R
```r
source("src/halflife_figures.R")
```

Or just Figure 2:
```r
source("src/utils.R")
source("src/funcs.R")
source("src/figS2.R")  # Creates tissue_meta
source("src/fig2.R")   # Loads Julia output
```

## Verification

The implementation ensures:
1. ✅ Same BQL value used in Julia and R
2. ✅ Same baseline subtraction algorithm in Julia and R
3. ✅ Models fit on baseline-subtracted data
4. ✅ Plots show baseline-subtracted data
5. ✅ Model predictions match the adjusted scale
6. ✅ LLQ appears at y=0 in plots

## Output Files

### data/crl/baseline_subtracted_data.tsv
Columns:
- animal: Animal ID
- day: Time point
- tissue: brain or colon
- light: Light peptide intensity
- heavy: Heavy peptide intensity (with BQL substitution)
- total: light + heavy
- prop_labeled: Baseline-subtracted proportion (LLQ removed, floored at 0)

### data/crl/mixture_params.tsv
Model fitting results including:
- Single-rate half-life
- Mixture model parameters (fast/slow pools, proportions)
- AIC/BIC for model comparison

## Notes

- The baseline subtraction is tissue-specific (each tissue has its own LLQ)
- LLQ is determined from day 0 samples where heavy == BQL
- This represents the maximum background signal for that tissue
- After adjustment, LLQ = 0 for all tissues
- Any negative values after subtraction are floored to 0
