# SCED: Source Contribution Error Distribution

**Probabilistic assessment of estimation errors in isotope-based root water uptake partitioning**

## 1. Overview

Bayesian isotope mixing models such as MixSIAR estimate the contributions of different soil water sources to plant water uptake and provide posterior credible intervals for those estimates. However, a credible interval's width does not directly indicate how far the estimated contribution is from its unknown true value.

The **Source Contribution Error Distribution (SCED)** model complements MixSIAR by predicting a *conditional probability distribution* for the absolute estimation error (AE) of each source contribution estimate. SCED combines the width of the 95% posterior credible interval, the number of candidate sources, and the geometry of the dual-isotope mixing space. It can report predicted error quantiles and the probability that AE exceeds a user-specified threshold. These are **model-predicted risks, not observed errors**.

This repository provides the ready-to-use SCED prediction tool, an apple-tree application dataset, the final source-level analysis data, and the code used to generate and select simulated isotope mixing scenarios.

## 2. Repository structure

```text
SCED/
├── 01_SCED_Model_and_Example/
│   ├── SCED_Prediction_Tool.R
│   ├── Wang_HESS_MixSIAR_all_48_cases_wide.xlsx
│   └── MixSIAR_run/
│       ├── Case_List.csv
│       ├── G001mix.csv, G001source.csv, G001_rjags.Rdata, ...
│       └── [other available case inputs and posterior files]
├── 02_Analysis_Data/
│   └── MixSIAR_3_2_Source_998.csv
├── 03_Simulation/
│   ├── SCED_generate_and_select_scenarios.R
│   └── inputs/
│       ├── contexts.csv
│       ├── pooled_fits.csv
│       ├── residuals.csv
│       ├── source_sd_pairs.csv
│       └── single_tree_sd_pairs.csv
└── README.md
```

The tree shows the intended release layout. Generated reports and intermediate simulation outputs do not need to be committed to the repository. See Section 8 for the licensing and citation status.

## 3. Using the SCED prediction tool

### 3.1. Required inputs

SCED **does not fit MixSIAR**. Each case must already have a completed, ungrouped MixSIAR posterior saved as an Rjags `.Rdata` file, together with the isotope inputs used for that fit.

Create a `Case_List.csv` manifest with these four columns:

| Column | Contents |
| --- | --- |
| `Case_ID` | Unique case identifier |
| `Mixture_File` | CSV with paired mixture-water (e.g., xylem-water) δ²H and δ¹⁸O observations |
| `Source_File` | CSV with source-water isotope measurements, or source means, SDs, and sample sizes |
| `Rjags_File` | Completed MixSIAR Rjags posterior (`.Rdata`) containing source contribution draws |

For example, if the manifest and the referenced files are in the same `MixSIAR_run/` directory:

```csv
Case_ID,Mixture_File,Source_File,Rjags_File
G016,G016mix.csv,G016source.csv,G016_rjags.Rdata
G017,G017mix.csv,G017source.csv,G017_rjags.Rdata
```

File paths are resolved relative to the **manifest's directory** (absolute paths are also supported). One manifest can include multiple independent cases. The supplied `MixSIAR_run/Case_List.csv` lists the available example runs.

By default, the code expects the hydrogen and oxygen tracer columns to be named `H` and `O`, respectively. If your MixSIAR files use different tracer names, edit `TRACER_D2H` and `TRACER_D18O` at the beginning of `SCED_Prediction_Tool.R`. The source CSV may contain either paired raw observations or MixSIAR-compatible source means, SDs, and sample sizes. The source labels and posterior must correspond to the **same MixSIAR fit**.

### 3.2. Run predictions

From the repository root, run the following command (shown on one line for Windows and other shells):

```bash
Rscript "01_SCED_Model_and_Example/SCED_Prediction_Tool.R" "01_SCED_Model_and_Example/MixSIAR_run/Case_List.csv" "01_SCED_Model_and_Example/SCED_Batch_Report"
```

To analyze your own MixSIAR cases, replace the manifest path with your `Case_List.csv` and choose an output directory. The default exceedance threshold is `AE_THRESHOLD <- 0.20`, defined near the top of the script; change it to another value strictly between 0 and 1 as needed. A threshold of 0.20 means an **absolute difference of 20 percentage points** in source contribution, not a relative error of 20%.

The SCED coefficients, predictor-standardization parameters, and fixed Beta precision are embedded in the script. **No SCED refitting is needed** for prediction.

### 3.3. Outputs and interpretation

A batch run writes:

```text
SCED_Batch_Report/
├── SCED_All_Source_Predictions.csv
├── SCED_Batch_Report.html
└── Cases/
    └── <Case_ID>/
        ├── SCED_Probability_Report.html
        ├── Mixing_Space.svg
        └── SCED_Probability_Curves.svg
```

The combined CSV contains one record per source for successfully predicted cases. Important fields are:

| Field | Interpretation |
| --- | --- |
| `MixSIAR_mean` | Posterior mean estimated source contribution |
| `MixSIAR_q025`, `MixSIAR_q975` | Limits of the MixSIAR 95% credible interval |
| `W95` | Width of that credible interval, on a 0–1 scale |
| `Predicted_mean_AE` | SCED-predicted conditional mean absolute error |
| `AE_Q90`, `AE_Q95` | SCED-predicted 90th and 95th percentiles of AE |
| `P_AE_gt_threshold` | Predicted probability that AE exceeds `AE_threshold` |
| `AE_threshold` | Error threshold used for the exceedance prediction |

For example, **`AE_Q90 = 0.20`** means that SCED assigns a 90% probability to an absolute contribution error **no greater than 20 percentage points** for that source estimate. It does **not** mean that the true error is known, that the estimate is guaranteed to meet the bound, or that the MixSIAR credible interval has a 90% coverage probability.

Failed or inapplicable cases are recorded in the batch HTML overview without stopping the remaining cases. Each successful case also receives an illustrated report showing the isotope mixing space and predicted cumulative AE distributions.

### 3.4. Model predictors

SCED uses the following predictors (calculated automatically from the supplied MixSIAR inputs and posterior):

| Predictor | Meaning |
| --- | --- |
| `NOS` | Number of candidate water sources (3–7) |
| `W95` | Width of the posterior 95% credible interval |
| `S_nearest` | Standardized isotopic distance from a source to its nearest other source; log10-transformed in the model |
| `Rpos` | Relative radial position of the mixture within the source isotope convex hull; 0 at the defined center and 1 at the boundary |
| `Dmix` | Isotope-range-standardized distance between a source and the mixture (`DFM_scaled` in the code) |

The SCED mean equation also includes the `W95 × Rpos` interaction. The final model uses a constant-precision Beta distribution; details of the definitions and fitted structure are provided in the associated manuscript.

## 4. Apple-tree application examples

`Wang_HESS_MixSIAR_all_48_cases_wide.xlsx` retains **48 apple-tree isotope mixing cases** derived from the study by Wang et al. (2023). The supplied `MixSIAR_run/` directory currently contains the MixSIAR input files and completed Rjags results for **24 cases (`G001`–`G024`)**, rather than all 48. The available prior SCED batch analysis produced predictions for **13 cases**; the other cases should not be interpreted as having valid SCED predictions merely because they appear in the workbook or have MixSIAR results.

The manuscript illustrates two examples, **G016 and G017**, where similar estimated water-source contributions and posterior interval widths correspond to different SCED-predicted error risks. The complete workbook is retained to make the broader application dataset available, rather than restricting the repository to the two illustrated cases.

Before redistributing any third-party observational data, verify the original data-sharing terms and acknowledge the primary source: Wang et al. (2023), *Hydrology and Earth System Sciences*, 27, 123–137, https://doi.org/10.5194/hess-27-123-2023.

## 5. Analysis data and scenario simulation

### 5.1. Final analysis dataset

`02_Analysis_Data/MixSIAR_3_2_Source_998.csv` contains **4,994 source-level estimates from 998 retained simulation scenarios**. The true source contributions are known by construction, allowing observed simulation errors to be calculated as:

`AE = abs(p_post_mean - p_true)`

Selected columns include:

| Column(s) | Contents |
| --- | --- |
| `scenario_id`, `run_id`, `NOS`, `source_index` | Scenario and source identifiers |
| `p_true`, `p_post_mean` | Simulated true contribution and MixSIAR posterior mean |
| `q025`, `q975`, `W95`, `AE` | Posterior interval endpoints, interval width, and absolute estimation error |
| `S_nearest`, `Rpos_obs`, `DFM_scaled` | Source and mixture geometry descriptors |
| `source_mean_d18O`, `source_mean_d2H` | Source isotope means |
| `source_SD_d18O`, `source_SD_d2H` | Source isotope standard deviations |
| `Rhat_max_p`, `ESS_min_p` | Summary convergence diagnostics |
| `scenario_weight_per_source` | Source-level weight assigning equal total weight to each scenario |

Additional columns are preserved from the analysis workflow. In particular, `S_min` is a legacy/auxiliary geometry variable; **the released SCED predictor uses `S_nearest`, not `S_min`**. No separate model A–C fitting scripts or cross-validation output files are included in this streamlined repository.

### 5.2. Generate and select simulation scenarios

The integrated script `03_Simulation/SCED_generate_and_select_scenarios.R` executes three original stages consecutively:

1. **Generation:** use empirical isotope parameters in `inputs/` to simulate 15,000 candidate mixing scenarios (3 reference datasets × 5 NOS levels × 1,000 candidates), with known source contributions and six paired mixture-water observations per scenario.
2. **Geometry audit:** calculate candidate-level isotope mixing-space descriptors and assess coverage, without changing the candidates or selecting cases.
3. **Selection:** select **1,000 scenarios** (200 per NOS level) using a random baseline of 140 scenarios per NOS followed by geometry-coverage augmentation; the selection does not use MixSIAR estimation errors.

The study's original simulation generated **14,996 complete candidate scenarios**, from which the 1,000 cases were selected. The program retains the original fixed generation and selection seeds (`20261023` and `20261024`, respectively). It creates its own `outputs_v3_2/`, `coverage_diagnostics/`, and `selection_final_V1/` output directories beside the script.

Run from the repository root:

```bash
Rscript "03_Simulation/SCED_generate_and_select_scenarios.R"
```

This simulation can be computationally demanding. The workflow **does not run MixSIAR**, assess MCMC convergence, re-create the posterior estimates for the 998 retained scenarios, or refit/evaluate SCED. The supplied analysis CSV contains the final data used in those later stages. Re-running the generator also requires the same input parameter files and a compatible R/software environment; newly generated scenarios should not automatically replace the published analysis dataset.

## 6. Applicability and limitations

The current SCED implementation is designed for:

- Paired δ²H and δ¹⁸O mixing data with **3–7 candidate sources**.
- **Ungrouped** MixSIAR `p.global` posteriors, not grouped or individual-effect models.
- **Zero isotope discrimination** (both mean and SD), with no concentration dependence. The tool does not inspect a discrimination file; the user must verify this assumption.
- MixSIAR settings comparable to those used in calibration, including the **process × residual** error structure. The training simulations used six paired mixture observations per scenario and nominal source sample size six.
- Mixture isotope means lying **inside the convex hull** defined by the source isotope means.

The tool checks basic input validity, geometric feasibility, selected posterior diagnostics, and potential extrapolation. Some differences from training settings produce warnings; some violations prevent prediction. These checks are not an independent validation of a new dataset.

SCED gives **probabilistic information about possible estimation error**. It cannot determine the actual error of an individual field estimate when the true source contribution is unknown. Predictions for data substantially outside the simulation conditions should be interpreted cautiously.

## 7. Software requirements

- **SCED prediction:** base R for a CSV manifest; the optional `readxl` package is needed for an Excel manifest. A completed MixSIAR `.Rdata` posterior is required, but SCED prediction does not re-run MixSIAR.
- **Simulation workflow:** R **3.6.0 or later**, with `readxl` and `writexl` installed.

Install the simulation dependencies with:

```r
install.packages(c("readxl", "writexl"))
```

The bundled R scripts can be run from a terminal with `Rscript`. The simulation script is intended to run as a single file; users do not need to `source()` its three original components separately.

## 8. Citation and license

**Associated manuscript:** *Probabilistic prediction of estimation errors in isotope-based root water uptake partitioning: Incorporating posterior uncertainty and isotope mixing-space geometry*. Full author, journal, and DOI details will be added when available.

If you use the SCED predictions or the simulation dataset in published work, please cite the associated manuscript once its bibliographic details are available. For the apple-tree observational material, also cite the original Wang et al. (2023) study identified in Section 4.

**License status:** A software license and a separate data-reuse statement must be selected and added **before public release**. This README does not itself grant permission to redistribute third-party observations. Do not assume the code, original empirical data, and derived application material have identical reuse rights.
