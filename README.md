# Imputation, calibration, and joint plausibility

Code and data accompanying the manuscript:

> **When Does Imputation Method Matter? Calibratability Governs a Conditional Trade-off Between Calibration and Joint Plausibility**
> S. Spektor and S. Banerjee
> *Submitted to* Statistics in Medicine

This repository contains everything needed to reproduce the analyses in the
paper: three benchmark datasets, the main analysis scripts for each dataset, and
the supporting experiments (semi-synthetic mechanism test, Mahalanobis joint
plausibility, multiple-imputation robustness).

## Contents

**Datasets** (public, originally from the UCI Machine Learning Repository):

- `diabetes.csv` — Pima Indians Diabetes (768 rows, 8 features + outcome)
- `Final-breast_cancer_data-2.csv` — Breast Cancer Wisconsin Diagnostic
- `heart_cleveland.csv` — Cleveland Heart Disease

**Main analysis scripts** (one per dataset):

- `Analysis_pima.R` — full pipeline on Pima
- `Analysis_bcw.R` — full pipeline on Breast Cancer Wisconsin
- `Analysis_heart.R` — full pipeline on Cleveland Heart Disease (mixed-type)

**Supporting analyses**:

- `semisynthetic_mechanism.R` — semi-synthetic experiment isolating
  calibratability vs data-type composition (Section: Mechanism)
- `joint_plausibility_mahalanobis.R` — joint-plausibility metric under two
  independent definitions (convex hull, Mahalanobis)
- `mice_m_pima.R`, `mice_m_bcw.R` — multiple-imputation robustness ($m=20$ vs
  $m=1$, with Rubin pooling)

## Requirements

- R ≥ 4.3 (developed on 4.3.3)
- Packages: `mice`, `glmnet`, `randomForest`, `pROC`, `geometry`, `dplyr`

Exact package versions used to produce the results are listed in
`sessionInfo.txt`.

## To reproduce

Each main analysis script is self-contained. From R or RStudio with this
repository as the working directory:

```r
source("Analysis_pima.R")    # ~10-20 min, 100 seeds
source("Analysis_bcw.R")     # ~15-25 min, 100 seeds
source("Analysis_heart.R")   # ~25-40 min, 100 seeds (mixed-type takes longer)
```

Each script writes its results as an `.rds` file containing the per-seed
metrics tables used in the paper.

Supporting analyses can be run after or independently:

```r
source("semisynthetic_mechanism.R")        # writes semisynthetic_results.csv
source("joint_plausibility_mahalanobis.R") # edit DATASET at top: pima/bcw/heart
source("mice_m_pima.R")                    # multiple-imputation robustness
source("mice_m_bcw.R")
```

## Reproducibility notes

- All Monte Carlo runs use `set.seed(s)` for `s` in `1:100`.
- Calibration slopes are summarized by the **median**, not the mean, because a
  small fraction of seeds on near-separable test splits produce extreme slope
  estimates. The median is robust to these; the bounded ECE corroborates the
  ranking. This is discussed in the manuscript.
- The semi-synthetic experiment introduces missingness only under MCAR at 20%
  to isolate the mechanism question from missingness-pattern effects.

## License

Code: MIT License (see `LICENSE`).

Data: the three datasets are redistributions of publicly available UCI
benchmarks. Original sources and licenses apply.

## Citation

If you use this code, please cite the paper (citation block to be added on
acceptance).

## Contact

S. Spektor / S. Banerjee, Department of Quantitative Sciences, Canisius
University.
