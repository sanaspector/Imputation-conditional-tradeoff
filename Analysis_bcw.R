# =============================================================================
#  Analysis_v2.R
#  Comparative imputation study: iterative LASSO vs missForest vs MICE-PMM
#  on the Pima Indians Diabetes dataset.
#
#  Key improvements over the original Analysis.R:
#    1. Iterative (chained) LASSO so the comparison with missForest is fair.
#    2. MICE-PMM added as the field-standard baseline.
#    3. 100-seed Monte Carlo replication with paired Wilcoxon tests
#       and 95% percentile bootstrap CIs on every metric.
#    4. Three missingness mechanisms (MCAR / MAR / MNAR) at 10%, 20%, 30%.
#    5. Expanded metrics:
#         - Imputation accuracy: NRMSE, MAE, Pearson r
#         - Univariate plausibility: hard-bound violations
#         - Joint plausibility: % of imputed (glucose, insulin) pairs
#             falling outside the observed bivariate convex hull
#         - Distributional fidelity: Wasserstein-1 distance
#         - Downstream: AUC, Brier score, calibration slope/intercept,
#             expected calibration error (ECE), decision-curve net benefit
#  Author: S. Spektor, May 2026
# =============================================================================

suppressPackageStartupMessages({
  required_packages <- c(
    "tidyverse", "glmnet", "randomForest", "mice",
    "pROC", "transport", "geometry", "boot"
  )
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
})

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
CFG <- list(
  data_path        = "Final-breast_cancer_data-2.csv",
  na_vars          = c("texture_mean", "smoothness_mean", "concavity_se",
                       "concave.points_se", "symmetry_se",
                       "fractal_dimension_se", "area_worst",
                       "symmetry_worst"),
  outcome_var      = "diagnosis",
  miss_rates       = c(0.10, 0.20, 0.30),
  mechanisms       = c("MCAR", "MAR", "MNAR"),
  n_seeds          = 5,
  max_iter_chain   = 10,
  conv_tol         = 1e-3,
  rf_ntree         = 100,
  dca_thresholds   = seq(0.05, 0.50, by = 0.05),
  # Data-driven "out-of-support" bounds: observed range with 5% expansion.
  # Reframed in the paper as out-of-support, NOT clinical plausibility,
  # because these are image-processing morphometric features, not values
  # a clinician reads directly off a report.
  clinical_bounds  = list(
    texture_mean         = c(9.22, 41.24),      # 9.71 - 5%  to  39.28 + 5%
    smoothness_mean      = c(0.0500, 0.1716),   # 0.0526 - 5%  to  0.1634 + 5%
    concavity_se         = c(0.0000, 0.4158),   # min stays at 0; 0.396 + 5%
    `concave.points_se`  = c(0.0000, 0.0554),   # min stays at 0; 0.0528 + 5%
    symmetry_se          = c(0.0075, 0.0830),
    fractal_dimension_se = c(0.0009, 0.0313),
    area_worst           = c(175.94, 4466.70),
    symmetry_worst       = c(0.1487, 0.6970)
  )
)


# -----------------------------------------------------------------------------
# 1. Load and clean data
# -----------------------------------------------------------------------------
load_bcw <- function(path = CFG$data_path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  complete_data <- df[complete.cases(df), ]
  cat(sprintf("Loaded %d rows, %d complete cases retained.\n",
              nrow(df), nrow(complete_data)))
  complete_data
}

# -----------------------------------------------------------------------------
# 2. Missingness simulation: MCAR, MAR, MNAR
# -----------------------------------------------------------------------------
simulate_missing <- function(data, vars, rate, mechanism, seed) {
  set.seed(seed)
  n   <- nrow(data)
  out <- data
  mask <- matrix(FALSE, nrow = n, ncol = length(vars),
                 dimnames = list(NULL, vars))

  for (v in vars) {
    n_miss <- floor(rate * n)
    if (mechanism == "MCAR") {
      idx <- sample.int(n, n_miss)
    } else if (mechanism == "MAR") {
      # Probability of missingness depends on Age and Pregnancies (observed).
      # Older / higher-parity patients are more likely to have a missing value -
      # plausible for retrospective EHR cohorts.
      lp <- scale(data$texture_mean)[, 1] + scale(data$area_worst)[, 1]
      p  <- plogis(lp - mean(lp))
      idx <- sample.int(n, n_miss, prob = p)
    } else if (mechanism == "MNAR") {
      # Probability of missingness depends on the variable's OWN value.
      # E.g. insulin tends to be missing when it is high (not measured in
      # asymptomatic patients) - the classic informative missingness pattern.
      x  <- data[[v]]
      p  <- plogis(scale(x)[, 1])
      idx <- sample.int(n, n_miss, prob = p)
    } else {
      stop("Unknown mechanism: ", mechanism)
    }
    out[[v]][idx]      <- NA
    mask[idx, v]       <- TRUE
  }
  list(data = out, mask = mask)
}

# -----------------------------------------------------------------------------
# 3. Imputation methods
# -----------------------------------------------------------------------------

# --- 3a. Iterative (chained) LASSO -------------------------------------------
# This is the critical fix vs the original Analysis.R.  We initialise with
# median imputation, then cycle through variables refitting cv.glmnet on
# the FULLY-FILLED matrix, until convergence.  This is the apples-to-apples
# analogue of missForest's iterative scheme.
impute_lasso_chained <- function(df_miss, target_vars,
                                 max_iter = CFG$max_iter_chain,
                                 tol = CFG$conv_tol) {
  df_work <- df_miss
  miss_idx <- lapply(target_vars, function(v) which(is.na(df_work[[v]])))
  names(miss_idx) <- target_vars

  # Median initialisation
  for (v in target_vars) {
    med <- median(df_work[[v]], na.rm = TRUE)
    df_work[[v]][miss_idx[[v]]] <- med
  }

  prev_vals <- lapply(target_vars, function(v) df_work[[v]][miss_idx[[v]]])
  names(prev_vals) <- target_vars

  for (it in seq_len(max_iter)) {
    max_change <- 0
    for (v in target_vars) {
      if (length(miss_idx[[v]]) == 0) next
      predictors <- setdiff(names(df_work), v)
      # Train on rows where v was originally observed
      train_rows <- setdiff(seq_len(nrow(df_work)), miss_idx[[v]])
      x_train <- as.matrix(df_work[train_rows, predictors, drop = FALSE])
      y_train <- df_work[[v]][train_rows]
      x_pred  <- as.matrix(df_work[miss_idx[[v]], predictors, drop = FALSE])

      fit <- tryCatch(
        cv.glmnet(x_train, y_train, alpha = 1,
                  nfolds = min(10, length(y_train) %/% 5)),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      new_vals <- as.numeric(predict(fit, newx = x_pred, s = "lambda.1se"))
      change <- max(abs(new_vals - prev_vals[[v]])) /
                (sd(df_work[[v]], na.rm = TRUE) + 1e-9)
      max_change <- max(max_change, change)
      df_work[[v]][miss_idx[[v]]] <- new_vals
      prev_vals[[v]] <- new_vals
    }
    if (max_change < tol) break
  }
  df_work
}

# --- 3b. missForest-style iterative Random Forest ----------------------------
impute_rf_chained <- function(df_miss, target_vars,
                              max_iter = CFG$max_iter_chain,
                              tol = CFG$conv_tol,
                              ntree = CFG$rf_ntree) {
  df_work <- df_miss
  miss_idx <- lapply(target_vars, function(v) which(is.na(df_work[[v]])))
  names(miss_idx) <- target_vars

  for (v in target_vars) {
    df_work[[v]][miss_idx[[v]]] <- median(df_work[[v]], na.rm = TRUE)
  }
  prev_vals <- lapply(target_vars, function(v) df_work[[v]][miss_idx[[v]]])
  names(prev_vals) <- target_vars

  for (it in seq_len(max_iter)) {
    max_change <- 0
    for (v in target_vars) {
      if (length(miss_idx[[v]]) == 0) next
      predictors <- setdiff(names(df_work), v)
      train_rows <- setdiff(seq_len(nrow(df_work)), miss_idx[[v]])
      fit <- tryCatch(
        randomForest(x = df_work[train_rows, predictors, drop = FALSE],
                     y = df_work[[v]][train_rows],
                     ntree = ntree),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      new_vals <- predict(fit,
                          newdata = df_work[miss_idx[[v]], predictors,
                                            drop = FALSE])
      change <- max(abs(new_vals - prev_vals[[v]])) /
                (sd(df_work[[v]], na.rm = TRUE) + 1e-9)
      max_change <- max(max_change, change)
      df_work[[v]][miss_idx[[v]]] <- new_vals
      prev_vals[[v]] <- new_vals
    }
    if (max_change < tol) break
  }
  df_work
}

# --- 3c. MICE with predictive mean matching (field standard) ----------------
impute_mice_pmm <- function(df_miss, m = 1, maxit = 10) {
  invisible(capture.output({
    imp <- mice(df_miss, m = m, method = "pmm",
                maxit = maxit, printFlag = FALSE)
  }))
  complete(imp, 1)
}

# -----------------------------------------------------------------------------
# 4. Evaluation metrics
# -----------------------------------------------------------------------------
accuracy_metrics <- function(imputed, truth, mask, vars) {
  rows <- list()
  for (v in vars) {
    idx <- which(mask[, v])
    if (length(idx) == 0) next
    t <- truth[[v]][idx]; i <- imputed[[v]][idx]
    ok <- is.finite(t) & is.finite(i)
    if (sum(ok) == 0) next
    t <- t[ok]; i <- i[ok]
    rmse  <- sqrt(mean((t - i)^2))
    nrmse <- rmse / sd(t)
    mae   <- mean(abs(t - i))
    r     <- suppressWarnings(cor(t, i))
    rows[[v]] <- data.frame(Variable = v, NRMSE = nrmse,
                            MAE = mae, Pearson = r)
  }
  do.call(rbind, rows)
}

plausibility_univariate <- function(imputed, mask, vars,
                                    bounds = CFG$clinical_bounds) {
  rows <- list()
  for (v in vars) {
    idx <- which(mask[, v])
    if (length(idx) == 0) next
    x <- imputed[[v]][idx]
    b <- bounds[[v]]
    viol <- mean(x < b[1] | x > b[2]) * 100
    neg  <- mean(x < 0) * 100  # the headline "negative insulin" stat
    rows[[v]] <- data.frame(Variable = v,
                            ViolationPct = viol,
                            NegativePct  = neg)
  }
  do.call(rbind, rows)
}

# Joint plausibility: % of imputed (Glucose, Insulin) points falling outside
# the convex hull of the originally observed (Glucose, Insulin) pairs.
# FIXED CODE FOR BREAST CANCER DATASET
plausibility_joint <- function(imputed, truth, mask,
                               v1 = "concavity_se", v2 = "concave.points_se") {
  # Safe boundary check: verify if the variables exist in the matrix columns
  if (!(v1 %in% colnames(mask)) || !(v2 %in% colnames(mask))) {
    return(data.frame(JointOutsidePct = NA_real_))
  }
  
  obs_rows <- which(!mask[, v1] & !mask[, v2])
  imp_rows <- which(mask[, v1] | mask[, v2])
  
  if (length(obs_rows) < 4 || length(imp_rows) == 0)
    return(data.frame(JointOutsidePct = NA_real_))
  
  obs_pts <- as.matrix(truth[obs_rows, c(v1, v2)])
  imp_pts <- as.matrix(imputed[imp_rows, c(v1, v2)])
  
  # Shield the entire geometric operation to prevent collinearity errors 
  # from halting the Monte Carlo simulation loop
  res_pct <- tryCatch({
    hull <- geometry::convhulln(obs_pts, options = "Tv")
    inside <- geometry::inhulln(hull, imp_pts)
    mean(!inside) * 100
  }, error = function(e) {
    NA_real_  # Gracefully returns NA if data points are degenerate on this seed
  })
  
  data.frame(JointOutsidePct = res_pct)
}

wasserstein_fidelity <- function(imputed, truth, vars) {
  rows <- list()
  for (v in vars) {
    a <- truth[[v]]; b <- imputed[[v]]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) == 0 || length(b) == 0) next
    w <- tryCatch(transport::wasserstein1d(a, b),
                  error = function(e) NA_real_)
    rows[[v]] <- data.frame(Variable = v, Wasserstein1 = w)
  }
  do.call(rbind, rows)
}

# Downstream prediction with logistic regression.
# Returns AUC, Brier, calibration slope/intercept, ECE, and net benefit.
downstream_metrics <- function(imputed_train, imputed_test,
                               outcome_var = CFG$outcome_var,
                               thresholds = CFG$dca_thresholds) {
  fmla <- as.formula(paste(outcome_var, "~ ."))
  fit <- suppressWarnings(
    glm(fmla, data = imputed_train, family = binomial())
  )
  p <- suppressWarnings(predict(fit, newdata = imputed_test, type = "response"))
  y <- imputed_test[[outcome_var]]
  ok <- is.finite(p) & is.finite(y)
  p <- pmin(pmax(p[ok], 1e-6), 1 - 1e-6); y <- y[ok]

  auc <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
  brier <- mean((p - y)^2)

  # Calibration via logistic recalibration on the logit of predicted prob
  logit_p <- log(p / (1 - p))
  cal_fit <- suppressWarnings(glm(y ~ logit_p, family = binomial()))
  cal_intercept <- coef(cal_fit)[1]
  cal_slope     <- coef(cal_fit)[2]

  # Expected calibration error (10 equal-width bins)
  bins <- cut(p, breaks = seq(0, 1, length.out = 11), include.lowest = TRUE)
  ece <- sum(tapply(seq_along(p), bins, function(ii) {
    if (length(ii) == 0) return(0)
    abs(mean(p[ii]) - mean(y[ii])) * length(ii) / length(p)
  }), na.rm = TRUE)

  # Decision-curve net benefit averaged over thresholds
  nb <- sapply(thresholds, function(th) {
    pred_pos <- p >= th
    if (sum(pred_pos) == 0) return(0)
    tp <- sum(pred_pos & y == 1) / length(y)
    fp <- sum(pred_pos & y == 0) / length(y)
    tp - fp * (th / (1 - th))
  })
  net_benefit <- mean(nb)

  data.frame(AUC = auc, Brier = brier,
             CalSlope = cal_slope, CalIntercept = cal_intercept,
             ECE = ece, NetBenefit = net_benefit)
}

# -----------------------------------------------------------------------------
# 5. One Monte-Carlo replicate
# -----------------------------------------------------------------------------
one_replicate <- function(complete_data, mechanism, rate, seed) {
  sim <- simulate_missing(complete_data, CFG$na_vars, rate, mechanism, seed)
  df_miss <- sim$data; mask <- sim$mask

  # Train/test split on the COMPLETE rows (the truth), then propagate the
  # missingness mask into the training split only.  This mimics the realistic
  # setting where the analyst imputes the training data and a clinician later
  # has fully observed test data (or vice versa - any consistent choice works).
  set.seed(seed + 1)
  n <- nrow(complete_data)
  train_idx <- sample.int(n, floor(0.7 * n))
  test_idx  <- setdiff(seq_len(n), train_idx)

  train_truth <- complete_data[train_idx, ]
  test_truth  <- complete_data[test_idx, ]
  train_miss  <- df_miss[train_idx, ]
  mask_train  <- mask[train_idx, , drop = FALSE]

  methods <- list(
    LASSO_chained = impute_lasso_chained(train_miss, CFG$na_vars),
    RF_chained    = impute_rf_chained(train_miss, CFG$na_vars),
    MICE_PMM      = impute_mice_pmm(train_miss)
  )

  results <- list()
  for (m_name in names(methods)) {
    imp <- methods[[m_name]]
    acc <- accuracy_metrics(imp, train_truth, mask_train, CFG$na_vars)
    plu <- plausibility_univariate(imp, mask_train, CFG$na_vars)
    jpl <- plausibility_joint(imp, train_truth, mask_train)
    wd  <- wasserstein_fidelity(imp, train_truth, CFG$na_vars)
    dwn <- downstream_metrics(imp, test_truth)

    results[[m_name]] <- list(
      acc = cbind(Method = m_name, Mechanism = mechanism, Rate = rate,
                  Seed = seed, acc),
      plu = cbind(Method = m_name, Mechanism = mechanism, Rate = rate,
                  Seed = seed, plu),
      jpl = cbind(Method = m_name, Mechanism = mechanism, Rate = rate,
                  Seed = seed, jpl),
      wd  = cbind(Method = m_name, Mechanism = mechanism, Rate = rate,
                  Seed = seed, wd),
      dwn = cbind(Method = m_name, Mechanism = mechanism, Rate = rate,
                  Seed = seed, dwn)
    )
  }
  results
}

# -----------------------------------------------------------------------------
# 6. Driver: Monte Carlo experiment
# -----------------------------------------------------------------------------
run_experiment <- function(complete_data,
                           mechanisms = CFG$mechanisms,
                           rates      = CFG$miss_rates,
                           n_seeds    = CFG$n_seeds) {
  acc_all <- plu_all <- jpl_all <- wd_all <- dwn_all <- list()
  conditions <- expand.grid(mech = mechanisms, rate = rates,
                            stringsAsFactors = FALSE)
  total <- nrow(conditions) * n_seeds
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  k <- 0
  for (ci in seq_len(nrow(conditions))) {
    mech <- conditions$mech[ci]; rate <- conditions$rate[ci]
    for (s in seq_len(n_seeds)) {
      k <- k + 1; setTxtProgressBar(pb, k)
      rep <- tryCatch(one_replicate(complete_data, mech, rate, seed = s),
                      error = function(e) {
                        message(sprintf(
                          "Replicate failed (mech=%s rate=%g seed=%d): %s",
                          mech, rate, s, conditionMessage(e)))
                        NULL
                      })
      if (is.null(rep)) next
      for (m_name in names(rep)) {
        acc_all[[length(acc_all) + 1]] <- rep[[m_name]]$acc
        plu_all[[length(plu_all) + 1]] <- rep[[m_name]]$plu
        jpl_all[[length(jpl_all) + 1]] <- rep[[m_name]]$jpl
        wd_all[[length(wd_all)  + 1]] <- rep[[m_name]]$wd
        dwn_all[[length(dwn_all)+ 1]] <- rep[[m_name]]$dwn
      }
    }
  }
  close(pb)
  list(
    accuracy      = do.call(rbind, acc_all),
    plausibility  = do.call(rbind, plu_all),
    joint_plaus   = do.call(rbind, jpl_all),
    wasserstein   = do.call(rbind, wd_all),
    downstream    = do.call(rbind, dwn_all)
  )
}

# -----------------------------------------------------------------------------
# 7. Summarisation: bootstrap CIs and paired Wilcoxon tests
# -----------------------------------------------------------------------------
summarise_bootstrap <- function(df, value_col, group_cols, B = 1000) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      Mean   = mean(.data[[value_col]], na.rm = TRUE),
      Median = median(.data[[value_col]], na.rm = TRUE),
      CI_low  = quantile(replicate(
                  B, mean(sample(.data[[value_col]],
                                 replace = TRUE), na.rm = TRUE)),
                  0.025, na.rm = TRUE),
      CI_high = quantile(replicate(
                  B, mean(sample(.data[[value_col]],
                                 replace = TRUE), na.rm = TRUE)),
                  0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

paired_tests <- function(df, value_col, ref_method = "RF_chained") {
  # Pairs replicates by (Mechanism, Rate, Seed, Variable if present) and
  # tests each non-reference method against the reference with a paired
  # Wilcoxon signed-rank test.
  has_var <- "Variable" %in% names(df)
  key <- c("Mechanism", "Rate", "Seed", if (has_var) "Variable")
  wide <- df %>%
    select(all_of(c(key, "Method", value_col))) %>%
    pivot_wider(names_from = Method, values_from = value_col)
  other <- setdiff(unique(df$Method), ref_method)
  out <- lapply(other, function(m) {
    a <- wide[[ref_method]]; b <- wide[[m]]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5) return(NULL)
    tt <- wilcox.test(a[ok], b[ok], paired = TRUE)
    data.frame(Reference = ref_method, Comparator = m,
               N_pairs = sum(ok),
               Median_diff = median(a[ok] - b[ok]),
               p_value = tt$p.value)
  })
  do.call(rbind, out)
}

# -----------------------------------------------------------------------------
# 8. Main
# -----------------------------------------------------------------------------
main <- function() {
  complete_data <- load_bcw()
  res <- run_experiment(complete_data)

  saveRDS(res, "imputation_results_raw.rds-BC")

  cat("\n===== Accuracy (NRMSE) summary =====\n")
  print(summarise_bootstrap(res$accuracy, "NRMSE",
                            c("Method", "Mechanism", "Rate", "Variable")))

  cat("\n===== Plausibility (% bound violations) =====\n")
  print(summarise_bootstrap(res$plausibility, "ViolationPct",
                            c("Method", "Mechanism", "Rate", "Variable")))

  cat("\n===== Joint plausibility (% outside (Glucose,Insulin) hull) =====\n")
  print(summarise_bootstrap(res$joint_plaus, "JointOutsidePct",
                            c("Method", "Mechanism", "Rate")))

  cat("\n===== Downstream AUC =====\n")
  print(summarise_bootstrap(res$downstream, "AUC",
                            c("Method", "Mechanism", "Rate")))

  cat("\n===== Downstream calibration slope (1.0 is ideal) =====\n")
  print(summarise_bootstrap(res$downstream, "CalSlope",
                            c("Method", "Mechanism", "Rate")))

  cat("\n===== Downstream ECE =====\n")
  print(summarise_bootstrap(res$downstream, "ECE",
                            c("Method", "Mechanism", "Rate")))

  cat("\n===== Paired Wilcoxon: NRMSE (RF_chained vs others) =====\n")
  print(paired_tests(res$accuracy, "NRMSE", "RF_chained"))

  cat("\n===== Paired Wilcoxon: Calibration slope =====\n")
  print(paired_tests(res$downstream, "CalSlope", "RF_chained"))

  cat("\n===== Paired Wilcoxon: ECE =====\n")
  print(paired_tests(res$downstream, "ECE", "RF_chained"))

  cat("\nDone. Raw results saved to imputation_results_raw.rds-BC\n")
  invisible(res)
}

if (sys.nframe() == 0) {
  main()
}
