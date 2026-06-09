# =============================================================================
#  Analysis_v3_heart_COMPLETE.R
#  Imputation comparison on the UCI Cleveland Heart Disease dataset.
#  Self-contained: type-aware imputers + full Monte Carlo driver + all metrics.
#
#  Mixed continuous/categorical handling:
#    Continuous: age, trestbps, chol, thalach, oldpeak
#    Categorical: sex, cp, fbs, restecg, exang, slope, ca, thal
#  Imputers branch on variable type. MICE assigns pmm/logreg/polyreg per column.
#
#  TO RUN: set CFG$data_path to your Cleveland CSV (must contain the raw `num`
#  column). Then:  Rscript Analysis_v3_heart_COMPLETE.R
#  Smoke test first with n_seeds = 5; then set n_seeds = 100 for the full run.
# =============================================================================
suppressPackageStartupMessages({
  pkgs <- c("glmnet","randomForest","mice","pROC","geometry")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  for (p in pkgs) library(p, character.only = TRUE)
  have_transport <- requireNamespace("transport", quietly = TRUE)
  if (have_transport) library(transport)
})

# -----------------------------------------------------------------------------
# 0. Configuration -- EDIT data_path
# -----------------------------------------------------------------------------
CFG <- list(
  data_path = "heart_cleveland.csv",
  continuous_vars  = c("age","trestbps","chol","thalach","oldpeak"),
  categorical_vars = c("sex","cp","fbs","restecg","exang","slope","ca","thal"),
  outcome_var = "num",
  na_vars = c("chol","thalach","oldpeak","cp","thal","ca"),
  miss_rates = c(0.10, 0.20, 0.30),
  mechanisms = c("MCAR","MAR","MNAR"),
  n_seeds = 5,                 # set to 100 for the full run
  max_iter_chain = 10,
  conv_tol = 1e-3,
  rf_ntree = 100,
  dca_thresholds = seq(0.05, 0.50, by = 0.05),
  mar_predictors = c("age","sex"),
  joint_pair = c("thalach","age"),       # most-correlated continuous pair
  out_of_support_quantile = c(0.005, 0.995),
  mice_m = 1                   # set to 20 for multiple-imputation robustness run
)

is_categorical <- function(v) v %in% CFG$categorical_vars

# -----------------------------------------------------------------------------
# 1. Load & prepare
# -----------------------------------------------------------------------------
load_heart <- function(path = CFG$data_path) {
  df <- read.csv(path, check.names = FALSE)
  df[[CFG$outcome_var]] <- ifelse(df[[CFG$outcome_var]] > 0, 1L, 0L)
  for (v in CFG$categorical_vars) df[[v]] <- as.factor(df[[v]])
  df[[CFG$outcome_var]] <- as.integer(df[[CFG$outcome_var]])
  cd <- df[complete.cases(df), ]
  cat(sprintf("Loaded %d rows, %d complete cases.\n", nrow(df), nrow(cd)))
  cd
}

# Build data-driven out-of-support bounds for continuous vars (once).
make_bounds <- function(full) {
  b <- list()
  for (v in CFG$continuous_vars) {
    q <- quantile(full[[v]], CFG$out_of_support_quantile, na.rm = TRUE)
    b[[v]] <- c(max(0, q[1]), q[2])
  }
  b
}

# One-hot design matrix for glmnet (factors -> dummies).
design_matrix <- function(df, predictors) {
  fmla <- as.formula(paste("~", paste(predictors, collapse = " + "), "- 1"))
  model.matrix(fmla, data = df)
}

# -----------------------------------------------------------------------------
# 2. Missingness simulation
# -----------------------------------------------------------------------------
simulate_missing <- function(data, vars, rate, mech, seed) {
  set.seed(seed); n <- nrow(data); out <- data
  mask <- matrix(FALSE, n, length(vars), dimnames = list(NULL, vars))
  for (v in vars) {
    nm <- floor(rate * n)
    if (mech == "MCAR") idx <- sample.int(n, nm)
    else if (mech == "MAR") {
      lp <- scale(as.numeric(data[[CFG$mar_predictors[1]]]))[,1] +
            scale(as.numeric(data[[CFG$mar_predictors[2]]]))[,1]
      idx <- sample.int(n, nm, prob = plogis(lp - mean(lp)))
    } else {
      x <- as.numeric(data[[v]])
      idx <- sample.int(n, nm, prob = plogis(scale(x)[,1]))
    }
    out[[v]][idx] <- NA; mask[idx, v] <- TRUE
  }
  list(data = out, mask = mask)
}

# -----------------------------------------------------------------------------
# 3a. Chained LASSO (type-aware)
# -----------------------------------------------------------------------------
impute_lasso_chained <- function(df_miss, target_vars,
                                 max_iter = CFG$max_iter_chain) {
  df <- df_miss
  miss <- lapply(target_vars, function(v) which(is.na(df[[v]]))); names(miss) <- target_vars
  for (v in target_vars) {
    if (is_categorical(v)) {
      lev <- names(sort(table(df[[v]]), decreasing = TRUE))[1]
      df[[v]][miss[[v]]] <- lev
      df[[v]] <- factor(df[[v]], levels = levels(df_miss[[v]]))
    } else df[[v]][miss[[v]]] <- median(df[[v]], na.rm = TRUE)
  }
  for (it in seq_len(max_iter)) {
    for (v in target_vars) {
      if (length(miss[[v]]) == 0) next
      preds <- setdiff(names(df), v)
      train <- setdiff(seq_len(nrow(df)), miss[[v]])
      x_tr <- design_matrix(df[train, , drop=FALSE], preds)
      x_pr <- design_matrix(df[miss[[v]], , drop=FALSE], preds)
      common <- intersect(colnames(x_tr), colnames(x_pr))
      if (length(common) < 1) next
      x_tr <- x_tr[, common, drop=FALSE]; x_pr <- x_pr[, common, drop=FALSE]
      if (is_categorical(v)) {
        y <- droplevels(df[[v]][train]); if (nlevels(y) < 2) next
        fam <- if (nlevels(y) == 2) "binomial" else "multinomial"
        fit <- tryCatch(cv.glmnet(x_tr, y, family = fam, alpha = 1), error=function(e) NULL)
        if (is.null(fit)) next
        pr <- predict(fit, newx = x_pr, s = "lambda.1se", type = "class")
        df[[v]][miss[[v]]] <- factor(as.character(pr), levels = levels(df[[v]]))
      } else {
        y <- df[[v]][train]
        fit <- tryCatch(cv.glmnet(x_tr, y, family = "gaussian", alpha = 1), error=function(e) NULL)
        if (is.null(fit)) next
        df[[v]][miss[[v]]] <- as.numeric(predict(fit, newx = x_pr, s = "lambda.1se"))
      }
    }
  }
  df
}

# -----------------------------------------------------------------------------
# 3b. Chained RF (type-aware; RF handles factors natively)
# -----------------------------------------------------------------------------
impute_rf_chained <- function(df_miss, target_vars,
                              max_iter = CFG$max_iter_chain, ntree = CFG$rf_ntree) {
  df <- df_miss
  miss <- lapply(target_vars, function(v) which(is.na(df[[v]]))); names(miss) <- target_vars
  for (v in target_vars) {
    if (is_categorical(v)) {
      lev <- names(sort(table(df[[v]]), decreasing = TRUE))[1]
      df[[v]][miss[[v]]] <- lev
      df[[v]] <- factor(df[[v]], levels = levels(df_miss[[v]]))
    } else df[[v]][miss[[v]]] <- median(df[[v]], na.rm = TRUE)
  }
  for (it in seq_len(max_iter)) {
    for (v in target_vars) {
      if (length(miss[[v]]) == 0) next
      preds <- setdiff(names(df), v)
      train <- setdiff(seq_len(nrow(df)), miss[[v]])
      yv <- df[[v]][train]; if (is_categorical(v)) yv <- droplevels(yv)
      if (is_categorical(v) && nlevels(yv) < 2) next
      fit <- tryCatch(randomForest(x = df[train, preds, drop=FALSE], y = yv, ntree = ntree),
                      error = function(e) NULL)
      if (is.null(fit)) next
      pr <- predict(fit, newdata = df[miss[[v]], preds, drop=FALSE])
      if (is_categorical(v)) df[[v]][miss[[v]]] <- factor(as.character(pr), levels=levels(df[[v]]))
      else df[[v]][miss[[v]]] <- as.numeric(pr)
    }
  }
  df
}

# -----------------------------------------------------------------------------
# 3c. MICE (per-type methods chosen automatically)
# -----------------------------------------------------------------------------
impute_mice <- function(df_miss, m = CFG$mice_m) {
  invisible(capture.output(imp <- mice(df_miss, m = m, maxit = 10, printFlag = FALSE)))
  imp
}

# -----------------------------------------------------------------------------
# 4. Metrics
# -----------------------------------------------------------------------------
accuracy_metrics <- function(imp, truth, mask, vars) {
  rows <- list()
  for (v in vars) {
    idx <- which(mask[,v]); if (length(idx)==0) next
    if (is_categorical(v)) {
      pfc <- mean(as.character(imp[[v]][idx]) != as.character(truth[[v]][idx]))
      rows[[v]] <- data.frame(Variable=v, Type="categorical", NRMSE=NA, PFC=pfc)
    } else {
      t <- truth[[v]][idx]; i <- as.numeric(imp[[v]][idx])
      ok <- is.finite(t)&is.finite(i); t<-t[ok]; i<-i[ok]
      rows[[v]] <- data.frame(Variable=v, Type="continuous",
                              NRMSE=sqrt(mean((t-i)^2))/sd(t), PFC=NA)
    }
  }
  do.call(rbind, rows)
}

plausibility_univariate <- function(imputed, truth, mask, vars, bounds) {
  rows <- list()
  for (v in vars) {
    idx <- which(mask[, v]); if (length(idx) == 0) next
    if (is_categorical(v)) {
      obs_levels <- unique(as.character(truth[[v]]))
      x <- as.character(imputed[[v]][idx])
      viol <- mean(!(x %in% obs_levels)) * 100
      rows[[v]] <- data.frame(Variable=v, Type="categorical",
                              ViolationPct=viol, NegativePct=NA_real_)
    } else {
      x <- as.numeric(imputed[[v]][idx]); b <- bounds[[v]]
      viol <- mean(x < b[1] | x > b[2]) * 100
      neg  <- mean(x < 0) * 100
      rows[[v]] <- data.frame(Variable=v, Type="continuous",
                              ViolationPct=viol, NegativePct=neg)
    }
  }
  do.call(rbind, rows)
}

plausibility_joint <- function(imputed, truth, mask,
                               v1 = CFG$joint_pair[1], v2 = CFG$joint_pair[2]) {
  if (!(v1 %in% colnames(truth)) || !(v2 %in% colnames(truth)))
    return(data.frame(JointOutsidePct = NA_real_))
  m1 <- if (v1 %in% colnames(mask)) mask[, v1] else rep(FALSE, nrow(truth))
  m2 <- if (v2 %in% colnames(mask)) mask[, v2] else rep(FALSE, nrow(truth))
  obs_rows <- which(!m1 & !m2); imp_rows <- which(m1 | m2)
  if (length(obs_rows) < 4 || length(imp_rows) == 0)
    return(data.frame(JointOutsidePct = NA_real_))
  obs_pts <- as.matrix(truth[obs_rows, c(v1, v2)])
  imp_pts <- as.matrix(imputed[imp_rows, c(v1, v2)])
  res_pct <- tryCatch({
    hull <- geometry::convhulln(obs_pts, options = "Tv")
    inside <- geometry::inhulln(hull, imp_pts)
    mean(!inside) * 100
  }, error = function(e) NA_real_)
  data.frame(JointOutsidePct = res_pct)
}

wasserstein_fidelity <- function(imputed, truth, vars) {
  rows <- list()
  for (v in vars) {
    if (is_categorical(v)) next  # Wasserstein-1 only meaningful for continuous
    a <- as.numeric(truth[[v]]); b <- as.numeric(imputed[[v]])
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a)==0 || length(b)==0) next
    w <- if (have_transport)
           tryCatch(transport::wasserstein1d(a,b), error=function(e) NA_real_)
         else NA_real_
    rows[[v]] <- data.frame(Variable=v, Wasserstein1=w)
  }
  if (length(rows)==0) return(data.frame(Variable=character(), Wasserstein1=numeric()))
  do.call(rbind, rows)
}

downstream_from_pred <- function(p, y, thr = CFG$dca_thresholds) {
  p <- pmin(pmax(p,1e-6),1-1e-6)
  auc <- as.numeric(pROC::auc(pROC::roc(y, p, quiet=TRUE)))
  brier <- mean((p-y)^2)
  cal <- suppressWarnings(glm(y ~ qlogis(p), family=binomial()))
  bins <- cut(p, seq(0,1,length.out=11), include.lowest=TRUE)
  ece <- sum(tapply(seq_along(p), bins, function(ii)
    if (length(ii)==0) 0 else abs(mean(p[ii])-mean(y[ii]))*length(ii)/length(p)), na.rm=TRUE)
  nb <- mean(sapply(thr, function(t){ pp <- p>=t; if (sum(pp)==0) return(0)
    sum(pp & y==1)/length(y) - sum(pp & y==0)/length(y)*(t/(1-t)) }))
  data.frame(AUC=auc, Brier=brier, CalSlope=unname(coef(cal)[2]),
             CalIntercept=unname(coef(cal)[1]), ECE=ece, NetBenefit=nb)
}

downstream_df <- function(train_df, test_df, outcome) {
  fmla <- as.formula(paste(outcome, "~ ."))
  fit <- suppressWarnings(glm(fmla, data=train_df, family=binomial()))
  p <- suppressWarnings(predict(fit, newdata=test_df, type="response"))
  y <- test_df[[outcome]]; ok <- is.finite(p)&is.finite(y)
  downstream_from_pred(p[ok], y[ok])
}

downstream_mice <- function(imp, test_df, outcome, m) {
  fmla <- as.formula(paste(outcome, "~ ."))
  P <- matrix(NA, nrow(test_df), m)
  for (k in seq_len(m)) {
    fit <- suppressWarnings(glm(fmla, data=complete(imp,k), family=binomial()))
    P[,k] <- suppressWarnings(predict(fit, newdata=test_df, type="response"))
  }
  y <- test_df[[outcome]]; p <- rowMeans(P); ok <- is.finite(p)&is.finite(y)
  downstream_from_pred(p[ok], y[ok])
}

# -----------------------------------------------------------------------------
# 5. One replicate
# -----------------------------------------------------------------------------
one_replicate <- function(full, bounds, mech, rate, seed) {
  sim <- simulate_missing(full, CFG$na_vars, rate, mech, seed)
  set.seed(seed + 1); n <- nrow(full)
  tr <- sample.int(n, floor(0.7*n)); te <- setdiff(seq_len(n), tr)
  train_truth <- full[tr,]; test_truth <- full[te,]
  train_miss <- sim$data[tr,]; mask_tr <- sim$mask[tr,,drop=FALSE]
  if (length(unique(test_truth[[CFG$outcome_var]])) < 2) return(NULL)

  imps <- list(
    LASSO_chained = impute_lasso_chained(train_miss, CFG$na_vars),
    RF_chained    = impute_rf_chained(train_miss, CFG$na_vars)
  )
  mice_obj <- impute_mice(train_miss, CFG$mice_m)

  acc<-plu<-jpl<-wd<-dwn<-list()
  for (m in names(imps)) {
    imp <- imps[[m]]
    acc[[m]] <- cbind(Method=m, Mechanism=mech, Rate=rate, Seed=seed,
                      accuracy_metrics(imp, train_truth, mask_tr, CFG$na_vars))
    plu[[m]] <- cbind(Method=m, Mechanism=mech, Rate=rate, Seed=seed,
                      plausibility_univariate(imp, train_truth, mask_tr, CFG$na_vars, bounds))
    jpl[[m]] <- cbind(Method=m, Mechanism=mech, Rate=rate, Seed=seed,
                      plausibility_joint(imp, train_truth, mask_tr))
    wd[[m]]  <- cbind(Method=m, Mechanism=mech, Rate=rate, Seed=seed,
                      wasserstein_fidelity(imp, train_truth, CFG$na_vars))
    dwn[[m]] <- cbind(Method=m, Mechanism=mech, Rate=rate, Seed=seed,
                      downstream_df(imp, test_truth, CFG$outcome_var))
  }
  # MICE: completed dataset 1 for accuracy/plausibility; Rubin pooling downstream
  mc <- complete(mice_obj, 1)
  acc[["MICE_PMM"]] <- cbind(Method="MICE_PMM", Mechanism=mech, Rate=rate, Seed=seed,
                             accuracy_metrics(mc, train_truth, mask_tr, CFG$na_vars))
  plu[["MICE_PMM"]] <- cbind(Method="MICE_PMM", Mechanism=mech, Rate=rate, Seed=seed,
                             plausibility_univariate(mc, train_truth, mask_tr, CFG$na_vars, bounds))
  jpl[["MICE_PMM"]] <- cbind(Method="MICE_PMM", Mechanism=mech, Rate=rate, Seed=seed,
                             plausibility_joint(mc, train_truth, mask_tr))
  wd[["MICE_PMM"]]  <- cbind(Method="MICE_PMM", Mechanism=mech, Rate=rate, Seed=seed,
                             wasserstein_fidelity(mc, train_truth, CFG$na_vars))
  dwn[["MICE_PMM"]] <- cbind(Method="MICE_PMM", Mechanism=mech, Rate=rate, Seed=seed,
                             downstream_mice(mice_obj, test_truth, CFG$outcome_var, CFG$mice_m))

  list(acc=do.call(rbind,acc), plu=do.call(rbind,plu), jpl=do.call(rbind,jpl),
       wd=do.call(rbind,wd), dwn=do.call(rbind,dwn))
}

# -----------------------------------------------------------------------------
# 6. Driver
# -----------------------------------------------------------------------------
run_experiment <- function(full) {
  bounds <- make_bounds(full)
  conds <- expand.grid(mech=CFG$mechanisms, rate=CFG$miss_rates, stringsAsFactors=FALSE)
  A<-P<-J<-W<-D<-list()
  total <- nrow(conds)*CFG$n_seeds; k<-0; pb<-txtProgressBar(0,total,style=3)
  for (ci in seq_len(nrow(conds))) for (s in seq_len(CFG$n_seeds)) {
    k<-k+1; setTxtProgressBar(pb,k)
    r <- tryCatch(one_replicate(full, bounds, conds$mech[ci], conds$rate[ci], s),
                  error=function(e){message(sprintf("\nfail mech=%s rate=%.2f seed=%d: %s",
                    conds$mech[ci],conds$rate[ci],s,conditionMessage(e))); NULL})
    if (is.null(r)) next
    A[[length(A)+1]]<-r$acc; P[[length(P)+1]]<-r$plu; J[[length(J)+1]]<-r$jpl
    W[[length(W)+1]]<-r$wd;  D[[length(D)+1]]<-r$dwn
  }
  close(pb)
  list(accuracy=do.call(rbind,A), plausibility=do.call(rbind,P),
       joint_plaus=do.call(rbind,J), wasserstein=do.call(rbind,W),
       downstream=do.call(rbind,D))
}

# -----------------------------------------------------------------------------
# 7. Main
# -----------------------------------------------------------------------------
main <- function() {
  full <- load_heart()
  res <- run_experiment(full)
  saveRDS(res, "imputation_results_raw.rds-HC")
  cat("\n--- Joint plausibility (median by method) ---\n")
  print(aggregate(JointOutsidePct ~ Method, res$joint_plaus, median, na.rm=TRUE))
  cat("\n--- Univariate violation % (mean by method) ---\n")
  print(aggregate(ViolationPct ~ Method, res$plausibility, mean, na.rm=TRUE))
  cat("\n--- Downstream calibration slope (median by method) ---\n")
  print(aggregate(CalSlope ~ Method, res$downstream, median, na.rm=TRUE))
  cat("\n--- Downstream AUC (median by method) ---\n")
  print(aggregate(AUC ~ Method, res$downstream, median, na.rm=TRUE))
  cat("\nSaved imputation_results_raw.rds-HC\n")
  invisible(res)
}

if (sys.nframe() == 0) main()
