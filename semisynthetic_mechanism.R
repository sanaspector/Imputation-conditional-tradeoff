# =============================================================================
#  semisynthetic_mechanism.R
#  Tests the boundary-condition mechanism by manipulating, on the Pima data
#  (where MICE-PMM's calibration advantage EXISTS), the two factors we
#  hypothesise control that advantage:
#     (A) SAMPLE SIZE  -> calibratability (smaller n => less headroom)
#     (B) DISCRETIZATION -> data-type composition (binning continuous vars
#                           into categories removes the marginal-fidelity
#                           channel MICE-PMM exploits)
#  If MICE-PMM's calibration edge shrinks as we (A) shrink n and (B) discretize,
#  that is near-causal evidence for the boundary, far stronger than 3
#  observational datasets alone.
# =============================================================================
options(warn=-1)
suppressPackageStartupMessages({library(mice); library(pROC)})

set.seed(1)
df <- read.csv("diabetes.csv")
for (v in c("Glucose","BloodPressure","SkinThickness","Insulin","BMI"))
  df[[v]][df[[v]]==0] <- NA
full0 <- df[complete.cases(df),]      # 392 complete cases
NAV <- c("Glucose","BloodPressure","SkinThickness","Insulin","BMI")

cal_slope <- function(p,y){p<-pmin(pmax(p,1e-6),1-1e-6)
  unname(coef(suppressWarnings(glm(y~qlogis(p),family=binomial())))[2])}

# downstream via MICE (m=1) and via simple median impute (proxy for a
# conditional-mean method) so we can measure the GAP between a donor method
# and a non-donor method.
mice_pred <- function(tm, te, oc){
  invisible(capture.output(imp<-mice(tm,m=1,method="pmm",maxit=10,printFlag=FALSE)))
  fit<-suppressWarnings(glm(as.formula(paste(oc,"~.")),data=complete(imp,1),family=binomial()))
  suppressWarnings(predict(fit,newdata=te,type="response"))
}
median_pred <- function(tm, te, oc){
  tmc <- tm
  for(v in NAV) tmc[[v]][is.na(tmc[[v]])] <- median(tmc[[v]],na.rm=TRUE)
  fit<-suppressWarnings(glm(as.formula(paste(oc,"~.")),data=tmc,family=binomial()))
  suppressWarnings(predict(fit,newdata=te,type="response"))
}

discretize <- function(d, vars, k){
  # replace continuous var by k-level factor (k=Inf => unchanged)
  if (is.infinite(k)) return(d)
  out <- d
  for(v in vars){
    br <- unique(quantile(d[[v]], seq(0,1,length.out=k+1), na.rm=TRUE))
    if(length(br)<3){next}
    out[[v]] <- as.numeric(cut(d[[v]], breaks=br, include.lowest=TRUE,
                               labels=FALSE))
  }
  out
}

run_condition <- function(n_target, k_levels, n_seeds=100, rate=0.2){
  gaps <- c()   # MICE cal slope - median cal slope (positive => MICE better)
  mslp <- c(); dslp <- c()
  for(s in 1:n_seeds){
    set.seed(s)
    # (A) subsample to n_target
    idx <- sample(nrow(full0), min(n_target, nrow(full0)))
    d <- full0[idx,]
    # (B) discretize continuous predictors to k levels
    d <- discretize(d, NAV, k_levels)
    n <- nrow(d)
    # impose MCAR missingness
    miss <- d
    for(v in NAV){ mi<-sample.int(n, floor(rate*n)); miss[[v]][mi]<-NA }
    tr <- sample.int(n, floor(0.7*n)); te <- setdiff(seq_len(n),tr)
    yt <- d$Outcome[te]; if(length(unique(yt))<2) next
    pm <- tryCatch(mice_pred(miss[tr,], d[te,], "Outcome"), error=function(e) NULL)
    pd <- tryCatch(median_pred(miss[tr,], d[te,], "Outcome"), error=function(e) NULL)
    if(is.null(pm)||is.null(pd)) next
    cm <- cal_slope(pm,yt); cd <- cal_slope(pd,yt)
    if(!is.finite(cm)||!is.finite(cd)) next
    mslp<-c(mslp,cm); dslp<-c(dslp,cd); gaps<-c(gaps, cm-cd)
  }
  data.frame(n_target=n_target, k_levels=ifelse(is.infinite(k_levels),"cont",k_levels),
             mice_cal=median(mslp,na.rm=TRUE), median_cal=median(dslp,na.rm=TRUE),
             gap=median(gaps,na.rm=TRUE), n_ok=length(gaps))
}

cat("=== (A) Sample-size sweep (continuous, vary n) ===\n")
resA <- do.call(rbind, lapply(c(392,250,150,80), function(n) run_condition(n, Inf)))
print(resA)

cat("\n=== (B) Discretization sweep (full n=392, vary #levels) ===\n")
resB <- do.call(rbind, lapply(c(Inf,5,3,2), function(k) run_condition(392, k)))
print(resB)

write.csv(rbind(resA,resB), "semisynthetic_results.csv", row.names=FALSE)
cat("\nInterpretation: 'gap' is MICE-PMM calibration slope minus a non-donor\n")
cat("(median-impute) baseline. A shrinking gap as n falls or as variables are\n")
cat("discretized is evidence that BOTH factors control the calibration advantage.\n")
