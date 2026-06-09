# =============================================================================
#  joint_plausibility_mahalanobis.R
#  Adds a SECOND, distribution-based definition of joint implausibility to
#  confirm the convex-hull result is not an artifact of the hull choice.
#
#  Two metrics, same coupled pair as the main paper:
#    1. Hull-based (as before): % of imputed pairs outside the convex hull
#       of observed pairs.  [reproduces the main-paper number]
#    2. Mahalanobis-based: % of imputed pairs whose Mahalanobis distance from
#       the observed pair distribution exceeds the 97.5% chi-square(2) cutoff.
#       This is a smooth, density-aware notion of "outside the joint support",
#       robust to the hard geometric boundary of the hull.
#
#  If MICE-PMM is the worst on BOTH metrics, the joint-implausibility finding
#  is definition-robust. Run on each dataset by editing the DATASET block.
# =============================================================================
options(warn = -1)
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(mice); library(geometry)
})

# ---- DATASET BLOCK: edit these three for pima / bcw / heart ----------------
DATASET <- "pima"   # "pima" | "bcw" | "heart"
cfgs <- list(
  pima = list(path="diabetes.csv", zero_na=TRUE,
              na_vars=c("Glucose","BloodPressure","SkinThickness","Insulin","BMI"),
              outcome="Outcome", pair=c("Glucose","Insulin"),
              cat=character(0)),
  bcw  = list(path="Final-breast_cancer_data-2.csv", zero_na=FALSE,
              na_vars=c("texture_mean","smoothness_mean","concavity_se",
                        "concave.points_se","symmetry_se","fractal_dimension_se",
                        "area_worst","symmetry_worst"),
              outcome="diagnosis", pair=c("concavity_se","concave.points_se"),
              cat=character(0)),
  heart= list(path="heart_disease-1.csv", zero_na=FALSE,
              na_vars=c("chol","thalach","oldpeak","cp","thal","ca"),
              outcome="num", pair=c("thalach","age"),
              cat=c("sex","cp","fbs","restecg","exang","slope","ca","thal"))
)
cfg <- cfgs[[DATASET]]
# ---------------------------------------------------------------------------

N_SEEDS <- 50; RATE <- 0.20  # representative single condition (MCAR, 20%)

load_data <- function(cfg){
  df <- read.csv(cfg$path, check.names=FALSE)
  if (cfg$zero_na) for (v in cfg$na_vars) df[[v]][df[[v]]==0] <- NA
  if (cfg$outcome=="num") df$num <- ifelse(df$num>0,1L,0L)
  for (v in cfg$cat) df[[v]] <- as.factor(df[[v]])
  df[complete.cases(df),]
}

is_cat <- function(v) v %in% cfg$cat

# Minimal chained imputers (continuous-only branch is enough since both pair
# members are continuous in every dataset; categorical predictors are handled
# by RF natively and one-hot for LASSO).
design <- function(d,preds){
  model.matrix(as.formula(paste("~",paste(preds,collapse="+"),"-1")),data=d)
}
impute <- function(dmiss, method){
  d <- dmiss
  miss <- lapply(cfg$na_vars, function(v) which(is.na(d[[v]]))); names(miss)<-cfg$na_vars
  for(v in cfg$na_vars){
    if(is_cat(v)){lev<-names(sort(table(d[[v]]),decreasing=TRUE))[1]
      d[[v]][miss[[v]]]<-lev; d[[v]]<-factor(d[[v]],levels=levels(dmiss[[v]]))}
    else d[[v]][miss[[v]]]<-median(d[[v]],na.rm=TRUE)
  }
  if(method=="MICE"){
    invisible(capture.output(im<-mice(dmiss,m=1,maxit=10,printFlag=FALSE)))
    return(complete(im,1))
  }
  for(it in 1:8) for(v in cfg$na_vars){
    if(length(miss[[v]])==0) next
    preds<-setdiff(names(d),v); tr<-setdiff(seq_len(nrow(d)),miss[[v]])
    if(method=="RF"){
      yv<-d[[v]][tr]; if(is_cat(v)) yv<-droplevels(yv)
      if(is_cat(v)&&nlevels(yv)<2) next
      fit<-tryCatch(randomForest(x=d[tr,preds,drop=FALSE],y=yv,ntree=100),error=function(e)NULL)
      if(is.null(fit)) next
      pr<-predict(fit,newdata=d[miss[[v]],preds,drop=FALSE])
      d[[v]][miss[[v]]]<-if(is_cat(v)) factor(as.character(pr),levels=levels(d[[v]])) else as.numeric(pr)
    } else { # LASSO
      x_tr<-design(d[tr,,drop=FALSE],preds); x_pr<-design(d[miss[[v]],,drop=FALSE],preds)
      cm<-intersect(colnames(x_tr),colnames(x_pr)); x_tr<-x_tr[,cm,drop=FALSE]; x_pr<-x_pr[,cm,drop=FALSE]
      if(is_cat(v)){
        y<-droplevels(d[[v]][tr]); if(nlevels(y)<2) next
        fam<-if(nlevels(y)==2)"binomial" else "multinomial"
        fit<-tryCatch(cv.glmnet(x_tr,y,family=fam,alpha=1),error=function(e)NULL); if(is.null(fit)) next
        pr<-predict(fit,newx=x_pr,s="lambda.1se",type="class")
        d[[v]][miss[[v]]]<-factor(as.character(pr),levels=levels(d[[v]]))
      } else {
        fit<-tryCatch(cv.glmnet(x_tr,d[[v]][tr],family="gaussian",alpha=1),error=function(e)NULL); if(is.null(fit)) next
        d[[v]][miss[[v]]]<-as.numeric(predict(fit,newx=x_pr,s="lambda.1se"))
      }
    }
  }
  d
}

hull_pct <- function(obs, imp){
  res<-tryCatch({h<-convhulln(obs,options="Tv"); mean(!inhulln(h,imp))*100},
                error=function(e) NA_real_); res
}
maha_pct <- function(obs, imp){
  mu<-colMeans(obs); S<-cov(obs)
  res<-tryCatch({
    d2<-mahalanobis(imp, mu, S)
    mean(d2 > qchisq(0.975, df=2))*100
  }, error=function(e) NA_real_)
  res
}

full <- load_data(cfg); v1<-cfg$pair[1]; v2<-cfg$pair[2]
out <- data.frame()
for(method in c("LASSO","RF","MICE")) for(s in 1:N_SEEDS){
  set.seed(s); n<-nrow(full); dmiss<-full
  for(v in cfg$na_vars){mi<-sample.int(n,floor(RATE*n)); dmiss[[v]][mi]<-NA}
  mask <- is.na(dmiss[, cfg$na_vars, drop=FALSE])
  imp <- tryCatch(impute(dmiss, method), error=function(e) NULL); if(is.null(imp)) next
  # rows where either pair member imputed; both pair members must be numeric cols
  m1 <- if(v1 %in% cfg$na_vars) is.na(dmiss[[v1]]) else rep(FALSE,n)
  m2 <- if(v2 %in% cfg$na_vars) is.na(dmiss[[v2]]) else rep(FALSE,n)
  obs_rows<-which(!m1 & !m2); imp_rows<-which(m1 | m2)
  if(length(obs_rows)<5 || length(imp_rows)==0) next
  obs<-as.matrix(full[obs_rows,c(v1,v2)]); imp_pts<-as.matrix(imp[imp_rows,c(v1,v2)])
  out<-rbind(out,data.frame(method=method, seed=s,
                            hull=hull_pct(obs,imp_pts), maha=maha_pct(obs,imp_pts)))
}
cat(sprintf("\n=== %s: joint implausibility, two definitions (median over %d seeds) ===\n",
            toupper(DATASET), N_SEEDS))
agg<-aggregate(cbind(hull,maha)~method, out, function(x) round(median(x,na.rm=TRUE),2))
print(agg)
write.csv(out, paste0("joint_two_metrics_",DATASET,".csv"), row.names=FALSE)
cat("\nIf MICE is worst on BOTH hull and maha columns, the finding is",
    "definition-robust.\n")
