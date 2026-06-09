options(warn=-1)
suppressPackageStartupMessages({library(mice); library(pROC)})
cal_slope <- function(p,y){p<-pmin(pmax(p,1e-6),1-1e-6);unname(coef(suppressWarnings(glm(y~qlogis(p),family=binomial())))[2])}
auc_val <- function(p,y) as.numeric(pROC::auc(pROC::roc(y,p,quiet=TRUE)))
mdown <- function(tm,tt,oc,m){invisible(capture.output(imp<-mice(tm,m=m,method="pmm",maxit=10,printFlag=FALSE)))
  f<-as.formula(paste(oc,"~.")); P<-matrix(NA,nrow(tt),m)
  for(k in 1:m){fit<-suppressWarnings(glm(f,data=complete(imp,k),family=binomial()));P[,k]<-suppressWarnings(predict(fit,newdata=tt,type="response"))}
  rowMeans(P)}
df<-read.csv("diabetes.csv");for(v in c("Glucose","BloodPressure","SkinThickness","Insulin","BMI"))df[[v]][df[[v]]==0]<-NA
full<-df[complete.cases(df),]; nav<-c("Glucose","BloodPressure","SkinThickness","Insulin","BMI")
res<-data.frame()
for(mech in c("MCAR","MAR","MNAR"))for(rate in c(0.1,0.2,0.3))for(s in 1:10){
  set.seed(s);n<-nrow(full);miss<-full
  for(v in nav){
    if(mech=="MCAR")idx<-sample.int(n,floor(rate*n))
    else if(mech=="MAR"){lp<-scale(full$Age)[,1]+scale(full$Pregnancies)[,1];idx<-sample.int(n,floor(rate*n),prob=plogis(lp-mean(lp)))}
    else{idx<-sample.int(n,floor(rate*n),prob=plogis(scale(full[[v]])[,1]))}
    miss[[v]][idx]<-NA}
  set.seed(s+1);tr<-sample.int(n,floor(0.7*n));te<-setdiff(seq_len(n),tr)
  yt<-full[te,"Outcome"];if(length(unique(yt))<2)next
  tryCatch({p1<-mdown(miss[tr,],full[te,],"Outcome",1);p20<-mdown(miss[tr,],full[te,],"Outcome",20)
    res<-rbind(res,data.frame(mech,rate,s,method=c("m1","m20"),
      CalSlope=c(cal_slope(p1,yt),cal_slope(p20,yt)),AUC=c(auc_val(p1,yt),auc_val(p20,yt))))},error=function(e){})
}
write.csv(res,"pima_m_results.csv",row.names=FALSE)
cat("PIMA median CalSlope by m:\n");print(aggregate(CalSlope~method,res,median))
cat("PIMA median AUC by m:\n");print(aggregate(AUC~method,res,median))
w<-reshape(res[,c("mech","rate","s","method","CalSlope")],idvar=c("mech","rate","s"),timevar="method",direction="wide")
d<-w$CalSlope.m20-w$CalSlope.m1;d<-d[is.finite(d)]
cat(sprintf("Paired m20-m1 CalSlope: median=%+.4f p=%.3g n=%d\n",median(d),suppressWarnings(wilcox.test(d))$p.value,length(d)))
