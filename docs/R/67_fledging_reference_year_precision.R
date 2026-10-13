#!/usr/bin/env Rscript
Sys.setlocale('LC_CTYPE','English_United States.utf8')
.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages(library(glmmTMB))
out <- 'manuscript_audit_2026-10-12/new_reference'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ref_path <- 'models/glmmTMB/fledging_bb_reference_year_precision.rds'
audited_path <- 'manuscript_audit_2026-10-08/precision_year_sensitivity/nlminb_refined.rds'
if(file.exists(ref_path))full<-readRDS(ref_path)else if(file.exists(audited_path))full<-readRDS(audited_path)else{
 raw<-read.csv('data_derived/fledging_model_data.csv')
 raw$year_f<-factor(raw$year_f);raw$female_id<-factor(raw$female_id)
 raw$age_class<-factor(raw$age_class,levels=c('YOUNG','OLD'))
 raw$failed<-raw$clutch-raw$fledged
 full<-glmmTMB(cbind(fledged,failed)~rel_LD*age_class+temp_mean_c+year_c+(1|year_f)+(1|female_id),
  dispformula=~rel_LD+temp_mean_c+year_c+age_class+(1|year_f),ziformula=~0,family=betabinomial(),data=raw,
  control=glmmTMBControl(optCtrl=list(iter.max=20000,eval.max=20000),parallel=4))
 stopifnot(full$fit$convergence==0L,isTRUE(full$sdr$pdHess));saveRDS(full,ref_path)
}
d <- full$frame;d$fledged <- d[[1]][,1];d$failed <- d[[1]][,2]
stopifnot(nrow(d)==9790L,length(unique(d$year_f))==39L,
          sum(d$fledged==0)==2315L,full$fit$convergence==0L,
          isTRUE(full$sdr$pdHess))
f <- cbind(fledged,failed)~rel_LD*age_class+temp_mean_c+year_c+(1|year_f)+(1|female_id)
disp0 <- ~1+(1|year_f)
start_from <- function(m,only_intercept=FALSE){
 pp<-m$fit$parfull
 keys<-c('beta','betadisp','theta','thetadisp','b','bdisp')
 ans<-setNames(lapply(keys,function(k)unname(pp[names(pp)==k])),keys)
 if(only_intercept)ans$betadisp<-unname(fixef(m)$disp['(Intercept)'])
 ans[lengths(ans)>0]
}
fits<-list()
for(method in c('nlminb','BFGS')){
 path<-file.path(out,paste0('constant_fixed_precision_',method,'.rds'))
 if(file.exists(path))m<-readRDS(path)else{
  control<-if(method=='nlminb')glmmTMBControl(optCtrl=list(iter.max=20000,eval.max=20000,rel.tol=1e-11,x.tol=1e-10),parallel=4)else glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS'),optCtrl=list(maxit=20000,reltol=1e-12),parallel=4)
  st<-if(method=='nlminb')start_from(full,TRUE)else start_from(fits$nlminb)
  m<-glmmTMB(f,dispformula=disp0,ziformula=~0,family=betabinomial(),data=d,start=st,control=control)
  saveRDS(m,path)
 }
 fits[[method]]<-m
 diag<-data.frame(method,convergence=m$fit$convergence,pdHess=m$sdr$pdHess,gradient=max(abs(m$sdr$gradient.fixed)),logLik=as.numeric(logLik(m)),df=attr(logLik(m),'df'),AIC=AIC(m),BIC=BIC(m))
 write.csv(diag,file.path(out,paste0(method,'_diagnostics.csv')),row.names=FALSE);print(diag)
 stopifnot(isTRUE(m$sdr$pdHess))
}
stopifnot(fits$BFGS$fit$convergence==0L)
refined_path<-file.path(out,'constant_fixed_precision_nlminb_warm.rds')
if(file.exists(refined_path))refined<-readRDS(refined_path)else{
 refined<-glmmTMB(f,dispformula=disp0,ziformula=~0,family=betabinomial(),data=d,
  start=start_from(fits$BFGS),control=glmmTMBControl(optCtrl=list(iter.max=20000,eval.max=20000,rel.tol=1e-10),parallel=4))
 saveRDS(refined,refined_path)
}
rd<-data.frame(method='nlminb_warm',convergence=refined$fit$convergence,pdHess=refined$sdr$pdHess,gradient=max(abs(refined$sdr$gradient.fixed)),logLik=as.numeric(logLik(refined)),df=attr(logLik(refined),'df'),AIC=AIC(refined),BIC=BIC(refined))
write.csv(rd,file.path(out,'nlminb_warm_diagnostics.csv'),row.names=FALSE);print(rd)
stopifnot(refined$fit$convergence==0L,isTRUE(refined$sdr$pdHess),
          abs(as.numeric(logLik(refined)-logLik(fits$BFGS)))<1e-4)
null<-refined
lrt<-anova(null,full)
metrics<-function(m,label)data.frame(model=label,n=nobs(m),df=attr(logLik(m),'df'),logLik=as.numeric(logLik(m)),AIC=AIC(m),BIC=BIC(m))
tab<-rbind(metrics(null,'Constant fixed precision + year RE'),metrics(full,'Predictor-dependent precision + year RE'))
tab$delta_AIC<-tab$AIC-min(tab$AIC);tab$delta_BIC<-tab$BIC-min(tab$BIC)
tab$LR_chisq<-c(NA,lrt$Chisq[2]);tab$LRT_df<-c(NA,lrt$`Chi Df`[2]);tab$LRT_p<-c(NA,lrt$`Pr(>Chisq)`[2])
write.csv(tab,file.path(out,'precision_comparison.csv'),row.names=FALSE);print(tab)
co<-do.call(rbind,lapply(c('cond','disp'),function(k){a<-summary(full)$coefficients[[k]];data.frame(component=if(k=='cond')'location'else'precision_log_phi',term=rownames(a),beta=a[,1],SE=a[,2],CI_low=a[,1]-1.96*a[,2],CI_high=a[,1]+1.96*a[,2],z=a[,3],p=a[,4])}))
write.csv(co,file.path(out,'coefficients.csv'),row.names=FALSE);print(co)
vc<-VarCorr(full)
sdtab<-data.frame(component=c('location','location','precision_log_phi'),group=c('year','female','year'),SD=c(attr(vc$cond$year_f,'stddev')[[1]],attr(vc$cond$female_id,'stddev')[[1]],attr(vc$disp$year_f,'stddev')[[1]]))
write.csv(sdtab,file.path(out,'random_effect_SD.csv'),row.names=FALSE);print(sdtab)
saveRDS(full,file.path(out,'reference_full.rds'))
cat('New reference and matched-random-structure LRT verified.\n')
