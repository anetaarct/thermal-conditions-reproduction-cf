#!/usr/bin/env Rscript
Sys.setlocale('LC_CTYPE','English_United States.utf8')
.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages(library(glmmTMB))
source('R/58_window_sensitivity_data.R')
out<-'manuscript_audit_2026-10-12/new_reference';tt<-read.csv('tables/window_sensitivity_temperatures.csv')
full<-readRDS(file.path(out,'reference_full.rds'));pp<-full$fit$parfull
keys<-c('beta','betadisp','theta','thetadisp','b','bdisp');st<-setNames(lapply(keys,function(k)unname(pp[names(pp)==k])),keys)
results<-list();diagnostics<-list()
for(win in window_specs$window){
 d<-prepare_window_data('data_derived/fledging_model_data.csv',tt,win)
 path<-file.path(out,paste0('window_',win,'_year_precision.rds'))
 if(win=='main')m<-full else if(file.exists(path))m<-readRDS(path)else{
  f<-cbind(fledged,failed)~rel_LD*age_class+temp_mean_c+year_c+(1|year_f)+(1|female_id)
  phi<-~rel_LD+temp_mean_c+year_c+age_class+(1|year_f)
  m<-glmmTMB(f,dispformula=phi,family=betabinomial(),ziformula=~0,data=d,start=st,
   control=glmmTMBControl(optCtrl=list(iter.max=20000,eval.max=20000,rel.tol=1e-9),parallel=4))
  if(m$fit$convergence!=0L||!isTRUE(m$sdr$pdHess)){
   p<-m$fit$parfull;ss<-setNames(lapply(keys,function(k)unname(p[names(p)==k])),keys)
   m<-glmmTMB(f,dispformula=phi,family=betabinomial(),ziformula=~0,data=d,start=ss,
    control=glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS'),optCtrl=list(maxit=20000,reltol=1e-11),parallel=4))
  }
  stopifnot(m$fit$convergence==0L,isTRUE(m$sdr$pdHess));saveRDS(m,path)
 }
 ll_difference<-NA_real_
 if(win!='main'){
  alt_path<-file.path(out,paste0('window_',win,'_BFGS_verified.rds'))
  if(file.exists(alt_path))alt<-readRDS(alt_path)else{
   p<-m$fit$parfull;ss<-setNames(lapply(keys,function(k)unname(p[names(p)==k])),keys)
   alt<-glmmTMB(cbind(fledged,failed)~rel_LD*age_class+temp_mean_c+year_c+(1|year_f)+(1|female_id),
    dispformula=~rel_LD+temp_mean_c+year_c+age_class+(1|year_f),ziformula=~0,family=betabinomial(),data=d,start=ss,
    control=glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS'),optCtrl=list(maxit=20000,reltol=1e-12),parallel=4))
   saveRDS(alt,alt_path)
  }
  ll_difference<-abs(as.numeric(logLik(m)-logLik(alt)))
  stopifnot(alt$fit$convergence==0L,isTRUE(alt$sdr$pdHess),ll_difference<1e-4)
 }
 diagnostics[[win]]<-data.frame(window=win,n=nobs(m),convergence=m$fit$convergence,pdHess=m$sdr$pdHess,max_gradient=max(abs(m$sdr$gradient.fixed)),scaled_gradient=max(abs(m$sdr$gradient.fixed)*sqrt(diag(m$sdr$cov.fixed))),optimizer_logLik_difference=ll_difference,AIC=AIC(m))
 results[[win]]<-do.call(rbind,lapply(c('cond','disp'),function(k){a<-summary(m)$coefficients[[k]]['temp_mean_c',];data.frame(window=win,component=k,estimate=a[1],SE=a[2],lower95=a[1]-1.96*a[2],upper95=a[1]+1.96*a[2],p=a[4],annual_temperature_SD=sd(tt$temperature[tt$window==win&tt$year%in%as.integer(as.character(d$year_f))]))}))
 write.csv(do.call(rbind,results),file.path(out,'window_sensitivity.csv'),row.names=FALSE)
 write.csv(do.call(rbind,diagnostics),file.path(out,'window_diagnostics.csv'),row.names=FALSE)
 print(results[[win]])
}
cat('All calendar sensitivities use year RE in both model components.\n')
