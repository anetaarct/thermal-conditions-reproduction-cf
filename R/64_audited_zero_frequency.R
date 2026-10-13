Sys.setlocale('LC_CTYPE','English_United States.utf8')
.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages(library(glmmTMB))
out<-'manuscript_audit_2026-10-13/zero_frequency_audit'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
models<-list(reference_year_precision=readRDS('models/glmmTMB/fledging_bb_reference_year_precision.rds'))
summaries<-list()
summ<-function(z,obs,model,kind,expected=NA_real_){data.frame(model,simulation=kind,n_sim=length(z),observed=obs,expected=expected,median=median(z),lower=quantile(z,.025),upper=quantile(z,.975),n_ge_observed=sum(z>=obs),upper_tail_plus_one=(sum(z>=obs)+1)/(length(z)+1))}
for(label in names(models)){
 m<-models[[label]];d<-m$frame;n<-rowSums(d[[1]]);obs<-mean(d[[1]][,1]==0);stopifnot(length(n)==9790,length(unique(d$year_f))==39)
 mu<-as.numeric(predict(m,type='response'));phi<-as.numeric(predict(m,type='disp'))
 p0<-exp(lbeta(mu*phi,(1-mu)*phi+n)-lbeta(mu*phi,(1-mu)*phi))
 set.seed(20261009)
 z<-replicate(1000,{p<-rbeta(length(n),mu*phi,(1-mu)*phi);mean(rbinom(length(n),n,p)==0)})
 summaries[[paste0(label,'_conditional')]]<-summ(z,obs,label,'conditional_fitted_modes',mean(p0))
 write.csv(data.frame(rep=seq_along(z),zero=z),file.path(out,paste0(label,'_conditional.csv')),row.names=FALSE)
 # Independent population-level random effects, shared by records of each group.
 X<-model.matrix(~rel_LD*age_class+temp_mean_c+year_c,d);b<-fixef(m)$cond;eta<-drop(X[,names(b),drop=FALSE]%*%b)
 Z<-model.matrix(~rel_LD+temp_mean_c+year_c+age_class,d);g<-fixef(m)$disp;eta_phi<-drop(Z[,names(g),drop=FALSE]%*%g)
 iy<-as.integer(d$year_f);iff<-as.integer(d$female_id);vc<-VarCorr(m)
 sy<-attr(vc$cond$year_f,'stddev')[[1]];sf<-attr(vc$cond$female_id,'stddev')[[1]]
 sp<-attr(vc$disp$year_f,'stddev')[[1]]
 set.seed(20261010)
 z<-replicate(1000,{u<-rnorm(nlevels(d$year_f),0,sy);v<-rnorm(nlevels(d$female_id),0,sf);a<-rnorm(nlevels(d$year_f),0,sp);mup<-plogis(eta+u[iy]+v[iff]);pp<-exp(eta_phi+a[iy]);p<-rbeta(length(n),mup*pp,(1-mup)*pp);mean(rbinom(length(n),n,p)==0)})
 summaries[[paste0(label,'_unconditional')]]<-summ(z,obs,label,'new_group_random_effects')
 write.csv(data.frame(rep=seq_along(z),zero=z),file.path(out,paste0(label,'_unconditional.csv')),row.names=FALSE)
 # Validate the response shape and the erroneous flattening on actual simulate() output.
 ss<-simulate(m,nsim=20,seed=20261011)
 check<-do.call(rbind,lapply(seq_along(ss),function(i){a<-ss[[i]];stopifnot(is.matrix(a),nrow(a)==length(n),ncol(a)==2,all(rowSums(a)==n));data.frame(rep=i,correct=mean(a[,1]==0),incorrect_flattened=mean(as.numeric(a)==0))}))
 write.csv(check,file.path(out,paste0(label,'_simulate_shape_check.csv')),row.names=FALSE)
 annual<-aggregate(cbind(observed=as.numeric(d[[1]][,1]==0),predicted=p0,temperature=d$temp_mean_c),list(year=as.character(d$year_f)),mean)
 counts<-table(d$year_f);annual$n<-as.numeric(counts[annual$year]);annual$residual<-annual$observed-annual$predicted
 write.csv(annual,file.path(out,paste0(label,'_annual_zero_calibration.csv')),row.names=FALSE)
}
tab<-do.call(rbind,summaries);write.csv(tab,file.path(out,'simulation_summary.csv'),row.names=FALSE);print(tab)
# A distributional sensitivity analysis, not an automatic replacement of the main model.
m<-models$year_precision;d<-m$frame;d$fledged<-d[[1]][,1];d$failed<-d[[1]][,2]
zi_file<-file.path(out,'year_precision_zero_inflated_intercept.rds')
if(file.exists(zi_file))zi<-readRDS(zi_file)else{
 zi<-glmmTMB(cbind(fledged,failed)~rel_LD*age_class+temp_mean_c+year_c+(1|year_f)+(1|female_id),
  dispformula=~rel_LD+temp_mean_c+year_c+age_class+(1|year_f),ziformula=~1,family=betabinomial(),data=d,
  control=glmmTMBControl(optCtrl=list(iter.max=20000,eval.max=20000),parallel=4))
 saveRDS(zi,zi_file)
}
diag<-data.frame(convergence=zi$fit$convergence,pdHess=zi$sdr$pdHess,gradient=max(abs(zi$sdr$gradient.fixed)),AIC=AIC(zi),logLik=as.numeric(logLik(zi)))
write.csv(diag,file.path(out,'zero_inflated_diagnostics.csv'),row.names=FALSE);print(diag)
if(zi$fit$convergence==0&&isTRUE(zi$sdr$pdHess)){
 a<-summary(zi)$coefficients$disp['temp_mean_c',];co<-data.frame(beta=a[1],SE=a[2],low=a[1]-1.96*a[2],high=a[1]+1.96*a[2],p=a[4],structural_zero_probability=plogis(fixef(zi)$zi[1]))
 write.csv(co,file.path(out,'zero_inflated_temperature_phi.csv'),row.names=FALSE);print(co)
}
