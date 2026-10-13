#!/usr/bin/env Rscript
Sys.setlocale('LC_CTYPE','English_United States.utf8')
.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages({library(glmmTMB);library(ggplot2);library(readr);library(dplyr);library(tibble);library(patchwork)})
out<-'manuscript_audit_2026-10-12/new_reference'
m<-readRDS(file.path(out,'reference_full.rds'))
null<-readRDS(file.path(out,'constant_fixed_precision_nlminb_warm.rds'))
d<-m$frame;d$year<-as.integer(as.character(d$year_f))
stopifnot(nobs(m)==9790L,isTRUE(m$sdr$pdHess),m$fit$convergence==0L)
saveRDS(m,'models/glmmTMB/fledging_bb_reference_year_precision.rds')
saveRDS(null,'models/glmmTMB/fledging_bb_no_fixed_precision_predictors_year_RE.rds')
cmp<-read.csv(file.path(out,'precision_comparison.csv'))
cmp$model<-c('No fixed precision predictors + year RE','Predictor-dependent precision + year RE')
write_csv(cmp,'tables/fledging_bb_precision_comparison.csv')
co<-read.csv(file.path(out,'coefficients.csv'));co$model<-'beta_binomial_reference_year_precision'
write_csv(co,'tables/fledging_bb_main_coefficients.csv')
sdtab<-read.csv(file.path(out,'random_effect_SD.csv'))
var<-data.frame(component=c('Location','Location','Precision'),quantity=c('Year random-intercept SD','Female random-intercept SD','Year random-intercept SD'),estimate=sdtab$SD,unit=c('logit scale','logit scale','log-phi scale'))
write_csv(var,'tables/fledging_variance_decomposition.csv')
write_csv(data.frame(model=c('reference_year_precision','no_fixed_precision_predictors_year_RE'),pdHess=c(m$sdr$pdHess,null$sdr$pdHess),convergence_code=c(m$fit$convergence,null$fit$convergence)),'tables/fledging_chapter_diagnostics.csv')
# Exact covariance-based age-specific slopes in the location component.
b<-fixef(m)$cond;V<-vcov(m)$cond
trends<-do.call(rbind,lapply(c('YOUNG','OLD'),function(age){w<-setNames(rep(0,length(b)),names(b));w['rel_LD']<-1;if(age=='OLD')w['rel_LD:age_classOLD']<-1;est<-sum(w*b);se<-sqrt(drop(t(w)%*%V%*%w));data.frame(age_class=age,rel_LD.trend=est,SE=se,asymp.LCL=est-1.96*se,asymp.UCL=est+1.96*se)}))
write_csv(trends,'tables/fledging_relLD_age_emtrends.csv')
# Fixed-effect predictions at zero random effects. CIs are constructed on the
# logit/log link scale and then back-transformed; they are not predictive bands.
pred<-function(grid,comp){
 if(comp=='location'){X<-model.matrix(~rel_LD*age_class+temp_mean_c+year_c,grid);bb<-fixef(m)$cond;vv<-vcov(m)$cond;inv<-plogis}
 else{X<-model.matrix(~rel_LD+temp_mean_c+year_c+age_class,grid);bb<-fixef(m)$disp;vv<-vcov(m)$disp;inv<-exp}
 X<-X[,names(bb),drop=FALSE];eta<-drop(X%*%bb);se<-sqrt(rowSums((X%*%vv)*X))
 transform(grid,estimate=inv(eta),SE_link=se,CI_low=inv(eta-1.96*se),CI_high=inv(eta+1.96*se))
}
grid<-function(n)data.frame(rel_LD=rep(0,n),age_class=factor(rep('YOUNG',n),levels=c('YOUNG','OLD')),temp_mean_c=rep(0,n),year_c=rep(0,n))
center<-read.csv('tables/fledging_bb_centring_constants.csv')$temp_mean[1]
tt<-grid(160);tt$temp_mean_c<-seq(quantile(d$temp_mean_c,.05),quantile(d$temp_mean_c,.95),length.out=160);tt$temp_mean<-tt$temp_mean_c+center
dt<-pred(tt,'location');dp<-pred(tt,'precision')
write_csv(dt[,c('temp_mean','estimate','SE_link','CI_low','CI_high')],'tables/fig5_1_success_temperature.csv')
write_csv(dp[,c('temp_mean','estimate','SE_link','CI_low','CI_high')],'tables/fig5_1b_phi_temperature.csv')
ld<-grid(320);ld$rel_LD<-rep(seq(quantile(d$rel_LD,.05),quantile(d$rel_LD,.95),length.out=160),2);ld$age_class<-factor(rep(c('YOUNG','OLD'),each=160),levels=c('YOUNG','OLD'))
dl<-pred(ld,'location');ds<-pred(ld,'precision')
write_csv(dl[,c('rel_LD','age_class','estimate','SE_link','CI_low','CI_high')],'tables/fig5_2_success_rellD_age.csv')
write_csv(ds[,c('rel_LD','age_class','estimate','SE_link','CI_low','CI_high')],'tables/fig5_3_phi_rellD.csv')
ag<-grid(2);ag$age_class<-factor(c('YOUNG','OLD'),levels=c('YOUNG','OLD'))
da<-pred(ag,'location');db<-pred(ag,'precision')
write_csv(da[,c('age_class','estimate','SE_link','CI_low','CI_high')],'tables/fig5_3a_success_age.csv')
write_csv(db[,c('age_class','estimate','SE_link','CI_low','CI_high')],'tables/fig5_3b_phi_age.csv')
years<-d[!duplicated(d$year),c('year','year_c')];years<-years[order(years$year),]
yr<-grid(nrow(years));yr$year<-years$year;yr$year_c<-years$year_c
dy<-pred(yr,'location');dz<-pred(yr,'precision')
write_csv(dy[,c('year','estimate','SE_link','CI_low','CI_high')],'tables/fig5_4a_success_year.csv')
write_csv(dz[,c('year','estimate','SE_link','CI_low','CI_high')],'tables/fig5_4b_phi_year.csv')
blue<-'#0072B2';orange<-'#D55E00';th<-theme_classic(base_size=11)
curve<-function(dat,x,ylab,col){ggplot(dat,aes(x=.data[[x]],y=estimate))+geom_ribbon(aes(ymin=CI_low,ymax=CI_high),fill=col,alpha=.2)+geom_line(colour=col,linewidth=1)+labs(x=switch(x,temp_mean='Mean spring temperature (°C)',rel_LD='Relative laying date (days)',year='Year'),y=ylab)+th}
pt<-curve(dt,'temp_mean','Predicted egg-to-fledging success',blue);pp<-curve(dp,'temp_mean','Consistency (φ)',orange)
ggsave('figures/Fig5_1_success_phi_temperature.png',pt+pp,width=10,height=4.6,dpi=400)
ggsave('figures/Fig5_1_success_temperature.png',pt,width=6.4,height=4.8,dpi=400)
ldcurve<-function(dat,ylab){ggplot(dat,aes(rel_LD,estimate,colour=age_class,fill=age_class))+geom_ribbon(aes(ymin=CI_low,ymax=CI_high),alpha=.16,colour=NA)+geom_line(linewidth=1)+scale_colour_manual(values=c(YOUNG=blue,OLD=orange))+scale_fill_manual(values=c(YOUNG=blue,OLD=orange))+labs(x='Relative laying date (days)',y=ylab,colour='Age class',fill='Age class')+th+theme(legend.position='top')}
pl<-ldcurve(dl,'Predicted egg-to-fledging success');ps<-ldcurve(ds,'Consistency (φ)')
ggsave('figures/Fig5_2_success_rellD_age.png',pl+ps,width=10,height=4.6,dpi=400)
ageplot<-function(dat,ylab){ggplot(dat,aes(age_class,estimate,colour=age_class))+geom_errorbar(aes(ymin=CI_low,ymax=CI_high),width=.12)+geom_point(size=3)+scale_colour_manual(values=c(YOUNG=blue,OLD=orange),guide='none')+labs(x='Female age class',y=ylab)+th}
ggsave('figures/Fig5_3b_phi_age.png',ageplot(da,'Predicted egg-to-fledging success')+ageplot(db,'Consistency (φ)'),width=10,height=4.6,dpi=400)
py<-curve(dy,'year','Predicted egg-to-fledging success',blue);pz<-curve(dz,'year','Consistency (φ)',orange)
ggsave('figures/Fig5_4_time_location_scale.png',py+pz,width=10,height=4.6,dpi=400)
# Combined manuscript figures retain primary Bayesian clutch predictions.
cs<-function(file)read.csv(paste0('tables/',file))
quad<-function(a,b,c,e,name){ggsave(paste0('figures/',name,'.png'),((a+b)/(c+e))+plot_annotation(tag_levels='A'),width=10,height=8,dpi=400)}
quad(curve(cs('fig4_1_clutch_temperature.csv'),'temp_mean','Predicted clutch size',blue),curve(cs('fig4_1b_sigma_temperature.csv'),'temp_mean','Inconsistency (σ)',orange),pt,pp,'Fig3_location_scale_predictions')
quad(curve(cs('fig4_4a_annual_clutch.csv'),'year','Predicted clutch size',blue),curve(cs('fig4_4b_annual_sigma.csv'),'year','Inconsistency (σ)',orange),py,pz,'Fig5_reproductive_time_location_scale')
quad(curve(cs('fig4_2a_clutch_rellD.csv'),'rel_LD','Predicted clutch size',blue),curve(cs('fig4_2b_sigma_rellD.csv'),'rel_LD','Inconsistency (σ)',orange),pl,ps,'Fig6_relative_laying_date_location_scale')
cat('Exported verified reference coefficients and updated predictions.\n')
