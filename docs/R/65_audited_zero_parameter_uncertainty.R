Sys.setlocale('LC_CTYPE','English_United States.utf8')
.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages(library(glmmTMB))
out<-'manuscript_audit_2026-10-13/zero_frequency_audit'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
mods<-list(reference_year_precision=readRDS('models/glmmTMB/fledging_bb_reference_year_precision.rds'),zero_augmented=readRDS('manuscript_audit_2026-10-08/zero_frequency_audit/year_precision_zero_inflated_intercept.rds'))
rows<-list()
set.seed(20261014)
for(label in names(mods)){
 m<-mods[[label]];d<-m$frame;n<-rowSums(d[[1]]);obs<-mean(d[[1]][,1]==0)
 par<-m$fit$par;V<-m$sdr$cov.fixed;stopifnot(length(par)==nrow(V),isTRUE(m$sdr$pdHess));L<-chol(V)
 pp<-sweep(matrix(rnorm(1000*length(par)),1000,length(par))%*%L,2,par,'+')
 X<-model.matrix(~rel_LD*age_class+temp_mean_c+year_c,d);X<-X[,names(fixef(m)$cond),drop=FALSE]
 Z<-model.matrix(~rel_LD+temp_mean_c+year_c+age_class,d);Z<-Z[,names(fixef(m)$disp),drop=FALSE]
 iy<-as.integer(d$year_f);iff<-as.integer(d$female_id)
 # The theta order follows the conditional random-effect groups used by glmmTMB.
 cn<-names(VarCorr(m)$cond);stopifnot(identical(cn,c('year_f','female_id')))
 z<-vapply(1:1000,function(i){p<-pp[i,];eta<-drop(X%*%p[names(par)=='beta']);ep<-drop(Z%*%p[names(par)=='betadisp']);th<-p[names(par)=='theta'];u<-rnorm(nlevels(d$year_f),0,exp(th[1]));v<-rnorm(nlevels(d$female_id),0,exp(th[2]));td<-p[names(par)=='thetadisp'];a<-if(length(td))rnorm(nlevels(d$year_f),0,exp(td[1]))else rep(0,nlevels(d$year_f));mu<-plogis(eta+u[iy]+v[iff]);phi<-exp(ep+a[iy]);prob<-rbeta(length(n),mu*phi,(1-mu)*phi);y<-rbinom(length(n),n,prob);zi<-p[names(par)=='betazi'];if(length(zi))y[rbinom(length(n),1,plogis(zi[1]))==1]<-0;mean(y==0)},numeric(1))
 rows[[label]]<-data.frame(model=label,observed=obs,median=median(z),lower=quantile(z,.025),upper=quantile(z,.975),n_ge=sum(z>=obs),nsim=1000)
 write.csv(data.frame(rep=1:1000,zero=z),file.path(out,paste0(label,'_parameter_uncertainty_simulations.csv')),row.names=FALSE)
}
tab<-do.call(rbind,rows);print(tab);write.csv(tab,file.path(out,'parameter_uncertainty_summary.csv'),row.names=FALSE)
