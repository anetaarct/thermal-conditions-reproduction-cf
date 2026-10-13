#!/usr/bin/env Rscript
Sys.setlocale('LC_CTYPE','English_United States.utf8')
.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages(library(ggplot2))
wins<-c('main','april','short','long','selected')
cs<-do.call(rbind,lapply(wins,function(w)read.csv(paste0('tables/window_sensitivity_clutch_',w,'.csv'))))
cs<-cs[cs$parameter%in%c('b[2]','b_sigma[2]'),]
cs$component<-ifelse(cs$parameter=='b[2]','Clutch mean (eggs / C)','Clutch log residual SD (/ C)')
cs$interval_type<-'95% posterior credible interval'
fs<-read.csv('tables/window_sensitivity_fledging.csv')
fs$component<-ifelse(fs$component=='cond','Fledging mean (logit / C)','Fledging precision (log phi / C)')
fs$interval_type<-'95% Wald confidence interval'
cols<-c('window','component','estimate','SE','lower95','upper95','interval_type')
all<-rbind(cs[,cols],fs[,cols]);write.csv(all,'tables/window_sensitivity_comparison.csv',row.names=FALSE)
labels<-c(main='26 Mar - 7 May (primary)',april='1 - 30 Apr',short='9 Apr - 7 May',long='5 Mar - 7 May',selected='Full selected weekly bins')
all$window_label<-factor(labels[all$window],levels=rev(labels))
p<-ggplot(all,aes(estimate,window_label))+geom_vline(xintercept=0,linetype=2,colour='grey55')+
 geom_segment(aes(x=lower95,xend=upper95,yend=window_label),linewidth=.7)+
 geom_point(aes(colour=window=='main'),size=2.6)+scale_colour_manual(values=c('FALSE'='#526678','TRUE'='#007C83'),guide='none')+
 facet_wrap(~component,scales='free_x',ncol=2)+labs(x='Temperature coefficient per 1 °C (95% interval)',y=NULL,
 title='Sensitivity to the temperature window',subtitle='Fledging reference: year random effects in both location and precision',
 caption='Clutch: posterior CrIs; fledging: Wald CIs. Intervals describe coefficient uncertainty, not raw outcome variation.')+
 theme_bw(base_size=11)+theme(panel.grid.minor=element_blank())
ggsave('figures/window_sensitivity_comparison.png',p,width=11,height=7,dpi=200)
cat('Updated all five window comparisons and plotting data.\n')
