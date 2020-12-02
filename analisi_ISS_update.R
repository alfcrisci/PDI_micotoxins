#########################################################################################################
library(readxl)
library(fitdistrplus)
library(pastecs)
library(psych)
#########################################################################################################

# setwd("/home/alf/Scrivania/lav_toscano_ISS/full_PDI/codice")


source("aux_analisi_ISS.R")

ex_adults_m_f=readRDS("final_data/ex_adults_m_f.rds")


number_biomarker=lapply(ex_adults_m_f,function(x) length(na.omit((x$UB))))

res_pooled_UB=list()
res_pooled_LB=list()
res_UB=list()
res_LB=list()


res_desc=lapply(ex_adults_m_f,function(x) {c(mean(x$Mean,na.rm=T),mean(x$SD,na.rm=T),
                                             mean(x$UB,na.rm=T),sd(x$UB,na.rm=T),
                                             mean(x$LB,na.rm=T),sd(x$UB,na.rm=T))})

###############################################################################################
# con più di 6 dati

moredata=ex_adults_m_f[c(which(as.numeric(unlist(number_biomarker))>6))]

###############################################################################################
# con meno di 6 dati

poordata=ex_adults_m_f[c(which(as.numeric(unlist(number_biomarker))<=6))]

###############################################################################################

res_poor=lapply(poordata,function(x) {c(mean(x$Mean,na.rm=T),
                                        mean(x$SD,na.rm=T),
                                        mean(x$"resLOD_Method LOD",na.rm=T))
                                        })

res_more=lapply(moredata,function(x) {c(mean(x$Mean,na.rm=T),
                                        mean(x$SD,na.rm=T),
                                        mean(x$"resLOD_Method LOD",na.rm=T))
                                     })
###############################################################################################

res_pooled_UB=lapply(ex_adults_m_f,function(x) {as.numeric(na.omit(x$UB[x$UB>0]))})
res_pooled_LB=lapply(ex_adults_m_f,function(x) {as.numeric(na.omit(ifelse(x$LB==0,0.0001,x$LB)))})
res_pooled_mean_poor=lapply(res_poor,function(x) resample_function_weib(x[1],x[2],N=500))
res_pooled_mean_more=lapply(res_more,function(x) resample_function_weib(x[1],x[2],N=500))
###############################################################################################



num_UB=lapply(res_pooled_UB,length)
num_LB=lapply(res_pooled_LB,length)
num_mean=lapply(res_pooled_mean_poor,length) #500
num_mean_more=lapply(res_pooled_mean_more,length) #500



############################################################################################################
# Qui si trattano i dati upper bound

res_UB=list()
res_dists_UB=list()
res_param_UB=list()
res_summary_UB=list()
res_names=names(res_pooled_UB)


for ( i in 1:length(res_pooled_UB)) {
  
  x=as.numeric(res_pooled_UB[[i]])
  xboot=sample(x,size=500, replace = TRUE)
  
  png(paste0("Upper_bound_",as.character(res_names[i]),".png"))
  res_UB[[i]]=descdist(xboot, boot = 1000)
  title(paste("\n\nUpper bound ",as.character(res_names[i])))
  dev.off()
  
  #######################################################
  
  dists=list()
  idlist=c("weibull","exp","norm")

  fw <- try(fitdist(xboot,"weibull"))
  fexp <- try(fitdist(xboot,"exp"))
  flnor <- try(fitdist(xboot,"norm"))
  
  dists[[1]]=fw;dists[[2]]= fexp;dists[[3]]=flnor
  options(warn=0)
  
  for ( jj in idlist ) {
    outfilepdf=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_UB.pdf"))
    outfilepng=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_UB.png"))
    
    ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"density"))
    qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"QQ-plot"))
    cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"CDF"))
    pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"PP-plot"))
    
    
    out=ggarrange(ds, qq, cd, pp)
    ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
    ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")
  }
  
  res_dists_UB[[i]]=dists
  res_param_UB[[i]]=c(dists[[1]]$estimate,dists[[1]]$sd,dists[[2]]$estimate,dists[[2]]$sd,dists[[3]]$estimate,dists[[3]]$sd);
  res_summary_UB[[i]]=summary_large(x)
}

##########################################################################################################################################
# Qui si trattano i dati lower bound

res_LB=list()
res_dists_LB=list()
res_param_LB=list()
res_summary_LB=list()
res_names=names(res_pooled_LB)

for ( i in 1:length(res_pooled_LB)) {
  
  x=as.numeric(res_pooled_LB[[i]])
  xboot=sample(x,size=500, replace = TRUE)
  
  png(paste0("Lower_bound_",as.character(res_names[i]),".png"))
  res_LB[[i]]=descdist(xboot, boot = 1000)
  title(paste("\n\nLower bound ",as.character(res_names[i])))
  dev.off()
  
  #######################################################
  
  dists=list()
  idlist=c("weibull","exp","norm")
  fw <- try(fitdist(xboot,"weibull"))
  fexp <- try(fitdist(xboot,"exp"))
  flnor <- try(fitdist(xboot,"norm"))
  
  dists[[1]]=fw;dists[[2]]= fexp;dists[[3]]=flnor
  
  options(warn=0)
  
  for ( jj in idlist ) {
    outfilepdf=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_LB.pdf"))
    outfilepng=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_LB.png"))
    
    ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"density"))
    qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"QQ-plot"))
    cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"CDF"))
    pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"PP-plot"))
    
    
    out=ggarrange(ds, qq, cd, pp)
    ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
    ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")
  }
  
  res_dists_LB[[i]]=dists
  res_param_LB[[i]]=c(dists[[1]]$estimate,dists[[1]]$sd,dists[[2]]$estimate,dists[[2]]$sd,dists[[3]]$estimate,dists[[3]]$sd);
  res_summary_LB[[i]]=summary_large(x)
  
}

##########################################################################################################################################
# lavoro su popolaizoni con meno di 7 dati

res_mean=list()
res_dists_mean=list()
res_param_mean=list()
res_summary_mean=list()
res_names=names(res_pooled_mean_poor)

for ( i in 1:length(res_pooled_mean_poor)) {
  
  x=as.numeric(res_pooled_mean_poor[[i]])
  xboot=sample(x,size=500, replace = TRUE)
  
  #######################################################
  png(paste0("Mean_under7_",as.character(res_names[i]),".png"))
  res_mean[[i]]=descdist(xboot, boot = 1000)
  title(paste("\n\nMean fit ",as.character(res_names[i])))
  dev.off()
  
  dists=list()
  idlist=c("weibull","exp","norm")
  fw <- try(fitdist(xboot,"weibull"))
  fexp <- try(fitdist(xboot,"exp"))
  flnor <- try(fitdist(xboot,"norm"))
  dists[[1]]=fw;dists[[2]]= fexp;dists[[3]]=flnor
  options(warn=0)
  
  for ( jj in idlist ) {
    outfilepdf=paste0("Plot_under7","_",paste0(as.character(res_names[i]),"_dists_mean.pdf"))
    outfilepng=paste0("Plot_under7","_",paste0(as.character(res_names[i]),"_dists_mean.png"))
    
    ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"density"))
    qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"QQ-plot"))
    cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"CDF"))
    pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"PP-plot"))
    
    out=ggarrange(ds, qq, cd, pp)
    ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
    ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")
  }
  
  res_dists_mean[[i]]=dists
  res_param_mean[[i]]=c(dists[[1]]$estimate,dists[[1]]$sd,dists[[2]]$estimate,dists[[2]]$sd,dists[[3]]$estimate,dists[[3]]$sd);
  res_summary_mean[[i]]=summary_large(x)
  
}

##########################################################################################################################################
# Qui si trattano con quelli con più di 6 dati

res_mean_more=list()
res_dists_mean_more=list()
res_param_mean_more=list()
res_summary_mean_more=list()
res_names=names(res_pooled_mean_more)

for ( i in 1:length(res_pooled_mean_more)) {
  
  x=as.numeric(res_pooled_mean_more[[i]])
  xboot=sample(x,size=500, replace = TRUE)
  
  #######################################################
  png(paste0("Mean_over7_",as.character(res_names[i]),".png"))
  res_mean_more[[i]]=descdist(x, boot = 1000)
  title(paste("\n\nMean fit ",as.character(res_names[i])))
  dev.off()
  
  dists=list()
  idlist=c("weibull","exp","norm")
  fw <- try(fitdist(xboot,"weibull"))
  fexp <- try(fitdist(xboot,"exp"))
  flnor <- try(fitdist(xboot,"norm"))
  dists[[1]]=fw;dists[[2]]= fexp;dists[[3]]=flnor
  options(warn=0)
  
  for ( jj in idlist ) {
    outfilepdf=paste0("Plot_over7","_",paste0(as.character(res_names[i]),"_dists_mean.pdf"))
    outfilepng=paste0("Plot_over7","_",paste0(as.character(res_names[i]),"_dists_mean.png"))
    
    ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"density"))
    qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"QQ-plot"))
    cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"CDF"))
    pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"PP-plot"))
    
    out=ggarrange(ds, qq, cd, pp)
    ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
    ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")
  }
  
  res_dists_mean_more[[i]]=dists
  res_param_mean_more[[i]]=c(dists[[1]]$estimate,dists[[1]]$sd,dists[[2]]$estimate,dists[[2]]$sd,dists[[3]]$estimate,dists[[3]]$sd);
  res_summary_mean_more[[i]]=summary_large(x)
}

##########################################################################################################################################

names_dist=c("myco","weibull_shape_mean","weibull_scale_mean","weibull_shape_stderr","weibull_scale_stderr","exp_rate_mean","exp_rate_stderr","norm_mean_mean","norm_mean_sd","norm_mean__stderr","norm_mean__stderr")

res_param_UB_df=data.frame(myco=names(res_pooled_UB),do.call("rbind",res_param_UB))
res_param_LB_df=data.frame(myco=names(res_pooled_LB),do.call("rbind",res_param_LB))
res_param_mean_poor_df=data.frame(myco=names(res_pooled_mean_poor),do.call("rbind",res_param_mean))
res_param_mean_more_df=data.frame(myco=names(res_pooled_mean_more),do.call("rbind",res_param_mean_more))

names(res_param_UB_df)=names_dist
names(res_param_LB_df)=names_dist
names(res_param_mean_poor_df)=names_dist
names(res_param_mean_more_df)=names_dist

res_summary_UB_df=data.frame(myco_class=names(res_pooled_UB),do.call("rbind",res_summary_UB))
res_summary_LB_df=data.frame(myco_class=names(res_pooled_LB),do.call("rbind",res_summary_LB))
res_summary_mean_poor_df=data.frame(myco_class=names(res_pooled_mean_poor),do.call("rbind",res_summary_mean))
res_summary_mean_more_df=data.frame(myco_class=names(res_pooled_mean_more),do.call("rbind",res_summary_mean_more))

##################################################################################################

file.remove(c("fitparams_UB.xls","fitparams_LB.xls","fitparams_mean.xls"))

XLConnect::writeWorksheetToFile("fitparams_UB.xls",res_summary_UB_df,"res_summary_UB")
XLConnect::writeWorksheetToFile("fitparams_UB.xls",res_param_UB_df,"res_param_UB")

XLConnect::writeWorksheetToFile("fitparams_LB.xls",res_summary_LB_df,"res_summary_LB")
XLConnect::writeWorksheetToFile("fitparams_LB.xls",res_param_LB_df,"res_param_LB")

XLConnect::writeWorksheetToFile("fitparams_mean.xls",res_summary_mean_poor_df,"summary_less_7_data")
XLConnect::writeWorksheetToFile("fitparams_mean.xls",res_param_mean_poor_df,"fit_less_7_data")

XLConnect::writeWorksheetToFile("fitparams_mean.xls",res_summary_mean_more_df,"summary_more_7_data")
XLConnect::writeWorksheetToFile("fitparams_mean.xls",res_param_mean_more_df,"fit_more_7_data")


#############################################################################################################################
# check fit model of distribution

#########################################################################
# [1] "Reference"           "Mycotoxin_INFO"      "paramType_Biomarker"
# [4] "Country"             "SampleSize"          "SampleType_INFO"    
# [7] "Age"                 "Mean"                "Median"             
# [10] "GeometricMean"       "SD"                  "Percentile"         
# [13] "Value_INFO"          "Min"                 "Max"                
# [16] "resLOQ_Method LOQ"   "resLOD_Method LOD"   "LB"                 
# [19] "UB"                 
#########################################################################


#########################################################################################
# references

# [1] https://math.stackexchange.com/questions/2487059/whats-the-kurtosis-of-exponential-distribution
# https://rstudio-pubs-static.s3.amazonaws.com/208109_47ee36af179343658be360532b3afde3.html
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
# https://stats.stackexchange.com/questions/100862/weibull-distribution-from-given-mean
# https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# https://stats.stackexchange.com/questions/159452/how-can-i-recreate-a-weibull-distribution-given-mean-and-standard-deviation-and
# https://stats.stackexchange.com/questions/60511/weibull-distribution-parameters-k-and-c-for-wind-speed-data
