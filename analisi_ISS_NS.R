#########################################################################################################
library(readxl)
library(fitdistrplus)
library(pastecs)
library(psych)
library(KScorrect)
library(ggplot2)
library(ggpubr)
#########################################################################################################


setwd("/home/alf/Scrivania/codice_dati_PDI/PDI_micotoxins")
source("aux_analisi_ISS.R")
set.seed(2)

ex_adults_ad=rec_excel("final_data/extract_adults_AGGREGATED_NS_CHK.xlsx")
saveRDS(ex_adults_ad,file="final_data/ex_adults_ad_NS.rds")


number_biomarker=lapply(ex_adults_ad,function(x) length(na.omit((x$UB))))

res_pooled_UB=list()
res_pooled_LB=list()
res_UB=list()
res_LB=list()

res_desc=lapply(ex_adults_ad,function(x) {c(mean(x$Mean,na.rm=T),mean(x$SD,na.rm=T),
                                             mean(x$UB,na.rm=T),sd(x$UB,na.rm=T),
                                             mean(x$LB,na.rm=T),sd(x$LB,na.rm=T))})

###############################################################################################
# con più di 6 dati

moredata=ex_adults_ad[c(which(as.numeric(unlist(number_biomarker))>6))]


###############################################################################################



res_more=lapply(moredata,function(x) {c(mean(x$Mean,na.rm=T),
                                        mean(x$SD,na.rm=T),
                                        mean(x$"resLOD_Method LOD",na.rm=T))
                                     })
###############################################################################################

res_pooled_UB=lapply(ex_adults_ad,function(x) {as.numeric(na.omit(x$UB[x$UB>0]))})
res_pooled_LB=lapply(ex_adults_ad,function(x) {as.numeric(na.omit(ifelse(x$LB==0,0.0001,x$LB)))})
res_pooled_mean_more=lapply(res_more,function(x) resample_function_weib(x[1],x[2],N=500))


###############################################################################################


num_UB=lapply(res_pooled_UB,length)
num_LB=lapply(res_pooled_LB,length)
num_mean_more=lapply(res_pooled_mean_more,length) #500

############################################################################################################
dir.create("plots_NS")
setwd("plots_NS")
############################################################################################################
# UB upper bound data

res_UB=list()
res_dists_UB=list()
res_param_UB=list()
res_summary_UB=list()
res_names=names(res_pooled_UB)
res_gof_dists_UB=list()
res_kgof_dists_UB=list()

for ( i in 1:length(res_pooled_UB)) {
  
  x=as.numeric(res_pooled_UB[[i]])
  xboot=sample(x,size=500, replace = TRUE)
  
  
  #######################################################
  png(paste0("Upper_bound_",as.character(res_names[i]),"_NS.png"))
  res_UB[[i]]=descdist(xboot, boot = 1000)
  title(paste("\n\nUpper bound ",as.character(res_names[i])))
  dev.off()
  #######################################################
  
   
  dists=list()
  kgof_dists=list()
  gof_dists=list()
  idlist=c("weibull","exp","norm")

  fw <- try(fitdist(xboot,"weibull"))
  fexp <- try(fitdist(xboot,"exp"))
  flnor <- try(fitdist(xboot,"norm"))
  gof_fw <- try(gofstat(fw))
  gof_fexp <- try(gofstat(fexp))
  gof_flnor <- try(gofstat(flnor))
  kgof_fw <-LcKS(x, "pweibull",nreps=1000,parallel = TRUE)
  kgof_fexp <- LcKS(x, "pexp",nreps=1000,parallel = TRUE)
  kgof_flnor <- LcKS(x, "pnorm",nreps=1000,parallel = TRUE)
  
  dists[[1]]=fw;dists[[2]]= fexp;dists[[3]]=flnor
  kgof_dists[[1]]=kgof_fw;kgof_dists[[2]]= kgof_fexp;kgof_dists[[3]]=kgof_flnor
  gof_dists[[1]]=gof_fw;gof_dists[[2]]= gof_fexp;gof_dists[[3]]=gof_flnor
  options(warn=0)
  
  for ( jj in idlist ) {
    outfilepdf=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_UB_NS.pdf"))
    outfilepng=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_UB_NS.png"))
    
    ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"density"))
    qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"QQ-plot"))
    cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"CDF"))
    pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"PP-plot"))
    
    
   out=ggarrange(ds, qq, cd, pp)
   ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
   ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")
  }
  
  res_dists_UB[[i]]=dists
  res_param_UB[[i]]=c(dists[[1]]$estimate[1],dists[[1]]$estimate[2],
                      dists[[1]]$sd[1],dists[[1]]$sd[2],
                      dists[[2]]$estimate,dists[[2]]$sd,
                      dists[[3]]$estimate,dists[[3]]$sd);
  
  
  res_summary_UB[[i]]=stat.desc(as.numeric(res_pooled_UB[[i]]))

  res_gof_dists_UB[[i]]=c(gof_fw$ks,gof_fw$kstest,gof_fw$aic,
                          gof_fexp$ks,gof_fexp$kstest,gof_fexp$aic,
                          gof_flnor$ks,gof_flnor$kstest,gof_flnor$aic);
  res_kgof_dists_UB[[i]]=c(kgof_fw$D.obs,kgof_fw$p.value,
                           kgof_fexp$D.obs,kgof_fexp$p.value,
                           kgof_flnor$D.obs,kgof_flnor$p.value);
}




##########################################################################################################################################
# LB lower bound data

res_LB=list()
res_dists_LB=list()
res_param_LB=list()
res_summary_LB=list()
res_names=names(res_pooled_LB)
res_gof_dists_LB=list()
res_kgof_dists_LB=list()

for ( i in 1:length(res_pooled_LB)) {
  
  x=as.numeric(res_pooled_LB[[i]])
  outliers <- boxplot(x)$out
  if ((length(outliers) >4) & (length(x) >10)) {x=x[-which(x %in% outliers)]}
  xboot=sample(x,size=500, replace = TRUE)
  
  #######################################################
  png(paste0("Lower_bound_",as.character(res_names[i]),"_NS.png"))
  res_LB[[i]]=descdist(xboot, boot = 1000)
  title(paste("\n\nLower bound ",as.character(res_names[i])))
  dev.off()
  #######################################################
  
  dists=list()
  kgof_dists=list()
  gof_dists=list()
  idlist=c("weibull","exp","norm")
  
  fw <- try(fitdist(xboot,"weibull"))
  fexp <- try(fitdist(xboot,"exp"))
  flnor <- try(fitdist(xboot,"norm"))
  gof_fw <- try(gofstat(fw))
  gof_fexp <- try(gofstat(fexp))
  gof_flnor <- try(gofstat(flnor))
  kgof_fw <-LcKS(x, "pweibull",nreps=1000,parallel = TRUE)
  kgof_fexp <- LcKS(x, "pexp",nreps=1000,parallel = TRUE)
  kgof_flnor <- LcKS(x, "pnorm",nreps=1000,parallel = TRUE)
  
  dists[[1]]=fw;dists[[2]]= fexp;dists[[3]]=flnor
  kgof_dists[[1]]=kgof_fw;kgof_dists[[2]]= kgof_fexp;kgof_dists[[3]]=kgof_flnor
  gof_dists[[1]]=gof_fw;gof_dists[[2]]= gof_fexp;gof_dists[[3]]=gof_flnor
  options(warn=0)
  
  
  for ( jj in idlist ) {
    outfilepdf=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_LB_NS.pdf"))
    outfilepng=paste0("Plot","_",paste0(as.character(res_names[i]),"_dists_LB_NS.png"))
    
    ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"density"))
    qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"QQ-plot"))
    cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"CDF"))
    pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"PP-plot"))
    
    
    out=ggarrange(ds, qq, cd, pp)
    ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
    ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")
  }
  
  res_dists_LB[[i]]=dists
  res_param_LB[[i]]=c(dists[[1]]$estimate[1],dists[[1]]$estimate[2],
                      dists[[1]]$sd[1],dists[[1]]$sd[2],
                      dists[[2]]$estimate,dists[[2]]$sd,
                      dists[[3]]$estimate,dists[[3]]$sd);
  res_summary_LB[[i]]=stat.desc(as.numeric(res_pooled_LB[[i]]))
  
  res_gof_dists_LB[[i]]=c(gof_fw$ks,gof_fw$kstest,gof_fw$aic,
                          gof_fexp$ks,gof_fexp$kstest,gof_fexp$aic,
                          gof_flnor$ks,gof_flnor$kstest,gof_flnor$aic);
  res_kgof_dists_LB[[i]]=c(kgof_fw$D.obs,kgof_fw$p.value,
                      kgof_fexp$D.obs,kgof_fexp$p.value,
                      kgof_flnor$D.obs,kgof_flnor$p.value);
  
  
}

##########################################################################################################################################
##########################################################################################################################################
# mean values more than 6  data

res_mean_more=list()
res_dists_mean_more=list()
res_gof_dists_mean_more=list()
res_kgof_dists_mean_more=list()
res_param_mean_more=list()
res_summary_mean_more=list()
res_names=names(res_pooled_mean_more)


for ( i in 1:length(res_pooled_mean_more)) {
  
  x=as.numeric(res_pooled_mean_more[[i]])
  xboot=sample(x,size=500, replace = TRUE)
  
  #######################################################
  png(paste0("Mean_over_",as.character(res_names[i]),"_NS.png"))
  res_mean_more[[i]]=descdist(x, boot = 1000)
  title(paste("\n\nMean fit ",as.character(res_names[i])))
  dev.off()
  #######################################################
  
  dists=list()
  kgof_dists=list()
  gof_dists=list()
  idlist=c("weibull","exp","norm")
  
  fw <- try(fitdist(xboot,"weibull"))
  fexp <- try(fitdist(xboot,"exp"))
  flnor <- try(fitdist(xboot,"norm"))
  gof_fw <- try(gofstat(fw))
  gof_fexp <- try(gofstat(fexp))
  gof_flnor <- try(gofstat(flnor))
  kgof_fw <-LcKS(x, "pweibull",nreps=1000,parallel = TRUE)
  kgof_fexp <- LcKS(x, "pexp",nreps=1000,parallel = TRUE)
  kgof_flnor <- LcKS(x, "pnorm",nreps=1000,parallel = TRUE)
  
  dists[[1]]=fw;dists[[2]]= fexp;dists[[3]]=flnor
  kgof_dists[[1]]=kgof_fw;kgof_dists[[2]]= kgof_fexp;kgof_dists[[3]]=kgof_flnor
  gof_dists[[1]]=gof_fw;gof_dists[[2]]= gof_fexp;gof_dists[[3]]=gof_flnor
  options(warn=0)
  
  
  for ( jj in idlist ) {
    outfilepdf=paste0("Plot_over","_",paste0(as.character(res_names[i]),"_dists_mean_NS.pdf"))
    outfilepng=paste0("Plot_over","_",paste0(as.character(res_names[i]),"_dists_mean_NS.png"))
    
    ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"density"))
    qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"QQ-plot"))
    cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"CDF"))
    pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot",main=paste(as.character(res_names[i]),"PP-plot"))
    
    out=ggarrange(ds, qq, cd, pp)
    ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
    ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")
  }
  
  res_dists_mean_more[[i]]=dists
  res_param_mean_more[[i]]=c(dists[[1]]$estimate[1],dists[[1]]$estimate[2],
                             dists[[1]]$sd[1],dists[[1]]$sd[2],
                             dists[[2]]$estimate,dists[[2]]$sd,
                             dists[[3]]$estimate,dists[[3]]$sd);
  res_summary_mean_more[[i]]=stat.desc(as.numeric(res_pooled_mean_more[[i]]))
  res_gof_dists_mean_more[[i]]=c(gof_fw$ks,gof_fw$kstest,gof_fw$aic,
                                gof_fexp$ks,gof_fexp$kstest,gof_fexp$aic,
                                gof_flnor$ks,gof_flnor$kstest,gof_flnor$aic);
  res_kgof_dists_mean_more[[i]]=c(kgof_fw$D.obs,kgof_fw$p.value,
                             kgof_fexp$D.obs,kgof_fexp$p.value,
                             kgof_flnor$D.obs,kgof_flnor$p.value);
}

##########################################################################################################################################
setwd("..")
##########################################################################################################################################
# Organize results

names_dist=c("myco",
             "weibull_shape_mean",
             "weibull_scale_mean",
             "weibull_shape_stderr",
             "weibull_scale_stderr",
             "exp_rate_mean",
             "exp_rate_stderr",
             "norm_mean_mean",
             "norm_mean_sd",
             "norm_mean_stderr",
             "norm_sd_stderr")



names_gof=c("myco",
            "weibull_D",
            "weibull_kstest",
            "weibull_aic",
            "exp_D",
            "exp_kstest",
            "exp_AIC",
            "norm_D",
            "norm_kstest",
            "norm_AIC")

names_kgof=c("myco",
             "weibull_D",
             "weibull_pvalue",
             "exp_D",
             "exp_pvalue",
             "norm_D",
             "norm_pvalue")

res_param_UB_df=data.frame(myco=names(res_pooled_UB),do.call("rbind",res_param_UB))
res_param_LB_df=data.frame(myco=names(res_pooled_LB),do.call("rbind",res_param_LB))
res_param_mean_more_df=data.frame(myco=names(res_pooled_mean_more),do.call("rbind",res_param_mean_more))

res_gof_dists_UB_df=data.frame(myco=names(res_pooled_UB),do.call("rbind",res_gof_dists_UB))
res_gof_dists_LB_df=data.frame(myco=names(res_pooled_LB),do.call("rbind",res_gof_dists_LB))
res_gof_dists_mean_poor_df=data.frame(myco=names(res_pooled_mean_poor),do.call("rbind",res_gof_dists_mean))
res_gof_dists_mean_more_df=data.frame(myco=names(res_pooled_mean_more),do.call("rbind",res_gof_dists_mean_more))

res_kgof_dists_UB_df=data.frame(myco=names(res_pooled_UB),do.call("rbind",res_kgof_dists_UB))
res_kgof_dists_LB_df=data.frame(myco=names(res_pooled_LB),do.call("rbind",res_kgof_dists_LB))
res_kgof_dists_mean_more_df=data.frame(myco=names(res_pooled_mean_more),do.call("rbind",res_kgof_dists_mean_more))


names(res_param_UB_df)=names_dist
names(res_param_LB_df)=names_dist
names(res_param_mean_more_df)=names_dist


names(res_gof_dists_UB_df)=names_gof
names(res_gof_dists_LB_df)=names_gof
names(res_gof_dists_mean_more_df)=names_gof

names(res_kgof_dists_UB_df)=names_kgof
names(res_kgof_dists_LB_df)=names_kgof
names(res_kgof_dists_mean_more_df)=names_kgof



res_summary_UB_df=data.frame(myco_class=names(res_pooled_UB),do.call("rbind",res_summary_UB))
res_summary_LB_df=data.frame(myco_class=names(res_pooled_LB),do.call("rbind",res_summary_LB))
res_summary_mean_more_df=data.frame(myco_class=names(res_pooled_mean_more),do.call("rbind",res_summary_mean_more))



##################################################################################################

dir.create("gofstats_NS")
setwd("gofstats_NS")

saveRDS(res_dists_UB,"res_dists_UB_NS.rds")
saveRDS(res_dists_LB,"res_dists_LB_NS.rds")
saveRDS(res_dists_mean_more,"res_dists_mean_more_NS.rds")

saveRDS(res_gof_dists_UB,"res_gof_dists_UB_NS.rds")
saveRDS(res_gof_dists_LB,"res_gof_dists_LB_NS.rds")
saveRDS(res_gof_dists_mean_more,"res_gof_dists_mean_more_NS.rds")

saveRDS(res_kgof_dists_UB,"res_kgof_dists_UB_NS.rds")
saveRDS(res_kgof_dists_LB,"res_kgof_dists_LB_NS.rds")
saveRDS(res_kgof_dists_mean_more,"res_kgof_dists_mean_more_NS.rds")

setwd("..")

##################################################################################################

file.remove(c("fitparams_UB_NS.xls","fitparams_LB_NS.xls","fitparams_mean_NS.xls"))
 
XLConnect::writeWorksheetToFile("fitparams_UB_NS.xls",res_summary_UB_df,"res_summary_UB_NS")
XLConnect::writeWorksheetToFile("fitparams_UB_NS.xls",res_param_UB_df,"res_param_UB_NS")
XLConnect::writeWorksheetToFile("fitparams_UB_NS.xls",res_kgof_dists_UB_df,"kgof_dists_UB_NS")
XLConnect::writeWorksheetToFile("fitparams_UB_NS.xls",res_gof_dists_UB_df,"gof_dists_UB_NS")

XLConnect::writeWorksheetToFile("fitparams_LB_NS.xls",res_summary_LB_df,"res_summary_LB_NS")
XLConnect::writeWorksheetToFile("fitparams_LB_NS.xls",res_param_LB_df,"res_param_LB_NS")
XLConnect::writeWorksheetToFile("fitparams_LB_NS.xls",res_kgof_dists_LB_df,"kgof_dists_LB_NS")
XLConnect::writeWorksheetToFile("fitparams_LB_NS.xls",res_gof_dists_LB_df,"gof_dists_LB_NS")



XLConnect::writeWorksheetToFile("fitparams_mean_NS.xls",res_summary_mean_more_df,"res_summary_mean_NS")
XLConnect::writeWorksheetToFile("fitparams_mean_NS.xls",res_param_mean_more_df,"res_param_mean_NS")
XLConnect::writeWorksheetToFile("fitparams_mean_NS.xls",res_kgof_dists_mean_more_df,"kgof_dists_mean_NS")
XLConnect::writeWorksheetToFile("fitparams_mean_NS.xls",res_gof_dists_mean_more_df,"gof_dists_mean_NS")

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
# library(magick)
# setwd("plots/all_table")
# imgs=list.files()
# img_ls=sapply(imgs, image_read)
# img1 <- c(img_ls[[1]], img_ls[[2]], img_ls[[3]], img_ls[[4]])
# img1a <- image_append(img1)
# img2 <- c(img_ls[[5]], img_ls[[6]], img_ls[[7]], img_ls[[8]])
# img2a <- image_append(img2)
# fin=image_append(c(img1a,img2a),stack=T)
# image_write(fin, path = "full.png", format = "png")
#########################################################################################
# code references

# [1] https://math.stackexchange.com/questions/2487059/whats-the-kurtosis-of-exponential-distribution
# [2] https://rstudio-pubs-static.s3.amazonaws.com/208109_47ee36af179343658be360532b3afde3.html
# [3] https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
# [4] ttps://stats.stackexchange.com/questions/100862/weibull-distribution-from-given-mean
# [5] https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# [6] https://stats.stackexchange.com/questions/159452/how-can-i-recreate-a-weibull-distribution-given-mean-and-standard-deviation-and
# [7]  https://stats.stackexchange.com/questions/60511/weibull-distribution-parameters-k-and-c-for-wind-speed-data
