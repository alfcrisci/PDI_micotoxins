library(readxl)
library(doBy)

setwd("/home/alf/Scrivania/codice_dati_PDI/PDI_micotoxins")

###########################################################################################################
# functions

pdi_func=function(x,escr=70,vol_urine=2,weight=70) {return(x*(vol_urine/weight)*(100/escr))}
sumfit=function(x) {return(data.frame(par=x$estimate[1],errpar=x$estimate[2],aic=x$aic,names=x$distname))}
set.seed(2)
###########################################################################################################
# data labeling

EOTA=0.5       # (Schlatter et al., 1996)
EDON=0.72      # (EFSA, 2017, Knutsen HK)
EFB1=0.003     # (Van der Westhuizen et al., 1999)
EZEN=0.368     # (Solfrizzo)
ET2HT2=0.6     # (Se si presenta l’HT2 oltre all’T2 consideriamo la somma) (EFSA 2016).
ENIV=0.01      # (EFSA 2017).
ECIT=0.75      # (ingorare OHCIT) (EFSA 2012).
EAFM1=0.16225  # (ingorare OHCIT) (EFSA 2012).
EAFM1M=0.01705 # 1.23-2.18% MALE Excretion Rate, (ZHU, 1987)
EAFM1F=0.0154  # 1.30-1.78% FEMALE Excretion Rate, (ZHU, 1987) 


###########################################################################################################
# Scenarios weight creation


# https://efsa.onlinelibrary.wiley.com/doi/pdf/10.2903/j.efsa.2012.2579

sL=function(x) {as.numeric(quantile(na.omit(x),c(0.05)))}
sM=function(x) {as.numeric(quantile(na.omit(x),c(0.5)))}
sU=function(x) {as.numeric(quantile(na.omit(x),c(0.95)))}


Wdf=read.csv("final_data/Scenarios_full.csv")

scen_F=Wdf[Wdf$cat=="Female",]
scen_M=Wdf[Wdf$cat=="Male",]
scen_A=Wdf[Wdf$cat=="Adults",]

##################################################################################################
data_PDI=as.data.frame(read_xls("final_data/PDI_params_table_EU.xls",1))
data_PDI$exrate=data_PDI$exrate*100



########################################################################################################################à

calculate_scenario=function(data_PDI,Wdf) {
  
res_pdi_exp_F=list()
res_pdi_weibull_F=list()
res_pdi_norm_F=list()


for ( j in 1:nrow(Wdf))  {
  
res_pdi_exp_w=list()
res_pdi_weibull_w=list()
res_pdi_norm_w=list()

for ( i in 1:nrow(data_PDI)) {

temp_exp_PDI=sapply(rexp(10000,data_PDI$exp_rate_mean[i]),
                    FUN=function(x){ pdi_func(x,escr = data_PDI$exrate[i],weight=Wdf$Wkg[j])})

temp_weibull_PDI=sapply(rweibull(10000, data_PDI$weibull_shape_mean[i], data_PDI$weibull_scale_mean[i]),
                        FUN=function(x){ pdi_func(x,escr = data_PDI$exrate[i],weight=Wdf$Wkg[j])})

temp_norm_PDI=sapply(rnorm(10000, data_PDI$norm_mean_mean[i], data_PDI$norm_mean_sd[i]),
                        FUN=function(x){ pdi_func(x,escr = data_PDI$exrate[i],weight=Wdf$Wkg[j])})


res_pdi_exp_w[[i]]=c(Wdf$scen[j],Wdf$cat[j],Wdf$mycos[j],data_PDI$category[i],data_PDI$myco[i],Wdf$Wkg[j],t.test(temp_exp_PDI)$estimate,t.test(temp_exp_PDI)$conf.int)
res_pdi_weibull_w[[i]]=c(Wdf$scen[j],Wdf$cat[j],Wdf$mycos[j],data_PDI$category[i],data_PDI$myco[i],Wdf$Wkg[j],t.test(temp_weibull_PDI)$estimate,t.test(temp_weibull_PDI)$conf.int)
res_pdi_norm_w[[i]]=c(Wdf$scen[j],Wdf$cat[j],Wdf$mycos[j],data_PDI$category[i],data_PDI$myco[i],Wdf$Wkg[j],t.test(temp_norm_PDI)$estimate,t.test(temp_norm_PDI)$conf.int)

}


res_pdi_exp_F[[j]]=do.call("rbind",res_pdi_exp_w)
res_pdi_weibull_F[[j]]=do.call("rbind",res_pdi_weibull_w)
res_pdi_norm_F[[j]]=do.call("rbind",res_pdi_norm_w)

}

res_pdi_exp_F_df=as.data.frame(do.call("rbind",res_pdi_exp_F))
res_pdi_weibull_F_df=as.data.frame(do.call("rbind",res_pdi_weibull_F))
res_pdi_norm_F_df=as.data.frame(do.call("rbind",res_pdi_norm_F))

names(res_pdi_exp_F_df)=c("Scenario peso","Scenario","Country","PDI_class","Myco_class","Weight","Mean","low_conf_int","Upper_conf_int")
names(res_pdi_weibull_F_df)=c("Scenario peso","Scenario","Country","PDI_class","Myco_class","Weight","Mean","low_conf_int","Upper_conf_int")
names(res_pdi_norm_F_df)=c("Scenario peso","Scenario","Country","PDI_class","Myco_class","Weight","Mean","low_conf_int","Upper_conf_int")

res=list(res_pdi_exp_F_df,
         res_pdi_weibull_F_df,
         res_pdi_norm_F_df)
return(res)

}

########################################################################################################################à
mean_only_scen=calculate_scenario(data_PDI[1:3,],scen_A)
mean_scen_Adults=calculate_scenario(data_PDI[4:8,],scen_A)
mean_UB_Adults=calculate_scenario(data_PDI[11:18,],scen_A)
mean_LB_Adults=calculate_scenario(data_PDI[21:28,],scen_A)

mean_female=calculate_scenario(data_PDI[9,],scen_F)
mean_male=calculate_scenario(data_PDI[c(10),],scen_M)

UB_female=calculate_scenario(data_PDI[19,],scen_F)
UB_male=calculate_scenario(data_PDI[20,],scen_M)

LB_female=calculate_scenario(data_PDI[29,],scen_F)
LB_male=calculate_scenario(data_PDI[30,],scen_M)

file.remove("PDI_exp_results_EU.xls")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",Wdf,"Scenari")

i=1

XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",mean_only_scen[[i]],"mean_only")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",mean_scen_Adults[[i]],"mean_scen_Adults")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",mean_UB_Adults[[i]],"mean_UB_Adults")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",mean_LB_Adults[[i]],"mean_LB_Adults")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",mean_male[[i]],"mean_male")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",mean_female[[i]],"mean_female")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",UB_male[[i]],"UB_male")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",UB_female[[i]],"UB_female")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",LB_male[[i]],"LB_male")
XLConnect::writeWorksheetToFile("PDI_exp_results_EU.xls",LB_female[[i]],"LB_female")

file.remove("PDI_weibull_results_EU.xls")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",Wdf,"Scenari")

i=2
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",mean_only_scen[[i]],"mean_only")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",mean_scen_Adults[[i]],"mean_scen_Adults")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",mean_UB_Adults[[i]],"mean_UB_Adults")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",mean_LB_Adults[[i]],"mean_LB_Adults")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",mean_male[[i]],"mean_male")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",mean_female[[i]],"mean_female")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",UB_male[[i]],"UB_male")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",UB_female[[i]],"UB_female")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",LB_male[[i]],"LB_male")
XLConnect::writeWorksheetToFile("PDI_weibull_results_EU.xls",LB_female[[i]],"LB_female")

file.remove("PDI_norm_results_EU.xls")

i=3
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",Wdf,"Scenari")

XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",mean_only_scen[[i]],"mean_only")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",mean_scen_Adults[[i]],"mean_scen_Adults")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",mean_UB_Adults[[i]],"mean_UB_Adults")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",mean_LB_Adults[[i]],"mean_LB_Adults")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",mean_male[[i]],"mean_male")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",mean_female[[i]],"mean_female")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",UB_male[[i]],"UB_male")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",UB_female[[i]],"UB_female")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",LB_male[[i]],"LB_male")
XLConnect::writeWorksheetToFile("PDI_norm_results_EU.xls",LB_female[[i]],"LB_female")


##########################################################################################################################

###################
# references

#weibull mean 
# https://stats.stackexchange.com/questions/361191/attempting-to-find-mean-of-weibull-function-in-r
# https://stackoverflow.com/questions/11817883/fitting-a-3-parameter-weibull-distribution
