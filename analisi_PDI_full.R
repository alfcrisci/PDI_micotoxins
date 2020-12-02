library(readxl)
library(doBy)
library(XLConnect)

###########################################################################################################
# functions

pdi_func=function(x,escr=70,vol_urine=2,weight=70) {return(x*(vol_urine/weight)*(100/escr))}
sumfit=function(x) {return(data.frame(par=x$estimate[1],errpar=x$estimate[2],aic=x$aic,names=x$distname))}

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

data_weigth=as.data.frame(read_xls("final_data/biometria_EU.xls",1))
data_weigth$Country=as.character(data_weigth$Country)
data_weigth_ITNOUK=data_weigth[which(data_weigth$Country %in% c("IT","NO","UK")==T),]
data_weigth_ITNOUK=data_weigth_ITNOUK[which(data_weigth$Sex %in% c("Male","Female")==T),]
data_weigth_ITNOUK$Weight=as.numeric(data_weigth_ITNOUK$Weight)


# https://efsa.onlinelibrary.wiley.com/doi/pdf/10.2903/j.efsa.2012.2579

sL=function(x) {as.numeric(quantile(na.omit(x),c(0.05)))}
sM=function(x) {as.numeric(quantile(na.omit(x),c(0.5)))}
sU=function(x) {as.numeric(quantile(na.omit(x),c(0.95)))}

Wdf=summaryBy(Weight ~ Sample_type+Sex+Country, data=data_weigth_ITNOUK, FUN=c(sL,sM,sU))[7:12,2:6]
Wdf_full=summaryBy(Weight ~ Sample_type+Country, data=data_weigth_ITNOUK, FUN=c(sL,sM,sU))[4:6,2:5]

Wdf=rbind(Wdf,data.frame(Sex=c("Adults"),Wdf_full))
Wdf=rbind(Wdf,c("Female","EU",50,67.2,90.7)) # efsa 
Wdf=rbind(Wdf,c("Male","EU",63.0,82.0,105)) # efsa 
Wdf=rbind(Wdf,c("Adults","EU",52.0,73.9,100.0)) # efsa 

scen=c(rep("Lower",12),rep("Medium",12),rep("Upper",12))
Wkg=as.numeric(c(Wdf[,3],Wdf[,4],Wdf[,5])) # 
cat=c(Wdf[,1],Wdf[,1],Wdf[,1])
mycos=c(Wdf[,2],Wdf[,2],Wdf[,2])


Wdf=data.frame(cat,
      scen,
      mycos,
      Wkg,stringsAsFactors = F)

write.csv(Wdf,file="Scenarios_full.csv",row.names = F)

scen_F=Wdf[Wdf$cat=="Female",]
scen_M=Wdf[Wdf$cat=="Male",]
scen_A=Wdf[Wdf$cat=="Adults",]

##################################################################################################
data_PDI=as.data.frame(read_xls("final_data/PDI_params_table.xls",1))
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
                    FUN=function(x){ pdi_func(x,escr = data_PDI$exrate[i],weight=Wkg[j])})

temp_weibull_PDI=sapply(rweibull(10000, data_PDI$weibull_shape_mean[i], data_PDI$weibull_scale_mean[i]),
                        FUN=function(x){ pdi_func(x,escr = data_PDI$exrate[i],weight=Wkg[j])})

temp_norm_PDI=sapply(rnorm(10000, data_PDI$norm_mean_mean[i], data_PDI$norm_mean_sd[i]),
                        FUN=function(x){ pdi_func(x,escr = data_PDI$exrate[i],weight=Wkg[j])})


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
mean_UB_Adults=calculate_scenario(data_PDI[13:17,],scen_A)
mean_LB_Adults=calculate_scenario(data_PDI[22:26,],scen_A)

mean_male=calculate_scenario(data_PDI[c(10,12),],scen_M)
mean_female=calculate_scenario(data_PDI[9:11,],scen_F)

UB_male=calculate_scenario(data_PDI[19:21,],scen_M)
UB_female=calculate_scenario(data_PDI[18:20,],scen_F)

LB_male=calculate_scenario(data_PDI[28:30,],scen_M)
LB_female=calculate_scenario(data_PDI[27:29,],scen_F)

file.remove("PDI_exp_results.xls")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",Wdf,"Scenari")

i=1

XLConnect::writeWorksheetToFile("PDI_exp_results.xls",mean_only_scen[[i]],"mean_only")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",mean_scen_Adults[[i]],"mean_scen_Adults")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",mean_UB_Adults[[i]],"mean_UB_Adults")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",mean_LB_Adults[[i]],"mean_LB_Adults")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",mean_male[[i]],"mean_male")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",mean_female[[i]],"mean_female")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",UB_male[[i]],"UB_male")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",UB_female[[i]],"UB_female")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",LB_male[[i]],"LB_male")
XLConnect::writeWorksheetToFile("PDI_exp_results.xls",LB_female[[i]],"LB_female")

file.remove("PDI_weibull_results.xls")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",Wdf,"Scenari")

i=2
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",mean_only_scen[[i]],"mean_only")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",mean_scen_Adults[[i]],"mean_scen_Adults")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",mean_UB_Adults[[i]],"mean_UB_Adults")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",mean_LB_Adults[[i]],"mean_LB_Adults")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",mean_male[[i]],"mean_male")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",mean_female[[i]],"mean_female")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",UB_male[[i]],"UB_male")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",UB_female[[i]],"UB_female")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",LB_male[[i]],"LB_male")
XLConnect::writeWorksheetToFile("PDI_weibull_results.xls",LB_female[[i]],"LB_female")

file.remove("PDI_norm_results.xls")

i=3
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",Wdf,"Scenari")

XLConnect::writeWorksheetToFile("PDI_norm_results.xls",mean_only_scen[[i]],"mean_only")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",mean_scen_Adults[[i]],"mean_scen_Adults")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",mean_UB_Adults[[i]],"mean_UB_Adults")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",mean_LB_Adults[[i]],"mean_LB_Adults")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",mean_male[[i]],"mean_male")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",mean_female[[i]],"mean_female")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",UB_male[[i]],"UB_male")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",UB_female[[i]],"UB_female")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",LB_male[[i]],"LB_male")
XLConnect::writeWorksheetToFile("PDI_norm_results.xls",LB_female[[i]],"LB_female")




###################
# References

# weibull mean 
# https://stats.stackexchange.com/questions/361191/attempting-to-find-mean-of-weibull-function-in-r
# https://stackoverflow.com/questions/11817883/fitting-a-3-parameter-weibull-distribution
