
library (AICcmodavg)
library (xlsx)

#Read table

bryophyte<-read.csv("nicheVars_v4.csv")

# tabela só com valores numéricos

bryo_select <- select(bryophyte, 3:68)
bryo_new <- na.omit(bryo_select)

#scale only predictors, excluding response variable

bryo_scale_pred <- bryo_new %>% mutate_at(.vars=c(-48,-49,-52,-53,-56,-57,-60,-61),  .funs = list(sc=scale))# pLoss - scenario he45: Linear models (between niche metrics, spatial metrics and loss of area)

# pLoss - scenario he45: Linear models (between niche metrics, spatial metrics and loss of area)

#simple models

mod_Lhe45_null<- lm (pLoss_he45~1,data=bryo_scale_pred)
summary(mod_Lhe45_null)
par(mfrow=c(2,2))
AIC(mod_Lhe45_null)

mod_Lhe45_lm_br <- lm (pLoss_he45~hv_svm_log_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_br)
AIC(mod_Lhe45_lm_br)
par(mfrow=c(2,2))
plot(mod_Lhe45_lm_br)

mod_Lhe45_lm_pa <- lm (pLoss_he45~distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_pa)
AIC (mod_Lhe45_lm_pa)
plot(mod_Lhe45_lm_pa)

mod_Lhe45_lm_no <- lm (pLoss_he45~avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_no)
AIC (mod_Lhe45_lm_no)
plot (mod_Lhe45_lm_no)

mod_Lhe45_lm_rs <- lm (pLoss_he45~suitArea_TSS_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_rs)
AIC (mod_Lhe45_lm_rs)
plot(mod_Lhe45_lm_rs)

mod_Lhe45_lm_sp <- lm (pLoss_he45~wspPosLog_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_sp)
AIC (mod_Lhe45_lm_sp)
plot(mod_Lhe45_lm_sp)

mod_Lhe45_lm_so <- lm (pLoss_he45~avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_so)
AIC (mod_Lhe45_lm_so)
plot(mod_Lhe45_lm_so)

#additive models

mod_Lhe45_lm_br_pa_no <- lm (pLoss_he45~hv_svm_log_sc + distToAvgCentr_sc + avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_br_pa_no)
AIC (mod_Lhe45_lm_br_pa_no)

mod_Lhe45_lm_br_pa <- lm (pLoss_he45~hv_svm_log_sc + distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_br_pa)
AIC (mod_Lhe45_lm_br_pa)

mod_Lhe45_lm_rs_sp_so <- lm (pLoss_he45~suitArea_TSS_sc + wspPosLog_sc + avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_rs_sp_so)
AIC (mod_Lhe45_lm_rs_sp_so)

# interaction models

mod_Lhe45_lm_br_pa_int <- lm (pLoss_he45~hv_svm_log_sc * distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lhe45_lm_br_pa_int)
AIC (mod_Lhe45_lm_br_pa_int)


# pLoss - scenario he85: Linear models (between niche metrics, spatial metrics and loss of area)

#simple models

mod_Lhe85_null<- lm (pLoss_he85~1,data=bryo_scale_pred)
summary(mod_Lhe85_null)
AIC(mod_Lhe85_null)

mod_Lhe85_lm_br <- lm (pLoss_he85~hv_svm_log_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_br)
AIC(mod_Lhe85_lm_br)
plot(mod_Lhe85_lm_br)

mod_Lhe85_lm_pa <- lm (pLoss_he85~distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_pa)
AIC (mod_Lhe85_lm_pa)
plot(mod_Lhe85_lm_pa)

mod_Lhe85_lm_no <- lm (pLoss_he85~avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_no)
AIC (mod_Lhe85_lm_no)
plot(mod_Lhe85_lm_no)

mod_Lhe85_lm_rs <- lm (pLoss_he85~suitArea_TSS_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_rs)
AIC (mod_Lhe85_lm_rs)
plot(mod_Lhe85_lm_rs)

mod_Lhe85_lm_sp <- lm (pLoss_he85~wspPosLog_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_sp)
AIC (mod_Lhe85_lm_sp)
plot(mod_Lhe85_lm_sp)

mod_Lhe85_lm_so <- lm (pLoss_he85~avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_so)
AIC (mod_Lhe85_lm_so)
plot(mod_Lhe85_lm_so)

#additive models

mod_Lhe85_lm_br_pa_no <- lm (pLoss_he85~hv_svm_log_sc + distToAvgCentr_sc + avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_br_pa_no)
AIC (mod_Lhe85_lm_br_pa_no)

mod_Lhe85_lm_br_pa <- lm (pLoss_he85~hv_svm_log_sc + distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_br_pa)
AIC (mod_Lhe85_lm_br_pa)

mod_Lhe85_lm_rs_sp_so <- lm (pLoss_he85~suitArea_TSS_sc + wspPosLog_sc + avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_rs_sp_so)
AIC (mod_Lhe85_lm_rs_sp_so)

# interaction models

mod_Lhe85_lm_br_pa_int <- lm (pLoss_he85~hv_svm_log_sc * distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lhe85_lm_br_pa_int)
AIC (mod_Lhe85_lm_br_pa_int)


# pLoss - scenario mp45:Linear models (between niche metrics, spatial metrics and loss of area)


#simple models

mod_Lmp45_null<- lm (pLoss_mp45~1,data=bryo_scale_pred)
summary(mod_Lmp45_null)
AIC(mod_Lmp45_null)

mod_Lmp45_lm_br <- lm (pLoss_mp45~hv_svm_log_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_br)
AIC(mod_Lmp45_lm_br)
plot(mod_Lmp45_lm_br)

mod_Lmp45_lm_pa <- lm (pLoss_mp45~distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_pa)
AIC (mod_Lmp45_lm_pa)
plot(mod_Lmp45_lm_pa)

mod_Lmp45_lm_no <- lm (pLoss_mp45~avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_no)
AIC (mod_Lmp45_lm_no)
plot(mod_Lmp45_lm_no)

mod_Lmp45_lm_rs <- lm (pLoss_mp45~suitArea_TSS_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_rs)
AIC (mod_Lmp45_lm_rs)
plot(mod_Lmp45_lm_rs)

mod_Lmp45_lm_sp <- lm (pLoss_mp45~wspPosLog_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_sp)
AIC (mod_Lmp45_lm_sp)
plot(mod_Lmp45_lm_sp)

mod_Lmp45_lm_so <- lm (pLoss_mp45~avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_so)
AIC (mod_Lmp45_lm_so)
plot(mod_Lmp45_lm_so)

#additive models

mod_Lmp45_lm_br_pa_no <- lm (pLoss_mp45~hv_svm_log_sc + distToAvgCentr_sc + avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_br_pa_no)
AIC (mod_Lmp45_lm_br_pa_no)

mod_Lmp45_lm_br_pa <- lm (pLoss_mp45~hv_svm_log_sc + distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_br_pa)
AIC (mod_Lmp45_lm_br_pa)

mod_Lmp45_lm_rs_sp_so <- lm (pLoss_mp45~suitArea_TSS_sc + wspPosLog_sc + avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_rs_sp_so)
AIC (mod_Lmp45_lm_rs_sp_so)

# interaction models

mod_Lmp45_lm_br_pa_int <- lm (pLoss_mp45~hv_svm_log_sc * distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lmp45_lm_br_pa_int)
AIC (mod_Lmp45_lm_br_pa_int)

# pLoss - scenario mp85:Linear models (between niche metrics, spatial metrics and loss of area) 

#simple models

mod_Lmp85_null<- lm (pLoss_mp85~1,data=bryo_scale_pred)
summary(mod_Lmp85_null)
AIC(mod_Lmp85_null)


mod_Lmp85_lm_br <- lm (pLoss_mp85~hv_svm_log_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_br)
AIC(mod_Lmp85_lm_br)
plot(mod_Lmp85_lm_br)

mod_Lmp85_lm_pa <- lm (pLoss_mp85~distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_pa)
AIC (mod_Lmp85_lm_pa)
plot(mod_Lmp85_lm_pa)

mod_Lmp85_lm_no <- lm (pLoss_mp85~avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_no)
AIC (mod_Lmp85_lm_no)
plot(mod_Lmp85_lm_no)

mod_Lmp85_lm_rs <- lm (pLoss_mp85~suitArea_TSS_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_rs)
AIC (mod_Lmp85_lm_rs)
plot(mod_Lmp85_lm_rs)

mod_Lmp85_lm_sp <- lm (pLoss_mp85~wspPosLog_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_sp)
AIC (mod_Lmp85_lm_sp)
plot(mod_Lmp85_lm_sp)

mod_Lmp85_lm_so <- lm (pLoss_mp85~avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_so)
AIC (mod_Lmp85_lm_so)
plot(mod_Lmp85_lm_so)

#additive models

mod_Lmp85_lm_br_pa_no <- lm (pLoss_mp85~hv_svm_log_sc + distToAvgCentr_sc + avg_ovlp_jacc_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_br_pa_no)
AIC (mod_Lmp85_lm_br_pa_no)

mod_Lmp85_lm_br_pa <- lm (pLoss_mp85~hv_svm_log_sc + distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_br_pa)
AIC (mod_Lmp85_lm_br_pa)

mod_Lmp85_lm_rs_sp_so <- lm (pLoss_mp85~suitArea_TSS_sc + wspPosLog_sc + avgOvlpJacc_SDM_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_rs_sp_so)
AIC (mod_Lmp85_lm_rs_sp_so)

# interaction models

mod_Lmp85_lm_br_pa_int <- lm (pLoss_mp85~hv_svm_log_sc * distToAvgCentr_sc,data=bryo_scale_pred)
summary(mod_Lmp85_lm_br_pa_int)
AIC (mod_Lmp85_lm_br_pa_int)


# Multimodel inference: pLoss - scenario he45: Linear models (between niche metrics, spatial metrics and loss of area) 


cand.mod <- list(mod_Lhe45_lm_br, mod_Lhe45_lm_pa, mod_Lhe45_lm_no, mod_Lhe45_lm_rs, mod_Lhe45_lm_sp, mod_Lhe45_lm_so,mod_Lhe45_lm_br_pa_no, mod_Lhe45_lm_rs_sp_so, mod_Lhe45_lm_br_pa, mod_Lhe45_lm_br_pa_int, mod_Lhe45_null)
Modnames <- c("mod_Lhe45_lm_br", "mod_Lhe45_lm_pa","mod_Lhe45_lm_no", "mod_Lhe45_lm_rs", "mod_Lhe45_lm_sp", "mod_Lhe45_lm_so", "mod_Lhe45_lm_br_pa_no", "mod_Lhe45_lm_rs_sp_so","mod_Lhe45_lm_br_pa","mod_Lhe45_lm_br_pa_int", "mod_Lhe45_null")

MMI_AIC_tab_pLoss_he45<- aictab(cand.set=cand.mod, modnames = Modnames, second.ord=TRUE, sort =TRUE)

MMI_AIC_tab_pLoss_he45.table<-as.data.frame(MMI_AIC_tab_pLoss_he45)
file<- paste ("C:/Users/Helena/Documents/Helena/Software_Rstudio/mountain_bryophytes/regression/tables/", "AIC_pLoss_he45_scale.xlsx")
write.xlsx(MMI_AIC_tab_pLoss_he45.table, file)

# Multimodel inference: pLoss - scenario he85: Linear models (between niche metrics, spatial metrics and loss of area) 


cand.mod <- list(mod_Lhe85_lm_br, mod_Lhe85_lm_pa, mod_Lhe85_lm_no, mod_Lhe85_lm_rs, mod_Lhe85_lm_sp, mod_Lhe85_lm_so, mod_Lhe85_lm_br_pa_no, mod_Lhe85_lm_rs_sp_so, mod_Lhe85_lm_br_pa, mod_Lhe85_lm_br_pa_int,mod_Lhe85_null)
Modnames <- c("mod_Lhe85_lm_br", "mod_Lhe85_lm_pa", "mod_Lhe85_lm_no", "mod_Lhe85_lm_rs", "mod_Lhe85_lm_sp", "mod_Lhe85_lm_so", "mod_Lhe85_lm_br_pa_no", "mod_Lhe85_lm_rs_sp_so", "mod_Lhe85_lm_br_pa","mod_Lhe85_lm_br_pa_int", "mod_Lhe85_null")
MMI_AIC_tab_pLoss_he85<- aictab(cand.set=cand.mod, modnames = Modnames, second.ord=TRUE, sort =TRUE)

MMI_AIC_tab_pLoss_he85.table<-as.data.frame(MMI_AIC_tab_pLoss_he85)
file<- paste ("C:/Users/Helena/Documents/Helena/Software_Rstudio/mountain_bryophytes/regression/tables/", "AIC_pLoss_he85_scale.xlsx")
write.xlsx(MMI_AIC_tab_pLoss_he85.table, file)


# Multimodel inference: pLoss - scenario mp45: Linear models (between niche metrics, spatial metrics and loss of area) 


cand.mod <- list(mod_Lmp45_lm_br, mod_Lmp45_lm_pa, mod_Lmp45_lm_no, mod_Lmp45_lm_rs, mod_Lmp45_lm_sp, mod_Lmp45_lm_so, mod_Lmp45_lm_br_pa_no, mod_Lmp45_lm_rs_sp_so, mod_Lmp45_lm_br_pa, mod_Lmp45_lm_br_pa_int, mod_Lmp45_null)
Modnames <- c("mod_Lmp45_lm_br", "mod_Lmp45_lm_pa", "mod_Lmp45_lm_no", "mod_Lmp45_lm_rs", "mod_Lmp45_lm_sp", "mod_Lmp45_lm_so", "mod_Lmp45_lm_br_pa_no", "mod_Lmp45_lm_rs_sp_so", "mod_Lmp45_lm_br_pa", "mod_Lmp45_lm_br_pa_int", "mod_Lmp45_null")
MMI_AIC_tab_pLoss_mp45<- aictab(cand.set=cand.mod, modnames = Modnames, second.ord=TRUE, sort =TRUE)

MMI_AIC_tab_pLoss_mp45.table<-as.data.frame(MMI_AIC_tab_pLoss_mp45)
file<- paste ("C:/Users/Helena/Documents/Helena/Software_Rstudio/mountain_bryophytes/regression/tables/", "AIC_pLoss_mp45_scale.xlsx")
write.xlsx(MMI_AIC_tab_pLoss_mp45.table, file)

# Multimodel inference: pLoss - scenario mp85: Linear models (between niche metrics, spatial metrics and loss of area) 


cand.mod <- list(mod_Lmp85_lm_br, mod_Lmp85_lm_pa, mod_Lmp85_lm_no, mod_Lmp85_lm_rs, mod_Lmp85_lm_sp, mod_Lmp85_lm_so, mod_Lmp85_lm_br_pa_no, mod_Lmp85_lm_rs_sp_so, mod_Lmp85_lm_br_pa, mod_Lmp85_lm_br_pa_int, mod_Lmp85_null)
Modnames <- c("mod_Lmp85_lm_br", "mod_Lmp85_lm_pa", "mod_Lmp85_lm_no", "mod_Lmp85_lm_rs", "mod_Lmp85_lm_sp", "mod_Lmp85_lm_so", "mod_Lmp85_lm_br_pa_no", "mod_Lmp85_lm_rs_sp_so", "mod_Lmp85_lm_br_pa","mod_Lmp85_lm_br_pa_int","mod_Lmp85_null")
MMI_AIC_tab_pLoss_mp85<- aictab(cand.set=cand.mod, modnames = Modnames, second.ord=TRUE, sort =TRUE)

MMI_AIC_tab_pLoss_mp85.table<-as.data.frame(MMI_AIC_tab_pLoss_mp85)
file<- paste ("C:/Users/Helena/Documents/Helena/Software_Rstudio/mountain_bryophytes/regression/tables/", "AIC_pLoss_mp85_scale.xlsx")
write.xlsx(MMI_AIC_tab_pLoss_mp85.table, file)
