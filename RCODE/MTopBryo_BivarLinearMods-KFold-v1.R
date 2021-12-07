

library(readxl)
library(raster)
library(dplyr)
library(tidyr)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(ggcorrplot)
library(caret)


## ------------------------------------------------------------------------------- ##


nicheVarsPath <- "./OUT/nicheVars_v4.xlsx"

nv_ <- read_excel(nicheVarsPath) %>% as.data.frame()

# pLoss - scenario he45: Linear models (between niche metrics, spatial metrics and loss of area)
# 
# nv <- nv_ %>% mutate_at(.vars=c(-48,-49,-52,-53,-56,-57,-60,-61),  
#                                           .funs = list(sc=scale))


nv <- data.frame(
  nv_[,1:2],
  scale(nv_[,c(3:44, 65:68)]), # 
  nv_[,45:64] # Response
)




print(colnames(nv))


## ------------------------------------------------------------------------------- ##


respVars <- c("pLoss_he45", 
              "pLoss_he85", 
              "pLoss_mp45",
              "pLoss_mp85")


# 
# predVars <- c(
#   niche_breath = "hv_svm_log",
#   niche_pos    = "avgSpDistCentr",
#   niche_ovlp   = "avg_ovlp_jacc",
#   
#   spatial_range = "suitArea_TSS",
#   spatial_pos   = "spPosLog",
#   spatial_ovlp  = "avgOvlpJacc_SDM"
# )

predVars <- c(
  
  # Simple models
  smp_niche_breath  = "hv_svm_log",
  smp_niche_pos     = "distToAvgCentr",
  smp_niche_ovlp    = "avg_ovlp_jacc",
  smp_spatial_range = "suitArea_TSS",
  smp_spatial_pos   = "wspPosLog",
  smp_spatial_ovlp  = "avgOvlpJacc_SDM",
  
  # Additive models
  add_niche_all     = "hv_svm_log + distToAvgCentr + avg_ovlp_jacc",
  add_niche_brd_pos = "hv_svm_log + distToAvgCentr",
  add_spatial_all   = "suitArea_TSS + wspPosLog + avgOvlpJacc_SDM",
  
  # Interaction
  int_niche_brd_pos = "hv_svm_log * distToAvgCentr"
  
)


## ------------------------------------------------------------------------------- ##


predVarNames <- data.frame(
  
  vname = predVars,
  vtype = names(predVars)
)

#set.seed(123456)

nfolds <- 5

rnds <- 100

indVecs <- sample(1:nrow(nv))

testGroups <- split(indVecs, 1:nfolds)


k <- 0

nt <- rnds * length(respVars) * length(predVars) * nfolds

pb <- txtProgressBar(1,nt,style=3)


## ------------------------------------------------------------------------------- ##


for(r in 1:rnds){
  
  indVecs <- sample(1:nrow(nv))
  
  testGroups <- split(indVecs, 1:nfolds)
  
  
  for(respVar in respVars){
    
    for(predVar in predVars){
      
      form <- as.formula(paste(respVar,"~",predVar,sep=""))
      #print(form)
      
      for(i in 1:nfolds){
        
        k <- k+1
        
        testIdx <- testGroups[[i]]
        
        testDF <- nv[testIdx, ]
        trainDF <- nv[-testIdx, ]
        
        lmod <- lm(form, data = nv)
        
        obsData <- as.data.frame(testDF)[,respVar]
        
        predData <- predict(lmod, 
                            newdata = testDF, 
                            type = "response")
        
        r2 <- cor(x = obsData, predData)^2
        N <- nrow(nv)
        p <- length(coefficients(lmod))-1
        r2_adj <- ((1 - r2) * (N - 1)) / (N - p - 1)
        
        tmp <- data.frame(
          respVar = respVar,
          predVar = predVar,
          cvRnd   = r,
          itRnd   = i,
          R2      = r2,
          R2_adj  = r2_adj,
          corSp   = cor(obsData, predData, method="spearman"),
          RMSE    = sqrt(mean((obsData - predData)^2)),
          NRMSE   = (sqrt(mean((obsData - predData)^2))) / mean(obsData),
          MAD     = median(abs(obsData - predData))
        )
        
        if(k==1){
          res <- tmp
        }else{
          res <- rbind(res, tmp)
        }
        
        #cat(r,respVar,predVar,i,"\n\n",sep = " | ")
        setTxtProgressBar(pb, k)
      }
    }
  }
}


## ------------------------------------------------------------------------------- ##


mergedRes <- 

res %>% 
  select(-itRnd) %>% 
  group_by(cvRnd,respVar,predVar) %>% 
  summarise_all(.funs=list(mean)) %>% 
  ungroup() %>% 
  select(-cvRnd) %>% 
  group_by(respVar,predVar) %>% 
  summarise_all(.funs=list(avg=mean)) %>% 
  ungroup %>% 
  left_join(predVarNames, by = c("predVar" = "vname")) %>% 
  select(1,2,9,3:8) %>% 
  arrange(respVar, R2_adj_avg)


write.csv(mergedRes,"./OUT/merged_KFold-CV_pLoss-vs-Niche_Sp_metrics-v2.csv", 
          row.names = FALSE)

