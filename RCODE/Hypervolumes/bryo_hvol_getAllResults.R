
## --------------------------------------------------------- ##
##
## Aggregate all results
##
## João Gonçalves & Helena Hespanhol
## CIBIO/InBIO, FCUP
## Porto, 10/2018
##
## --------------------------------------------------------- ##


library(dplyr)
library(hypervolume)


# Load hypervolume objects from previous analyses -----
load("./OUT/HyperVolumeBySpecies-v3-20181101.RData")

nicheBreadthBySp <- data.frame(spName=spCodesAll,
           hv_box = hv_box_vols,
           hv_gau = hv_gauss_vols,
           hv_svm = hv_svm_vols,
           hv_box_log = log10(hv_box_vols),
           hv_gau_log = log10(hv_gauss_vols),
           hv_svm_log = log10(hv_svm_vols),
           stringsAsFactors = FALSE)


# ----------------------------------------------------------------------------------- #

# Load hypervolume overlap results
load("./OUT/NicheOvlpDistances-NewVars_v3.RData")

## Analyze hypervolume results ----
## Niche position/breadth vs. niche overlap

# Create full distance matrix based on Jaccard similarity of niches
dd_jacc <- as.dist(ovlp_jacc, upper=TRUE)
ddmat_jacc <- as.matrix(dd_jacc)
diag(ddmat_jacc) <- NA

dd_sors <- as.dist(ovlp_sors, upper=TRUE)
ddmat_sors <- as.matrix(dd_sors)
diag(ddmat_sors) <- NA


spCodesAll <- attr(dd_jacc,"Labels")

# Calculate the average overlap and get Hvolume
nicheOvlpBySp <- data.frame(spName = spCodesAll,
                 avg_ovlp_jacc = apply(ddmat_jacc, 1, mean, na.rm=TRUE), 
                 std_ovlp_jacc = apply(ddmat_jacc, 1, sd, na.rm=TRUE),
                 
                 avg_ovlp_sors = apply(ddmat_sors, 1, mean, na.rm=TRUE), 
                 std_ovlp_sors = apply(ddmat_sors, 1, sd, na.rm=TRUE),
                 stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------- #

rstBinDF <- readRDS("./OUT/rstBinDF_TSS_v3.RData")

distBinJacc <- ade4::dist.binary(t(rstBinDF), diag = FALSE, method = 1)
distBinJacc <- 1 - as.matrix(distBinJacc)
diag(distBinJacc) <- NA

distBinSors <- ade4::dist.binary(t(rstBinDF), diag = FALSE, method = 5)
distBinSors <- 1 - as.matrix(distBinSors)
diag(distBinSors) <- NA

spCodesAll <- colnames(distBinJacc)

# Calculate the average overlap and get Hvolume
nicheOvlpBySp_SDM <- data.frame(spName = spCodesAll,
                            avgOvlpJacc_SDM = apply(distBinJacc, 1, mean, na.rm=TRUE), 
                            stdOvlpJacc_SDM = apply(distBinJacc, 1, sd, na.rm=TRUE),
                            
                            avgOvlpSors_SDM = apply(distBinSors, 1, mean, na.rm=TRUE), 
                            stdOvlpSors_SDM = apply(distBinSors, 1, sd, na.rm=TRUE),
                            stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------- #

dismoIndI <- readRDS("./OUT/dismo_nicheOverlap_indI_v3.rds")
dismoIndD <- readRDS("./OUT/dismo_nicheOverlap_indD_v3.rds")

spCodesAll <- colnames(dismoIndI)

dismoIndI <- as.matrix(as.dist(dismoIndI,upper = TRUE, diag=FALSE))
diag(dismoIndI) <- NA

dismoIndD <- as.matrix(as.dist(dismoIndD,upper = TRUE, diag=FALSE))
diag(dismoIndD) <- NA

# Calculate the average overlap and get Hvolume
nicheOvlpBySp_dismoSDM <- data.frame(spName = spCodesAll,
                                avgOvlpInd_I = apply(dismoIndI, 1, mean, na.rm=TRUE), 
                                stdOvlpInd_I = apply(dismoIndI, 1, sd, na.rm=TRUE),
                                
                                avgOvlpInd_D = apply(dismoIndD, 1, mean, na.rm=TRUE), 
                                stdOvlpInd_D = apply(dismoIndD, 1, sd, na.rm=TRUE),
                                stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------- #

## Load env and species data -----
rstStack <- raster::stack("./DATA/RASTER/_VARS/current.tif")
names(rstStack) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")

# Load shapefile
spData <- readOGR("./DATA/VECTOR/Bryophyte_dataset","And_Gri_Rac_PI_all_2", 
                  stringsAsFactors = FALSE)

# Extract env data to points and remove NA's
# Scale the entire dataset across species
#
spDataVars <- data.frame(spCode = as.character(spData$Cod_esp), 
                         #scale(raster::extract(y = spData, x = rstStack)),
                         raster::extract(y = spData, x = rstStack),
                         stringsAsFactors = FALSE) %>% na.omit

aggEnvData <- spDataVars %>% 
  group_by(spCode) %>% 
  summarize_all(.funs = c("mean", "sd")) %>% 
  rename(spName=spCode)

spDataVars <- data.frame(spCode = as.character(spData$Cod_esp), 
                         scale(raster::extract(y = spData, x = rstStack)),
                         #raster::extract(y = spData, x = rstStack),
                         stringsAsFactors = FALSE) %>% na.omit

# ----------------------------------------------------------------------------------- #

# Get overall centroid for all the species records
pointCentroids <- apply(spDataVars %>% select(-1) %>% na.omit, 2, mean)

# Get centroids for the niche hypervolume of each species
for(sp in spCodesAll){
  pointCentroids <- rbind(pointCentroids,
                          hypervolume::get_centroid(hvObj_BySpecies[[name=paste("hv_svm",sp,sep="_")]]))
}

rownames(pointCentroids) <- c("AvgCenter",as.character(spCodesAll))
dd <- dist(pointCentroids) %>% as.matrix

spDistToAvg <- dd[-1,1]

# Keep only distance between species
dd1 <- dd[-1,-1]

for(i in 1:ncol(dd1)){
  
  tmp <- data.frame(
    spName = colnames(dd1)[i],
    avgSpDistCentr = mean(dd1[i,-i]),
    stdSpDistCentr = sd(dd1[i,-i]),
    Q05SpDistCentr = quantile(dd1[i,-i], probs=0.05),
    Q50SpDistCentr = quantile(dd1[i,-i], probs=0.50),
    Q95SpDistCentr = quantile(dd1[i,-i], probs=0.95),
    stringsAsFactors = FALSE)
  
  if(i==1){
    distToSpCentrs <- tmp
  }else{
    distToSpCentrs <- rbind(distToSpCentrs,tmp)
  }
}

distToSpCentrs <- distToSpCentrs %>% mutate(distToAvgCentr=spDistToAvg)

cor(distToSpCentrs[,-1],method="spearman") %>% round(2) 


# ----------------------------------------------------------------------------------- #

distribSize_current <- read.csv("./OUT/ProjectedAreaBySpecies_proj_current-v1.csv") %>% 
  rename(spName=spNames)

distribSize_he45bi50 <- read.csv("./OUT/ProjectedAreaBySpecies_proj_he45bi50-v1.csv") %>% 
  rename(spName=spNames)
distribSize_he85bi50 <- read.csv("./OUT/ProjectedAreaBySpecies_proj_he85bi50-v1.csv") %>% 
  rename(spName=spNames)
distribSize_mp45bi50 <- read.csv("./OUT/ProjectedAreaBySpecies_proj_mp45bi50-v1.csv") %>% 
  rename(spName=spNames)
distribSize_mp85bi50 <- read.csv("./OUT/ProjectedAreaBySpecies_proj_mp85bi50-v1.csv") %>% 
  rename(spName=spNames)

perChangeClim <- data.frame(spName=distribSize_current$spName,

pCh_he45bi50 = ((distribSize_he45bi50$suitArea_TSS - distribSize_current$suitArea_TSS)/
                  distribSize_current$suitArea_TSS)*100,

pCh_he85bi50 = ((distribSize_he85bi50$suitArea_TSS - distribSize_current$suitArea_TSS)/
                  distribSize_current$suitArea_TSS)*100,

pCh_mp45bi50 = ((distribSize_mp45bi50$suitArea_TSS - distribSize_current$suitArea_TSS)/
                  distribSize_current$suitArea_TSS)*100,

pCh_mp85bi50 = ((distribSize_mp85bi50$suitArea_TSS - distribSize_current$suitArea_TSS)/
                  distribSize_current$suitArea_TSS)*100
)


# Load relative habitat change data

relHabChange <- read.csv("./OUT/relHabitatChange_BySpScenario_v3.csv",stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------------- #

nicheVars <- 
aggEnvData %>% # Env data by species
  left_join(nicheBreadthBySp, by="spName") %>% # Niche size/breadth
  left_join(distToSpCentrs, by="spName") %>% # Niche position
  left_join(nicheOvlpBySp, by="spName") %>%  # Niche overlap
  left_join(nicheOvlpBySp_SDM, by="spName") %>%  # Niche overlap SDM Jacc/Sors
  left_join(nicheOvlpBySp_dismoSDM, by="spName") %>% # Niche overlap w/ dismo
  left_join(distribSize_current, by="spName") %>% 
  left_join(perChangeClim, by="spName") %>% 
  left_join(relHabChange, by="spName")


# Write data
write.csv(nicheVars,"./OUT/nicheVars_v3.csv", row.names = FALSE)

#cm <- cor(nicheVars[,-c(1:13)],method="spearman") %>% round(2) 
cm <- cor(nicheVars[,-1],method="spearman") %>% round(2) 

# Write data
write.csv(cm,"./OUT/nicheVars_SpearmanCorrelation_v3.csv", row.names = FALSE)

#cm <- cor(nicheVars[,-c(1:13)],method="spearman") %>% round(2) 
cmPr <- cor(nicheVars[,-1]) %>% round(2) 

# Write data
write.csv(cmPr,"./OUT/nicheVars_PearsonCorrelation_v3.csv", row.names = FALSE)


View(cm)

cm[33:48, ]

cm["netHab_he45",-c(13:48)] %>% sort
cm["pLoss_he45",-c(13:48)] %>% sort
cm["pStable_he45",-c(13:48)] %>% sort

cm["netHab_he85",-c(13:48)] %>% sort
cm["pLoss_he85",-c(13:48)] %>% sort
cm["pStable_he85",-c(13:48)] %>% sort

cm["netHab_mp45",-c(13:48)] %>% sort
cm["pLoss_mp45",-c(13:48)] %>% sort
cm["pStable_mp45",-c(13:48)] %>% sort

cm["netHab_mp85",-c(13:48)] %>% sort
cm["pLoss_mp85",-c(13:48)] %>% sort
cm["pStable_mp85",-c(13:48)] %>% sort

