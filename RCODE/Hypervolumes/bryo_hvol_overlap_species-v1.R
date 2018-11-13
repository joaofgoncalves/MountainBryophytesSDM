
## --------------------------------------------------------- ##
##
## Performs an overlap analysis between species hypervolumes
##
## João Gonçalves & Helena Hespanhol
## CIBIO/InBIO, FCUP
## Porto, 10/2018
##
## --------------------------------------------------------- ##

library(hypervolume)
library(raster)
library(dplyr)
library(magrittr)
library(rgdal)
library(ggplot2)
library(ggdendro)

# Wrapper function to calculate overlap statistics
hv_overlap <- function(hv1, hv2, verbose=FALSE, ...){
  hv_set <- hypervolume_set(hv1, hv2, check.memory = FALSE, verbose=verbose, ...)
  hv_ovlp_stats <- hypervolume_overlap_statistics(hv_set)
  return(hv_ovlp_stats)
}

# ------------------------------------------------------------------------- #
# Load data ----
# ------------------------------------------------------------------------- #

# Load hypervolume objects from previous analyses 
load("./OUT/HyperVolumeBySpecies-v3-20181101.RData")

# Load shapefile
spData <- readOGR("./DATA/VECTOR/Bryophyte_dataset","And_Gri_Rac_PI_all_2")


# ------------------------------------------------------------------------- #
# Perform overlap analysis by species pairs (lower-tri matrix only) ----
# ------------------------------------------------------------------------- #


spCodesAll <- unique(spData$Cod_esp)
len <- length(spCodesAll)

ltrimat <- matrix(1:len^2,len,len) %>% lower.tri
ovlp_jacc <- matrix(NA,len,len,dimnames = list(spCodesAll,spCodesAll))
ovlp_sors <- matrix(NA,len,len,dimnames = list(spCodesAll,spCodesAll))

tot2run <- sum(ltrimat)
pb <- txtProgressBar(1,tot2run,style = 3)

k <- 0
for(i in 1:len){
  for(j in 1:len){
    
    if(!ltrimat[i,j]){
      next 
    }
    
    k <- k+1
    sp1 <- as.character(spCodesAll[i])
    sp2 <- as.character(spCodesAll[j])
    
    hv1 <- hvObj_BySpecies[[paste("hv_svm_",sp1,sep="")]]
    hv2 <- hvObj_BySpecies[[paste("hv_svm_",sp2,sep="")]]
    hv_ovlp_ind <- hv_overlap(hv1, hv2)
    
    ovlp_jacc[i,j] <- hv_ovlp_ind[1]
    ovlp_sors[i,j] <- hv_ovlp_ind[2]
    
    setTxtProgressBar(pb, k)
  }
}

save(ovlp_jacc, ovlp_sors, file = "./OUT/NicheOvlpDistances-NewVars_v3.RData")


# --------------------------------------------------------------------------- #
# Make dendrogram of species hiche overlap distances ----
# --------------------------------------------------------------------------- #


hc_jacc <- hclust((1-as.dist(ovlp_jacc)), method="complete")
hc_sors <- hclust((1-as.dist(ovlp_sors)), method="complete")

#plot(hc_jacc, horiz=TRUE, hang=-1)
ggd <- ggdendrogram(hc_jacc, rotate = TRUE, size = 3) +
  labs(title="Dendrogram of species niche overlap",
       subtitle="Jaccard distance between hypervolumes")

ggsave("./OUT/DendroSpeciesNicheOvlp_Jacc-NewVars_v3.png",ggd,height = 7, width=9)

ggd <- ggdendrogram(hc_sors, rotate = TRUE, size = 3) +
  labs(title="Dendrogram of species niche overlap",
       subtitle="Sorensen distance between hypervolumes")

ggsave("./OUT/DendroSpeciesNicheOvlp_Sors-NewVars_v3.png",ggd,height = 7, width=9)


