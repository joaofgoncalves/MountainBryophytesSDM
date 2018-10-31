
## --------------------------------------------------------- ##
##
## Fit, calculate and plot hypervolumes
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

#rm(list=ls())

## ----------------------------------------------------------------------------------------------- ##
## Load env and species data -----
## ----------------------------------------------------------------------------------------------- ##

# Data paths to raster variables
climVarsPaths <- list.files("./DATA/RASTER/WorldClim_1km/a2000", pattern=".tif$", full.names = TRUE)
topoVarsPaths <- list.files("./DATA/RASTER/TopoVars", pattern=".tif$", full.names = TRUE)[c(2,5)]
soilVarsPaths <- list.files("./DATA/RASTER/Soil", pattern=".tif$", full.names = TRUE)[c(1,8)]

# Load raster stack with climatic and topographic data
rstStack <- stack(c(climVarsPaths,topoVarsPaths,soilVarsPaths))
names(rstStack) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOAWC","SOIPH")

# Load shapefile
spData <- readOGR("./DATA/VECTOR/Bryophyte_dataset","And_Gri_Rac_PI_all_2")

# Extract env data to points and remove NA's
# Scale the entire dataset across species
#
spDataVars <- data.frame(spCode = as.character(spData$Cod_esp), 
                         scale(raster::extract(y = spData, x = rstStack)),
                         #raster::extract(y = spData, x = rstStack),
                         stringsAsFactors = FALSE) %>% na.omit

## ----------------------------------------------------------------------------------------------- ##
## Perform the hypervolume analysis per species -----
## ----------------------------------------------------------------------------------------------- ##


# Init data holders for hv volume
hv_gauss_vols <- vector(mode="numeric", length = length(spCodesAll))
hv_svm_vols <- vector(mode="numeric", length = length(spCodesAll))
hv_box_vols <- vector(mode="numeric", length = length(spCodesAll))

# Gather hv objects
hvObj_BySpecies <- list()

# Species code
spCodesAll <- unique(spData$Cod_esp)

if(file.exists("hv_spCodesCheck.csv")){
  
  spCodesCheck <- read.csv("hv_spCodesCheck.csv", stringsAsFactors = FALSE)
  
}else{
  # Make verification table by species
  spCodesCheck <- data.frame(SPECIES = spCodesAll, 
                             #STATUS = c(rep("DONE",27),rep("NOT DONE",12)),
                             STATUS = c(rep("NOT DONE",length(spCodesAll))),
                             stringsAsFactors = FALSE)
  # Write table if it does not exists
  write.csv(spCodesCheck, "hv_spCodesCheck.csv", row.names = FALSE)
}


## Iterate by species

for(sp in spCodesAll){
  
  #i <- i+1
  i <- (1:length(spCodesAll))[spCodesAll==sp]
  
  cat("[",i," /",length(spCodesAll),"] -> Processing species:",sp,"........\n")
  
  if(spCodesCheck[i,"STATUS"] == "DONE"){
    cat("Skipping completed species.\n\n\n")
    next
  }
  
  # Filter / Extract data for rach species
  spDataBySpecies <- spDataVars %>% 
                      filter(spCode == sp) %>% 
                      select(-spCode) #%>% 
                      #scale # not used at species level
  
  hv_gauss <- NULL
  hv_svm <- NULL
  hv_box <- NULL
  
  # Make the hypervolume using different available methods
  hv_gauss <- try(hypervolume_gaussian(spDataBySpecies, name=paste("hv_gauss",sp,sep="_"), verbose = FALSE))
  hv_svm   <- try(hypervolume_svm(spDataBySpecies, name=paste("hv_svm",sp,sep="_"), verbose = FALSE))
  hv_box   <- try(hypervolume_box(spDataBySpecies, name=paste("hv_box",sp,sep="_"), verbose = FALSE))
  
  # Extract volume data for the hv object
  hv_gauss_vols[i] <- ifelse(is.null(hv_gauss) | inherits(hv_gauss,"try-error"), NA, hv_gauss@Volume)
  hv_svm_vols[i] <- ifelse(is.null(hv_svm) | inherits(hv_svm,"try-error"), NA, hv_svm@Volume)
  hv_box_vols[i] <- ifelse(is.null(hv_box) | inherits(hv_box,"try-error"), NA, hv_box@Volume)
  
  hvObj_BySpecies[[name=paste("hv_gauss",sp,sep="_")]] <- hv_gauss
  hvObj_BySpecies[[name=paste("hv_svm",sp,sep="_")]] <- hv_svm
  hvObj_BySpecies[[name=paste("hv_box",sp,sep="_")]] <- hv_box
  
  # Update status
  spCodesCheck[i,"STATUS"] <- "DONE"
  
  cat("\nHv volumes in Log10\n")
  cat("\nHv Gauss =",round(log10(hv_gauss_vols[i]),2),
      "\nHv SVM =",round(log10(hv_svm_vols[i]),2),
      "\nHv Box =",round(log10(hv_box_vols[i]),2),"\n")
  
  write.csv(spCodesCheck, "hv_spCodesCheck.csv", row.names = FALSE)
  
  cat("\nFinished processing species! ------------------------- \n\n\n")
  
}

save.image("HyperVolumeBySpecies-v2.RData")



## ----------------------------------------------------------------------------------------------- ##
## Make dotplot of volumes per species -----
## ----------------------------------------------------------------------------------------------- ##


hv_sp <- data.frame(spCodes=spCodesAll,
                    hv_gauss = hv_gauss_vols,
                    hv_svm = hv_svm_vols,
                    hv_box = hv_box_vols)

hv_sp_log10 <- data.frame(spCodes=spCodesAll,
                    hv_gauss_log = log10(hv_gauss_vols),
                    hv_svm_log = log10(hv_svm_vols),
                    hv_box_log = log10(hv_box_vols))

write.csv(hv_sp,"./RESULTS_2_SHARE/hvolume_bySpecies-v2.csv",row.names = FALSE)
write.csv(hv_sp_log10,"./RESULTS_2_SHARE/log10_hvolume_bySpecies-v2.csv",row.names = FALSE)

# Plot ----------------------------------------------------------------------------- 

hv_sp_ord <- hv_sp %>% arrange(hv_gauss)
hv_sp_ord[,"spCodes"] <- factor(as.character(hv_sp_ord[,"spCodes"]),hv_sp_ord[,"spCodes"],
                                hv_sp_ord[,"spCodes"])


g <- ggplot(hv_sp_ord, aes(x=spCodes, y=hv_gauss)) + 
  geom_point(col="tomato2", size=3) +   # Draw points
  geom_segment(aes(x=spCodes, 
                   xend=spCodes, 
                   y=min(hv_gauss), 
                   yend=max(hv_gauss)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Dot plot of hypervolume sizes", 
       subtitle="Species vs. size of hypervolume") +  
  xlab("Species codes") +
  ylab("Hypervolume size (Hv Gauss)\n← Specialist | Generalist →") +
  coord_flip()

plot(g)

ggsave("./OUT/HypervolumeSizesPerSpecies-HvGauss-v2.png", g, width = 5.5, height = 6)


# Plot  ----------------------------------------------------------------------------- 

hv_sp_ord <- hv_sp %>% arrange(hv_svm)
hv_sp_ord[,"spCodes"] <- factor(as.character(hv_sp_ord[,"spCodes"]),hv_sp_ord[,"spCodes"],
                                hv_sp_ord[,"spCodes"])


g <- ggplot(hv_sp_ord, aes(x=spCodes, y=hv_svm)) + 
  geom_point(col="tomato2", size=3) +   # Draw points
  geom_segment(aes(x=spCodes, 
                   xend=spCodes, 
                   y=min(hv_svm), 
                   yend=max(hv_svm)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Dot plot of hypervolume sizes", 
       subtitle="Species vs. size of hypervolume") +  
  xlab("Species codes") +
  ylab("Hypervolume size (Hv SVM)\n← Specialist | Generalist →") +
  coord_flip()

plot(g)

ggsave("./OUT/HypervolumeSizesPerSpecies-HvSVM-v2.png", g, width = 5.5, height = 6)


# Plot  ----------------------------------------------------------------------------- 

hv_sp_ord <- hv_sp %>% arrange(hv_box)
hv_sp_ord[,"spCodes"] <- factor(as.character(hv_sp_ord[,"spCodes"]),hv_sp_ord[,"spCodes"],
                                hv_sp_ord[,"spCodes"])

g <- ggplot(hv_sp_ord, aes(x=spCodes, y=hv_box)) + 
  geom_point(col="tomato2", size=3) +   # Draw points
  geom_segment(aes(x=spCodes, 
                   xend=spCodes, 
                   y=min(hv_box), 
                   yend=max(hv_box)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Dot plot of hypervolume sizes", 
       subtitle="Species vs. size of hypervolume") +  
  xlab("Species codes") +
  ylab("Hypervolume size (Hv Box)\n← Specialist | Generalist →") +
  coord_flip()

plot(g)

ggsave("./OUT/HypervolumeSizesPerSpecies-HvBox-v2.png", g, width = 5.5, height = 6)





