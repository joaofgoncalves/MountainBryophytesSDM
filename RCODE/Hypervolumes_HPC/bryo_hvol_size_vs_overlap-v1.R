
library(dplyr)
library(ggplot2)




## ----------------------------------------------------------------------------------------------- ##
## Load env and species data -----
## ----------------------------------------------------------------------------------------------- ##

# Data paths to raster variables
climVarsPaths <- list.files("./DATA/RASTER/WorldClim_1km/a2000", pattern=".tif$", full.names = TRUE)
topoVarsPaths <- list.files("./DATA/RASTER/TopoVars", pattern=".tif$", full.names = TRUE)[c(2,3,5)]

# Load raster stack with climatic and topographic data
rstStack <- stack(c(climVarsPaths,topoVarsPaths))
names(rstStack) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","ELEV","TOPRI")

# Load shapefile
spData <- readOGR("./DATA/VECTOR/Bryophyte_dataset","And_Gri_Rac_PI_all_2")

# Extract env data to points and remove NA's
# Scale the entire dataset across species
#
spDataVars <- data.frame(spCode = as.character(spData$Cod_esp), 
                         #scale(raster::extract(y = spData, x = rstStack)),
                         raster::extract(y = spData, x = rstStack),
                         stringsAsFactors = FALSE) %>% 
  na.omit


## ----------------------------------------------------------------------------------------------- ##
## Analyze hypervolume results -----
## ----------------------------------------------------------------------------------------------- ##

# Load pre-computed results
load("./OUT/HyperVolumeBySpecies-v1.RData")
load("./OUT/NicheOvlpDistances_v1.RData")

# Create full distance matrix based on Jaccard similarity of niches
dd <- as.dist(ovlp_jacc, upper=TRUE)

# Calculate the average overlap and get Hvolume
DF <- data.frame(spCode = spCodesAll, 
                 avg_ovlp=apply(as.matrix(dd), 1, mean), 
                 std_ovlp=apply(as.matrix(dd), 1, sd),
                 hv_svm=hv_svm_vols, stringsAsFactors=FALSE)

# Plot the hypervolume size vs the average overlap between species
plot(DF$avg_ovlp, log10(DF$hv_svm))

# Get the average elevation per species
avgElevs <- spDataVars %>% 
  group_by(spCode) %>% 
  summarize(avgElev = mean(ELEV), stdElev = sd(ELEV)) %>% 
  as.data.frame

# Merge data
tb <- merge(DF,avgElevs, by="spCode")

# Plot hv size vs elevation
plot(log10(tb$hv_svm), log10(tb$avgElev))
plot(tb$hv_svm, tb$avgElev)

# Perform a hierarchical cluster analysis based on these two variables
tbb <- scale(tb[,c(3:4)])
rownames(tbb) <- tb$spCode 
ddd <- dist(tbb)
hc <- hclust(ddd)
plot(hc, labels=tb$spCode)
