
library(dplyr)
library(ggplot2)
library(raster)
library(rgdal)
library(hypervolume)

# Load pre-computed results
load("./OUT/HyperVolumeBySpecies-v2.RData")
load("./OUT/NicheOvlpDistances-NewVars_v2.RData")

hvols_DF <- readxl::read_excel("RESULTS_2_SHARE/hvolumes_marginDist-v1.xlsx") %>% as.data.frame

## ----------------------------------------------------------------------------------------------- ##
## Load env and species data -----
## ----------------------------------------------------------------------------------------------- ##

# Data paths to raster variables
climVarsPaths <- list.files("./DATA/RASTER/WorldClim_1km/a2000", pattern=".tif$", full.names = TRUE)
topoVarsPaths <- list.files("./DATA/RASTER/TopoVars", pattern=".tif$", full.names = TRUE)[c(2,3,5)]
soilVarsPaths <- list.files("./DATA/RASTER/Soil", pattern=".tif$", full.names = TRUE)[c(1,8)]


# Load raster stack with climatic and topographic data
rstStack <- stack(c(climVarsPaths,topoVarsPaths,soilVarsPaths))
names(rstStack) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","ELEV","TOPRI","SOAWC","SOIPH")

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


spDataVarsElev <- spDataVars

## ----------------------------------------------------------------------------------------------- ##
## Analyze hypervolume results ----
## Niche position/breadth vs. niche overlap
## ----------------------------------------------------------------------------------------------- ##

# Create full distance matrix based on Jaccard similarity of niches
dd <- as.dist(ovlp_jacc, upper=TRUE)
ddmat <- as.matrix(dd)
diag(ddmat) <- NA

spCodesAll <- attr(dd,"Labels")

# Calculate the average overlap and get Hvolume
DF <- data.frame(hvols_DF,
                 avg_ovlp=apply(ddmat, 1, mean, na.rm=TRUE), 
                 std_ovlp=apply(ddmat, 1, sd, na.rm=TRUE),
                 stringsAsFactors=FALSE)

# Plot the hypervolume size vs the average overlap between species
plot(DF$hv_svm_log, DF$avg_ovlp)
cor.test(DF$hv_svm_log, DF$avg_ovlp)
cor.test(DF$hv_svm_log, DF$avg_ovlp, method="spearman")

plot(DF$hv_svm_log, DF$marginDistance)
cor.test(DF$hv_svm_log, DF$marginDistance)
cor.test(DF$hv_svm_log, DF$marginDistance, method="spearman")

g <- ggplot(DF %>% filter(spCode!="RACFAS"),
            aes(x    = hv_svm_log, 
                y    = avg_ovlp,
                ymin = avg_ovlp - 0.5*std_ovlp, 
                ymax = avg_ovlp + 0.5*std_ovlp)) + 
      geom_errorbar() + 
      geom_point(size=2.5, color="tomato") + 
     geom_smooth(method="gam", formula = y ~ s(x,k=2)) +
     xlab("Niche breadth (Log 10 hypervolume size)") + 
     ylab("Average niche overlap (Jaccard similarity)")
     #geom_quantile()

ggsave(filename = "./RESULTS_2_SHARE/NicheBreadth-vs-AvgNicheOverlap-v1.png", plot=g)

plot(g)

g <- ggplot(DF %>% filter(spCode!="RACFAS"),
            aes(x    = marginDistance, 
                y    = avg_ovlp,
                ymin = avg_ovlp - 0.5*std_ovlp, 
                ymax = avg_ovlp + 0.5*std_ovlp)) + 
  geom_errorbar() + 
  geom_point(size=2.5, color="tomato") + 
  geom_smooth(method="gam", formula = y ~ s(x,k=2)) +
  xlab("Niche position (centroid-to-centroid distance)") + 
  ylab("Average niche overlap (Jaccard similarity)")
#geom_quantile()

plot(g)

ggsave(filename = "./RESULTS_2_SHARE/NichePosition-vs-AvgNicheOverlap-v1.png", plot=g)


g1 <- mgcv::gam(avg_ovlp~s(hv_svm_log),data=DF)
summary(g1)
AIC(g1)

g2 <- mgcv::gam(avg_ovlp~s(marginDistance),data=DF)
summary(g2)
AIC(g2)


## ----------------------------------------------------------------------------- ##

# Get the average elevation per species
avgElevs <- spDataVars %>%
  group_by(spCode) %>%
  summarize(avgElev = mean(ELEV), stdElev = sd(ELEV), 
            nRecords=length(ELEV)) %>%
  as.data.frame

# Merge data
tb <- merge(DF,avgElevs, by="spCode")

#cor.test(tb$hv_svm, tb$nRecords, method="spearman")
cor.test(tb$hv_svm, tb$nRecords, method="pearson")
plot(hv_svm ~nRecords, data=tb)


# Plot hv size vs elevation
plot(log10(tb$hv_svm), log10(tb$avgElev))
plot(tb$hv_svm, tb$avgElev)

# Perform a hierarchical cluster analysis based on these two variables
tbb <- scale(tb[,c(3:4)])
rownames(tbb) <- tb$spCode 
ddd <- dist(tbb)
hc <- hclust(ddd)
plot(hc, labels=tb$spCode)

## ----------------------------------------------------------------------------- ##

# Get the average elevation per species
hv_DF <- spDataVarsElev %>%
  group_by(spCode) %>%
  summarize(avgElev = mean(ELEV), stdElev = sd(ELEV), 
            nRecords=length(ELEV)) %>%
  as.data.frame %>% 
  merge(hv_sp, by.x="spCode",by.y="spCodes") %>% 
  merge(hv_sp_log10, by.x="spCode",by.y="spCodes") 


cor.test(hv_DF$avgElev, hv_DF$hv_svm_log)
cor.test(hv_DF$stdElev, hv_DF$hv_svm_log)

cor.test(hv_DF$avgElev, hv_DF$marginDistance)
cor.test(hv_DF$stdElev, hv_DF$marginDistance)

## ----------------------------------------------------------------------------- ##

g <- ggplot(hv_DF,aes(y=avgElev,x=hv_svm_log)) +
  geom_errorbar(aes(ymin=avgElev - 0.5*stdElev, ymax=avgElev+0.5*stdElev)) + 
  geom_point(color="tomato", size=2.5) + 
  geom_smooth(method="lm") + 
  ylab("Elevation (m)") + 
  xlab("Niche breadth (Log10 hypervolume size)")
  
plot(g)

ggsave(filename = "./RESULTS_2_SHARE/NicheBreadth-Log10Hvolume_vs_Elevation-v1.png", plot=g)

## ----------------------------------------------------------------------------- ##

g <- ggplot(hv_DF,aes(y=avgElev,x=marginDistance)) +
  geom_errorbar(aes(ymin=avgElev - 0.5*stdElev, ymax=avgElev+0.5*stdElev)) + 
  geom_point(color="tomato", size=2.5) + 
  geom_smooth(method="lm") + 
  ylab("Elevation (m)") + 
  xlab("Niche position (centroid distance)")

plot(g)

ggsave(filename = "./RESULTS_2_SHARE/NichePos-MarginDist_vs_Elevation-v1.png", plot=g)


## ----------------------------------------------------------------------------- ##

g <- ggplot(hv_DF,aes(y=stdElev,x=hv_svm_log)) +
  #geom_errorbar(aes(ymin=avgElev - 0.5*stdElev, ymax=avgElev+0.5*stdElev)) + 
  geom_point(color="tomato", size=2.5) + 
  geom_smooth(method="lm") + 
  ylab("Elevation std-deviation (m)") + 
  xlab("Log10 hypervolume size")

plot(g)

ggsave(filename = "./RESULTS_2_SHARE/Log10Hvolume_vs_StdElevation-v1.png", plot=g)



