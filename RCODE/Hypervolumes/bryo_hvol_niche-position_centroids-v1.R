
## --------------------------------------------------------- ##
##
## Extract centroids and calculate niche positions
##
## João Gonçalves & Helena Hespanhol
## CIBIO/InBIO, FCUP
## Porto, 10/2018
##
## --------------------------------------------------------- ##

library(dplyr)
library(hypervolume)

# Get overall centroid for all the species records
pointCentroids <- apply(spDataVars %>% select(-1) %>% na.omit, 2, mean)

# Get centroids for the niche hypervolume of each species
for(sp in spCodesAll){
  pointCentroids <- rbind(pointCentroids,
  hypervolume::get_centroid(hvObj_BySpecies[[name=paste("hv_svm",sp,sep="_")]]))
}

rownames(pointCentroids) <- c("AvgCenter",as.character(spCodesAll))

dd <- dist(pointCentroids) %>% as.matrix

dd[,1] %>% sort


hv_DF <- hv_DF %>% mutate(marginDistance = dd[-1,1])


cor.test(hv_DF$hv_svm_log, hv_DF$marginDistance)
cor.test(hv_DF$hv_svm_log[-34], hv_DF$marginDistance[-34]) # Remove RACFAS??
cor.test(hv_DF$hv_svm_log, hv_DF$marginDistance, method="spearman")

g <- ggplot(hv_DF %>% filter(spCode!="RACFAS"),aes(y = marginDistance, x = hv_svm_log)) +
  geom_point(color="tomato", size=2.5) + 
  geom_smooth(method="lm") 
  ylab("Niche position (centroid-to-centroid distance)") + 
  xlab("Niche breadth (Log10 hypervolume size)")

plot(g)

ggsave(filename = "./RESULTS_2_SHARE/Log10Hvolume_vs_MarginalDistance-v1.png", plot=g)

write.csv(hv_DF, "./RESULTS_2_SHARE/hvolumes_marginDist-v1.csv", row.names=FALSE)



