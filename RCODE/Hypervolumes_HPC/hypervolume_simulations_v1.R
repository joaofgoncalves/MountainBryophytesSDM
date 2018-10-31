
library(geometry)
library(hypervolume)
library(cluster)

pts1 <- matrix(rnorm(1000, mean = 100, sd = 1),  nrow = 200, ncol = 5)
pts2 <- matrix(rnorm(1000, mean = 100, sd = 10), nrow = 200, ncol = 5)

qhull1 <- convhulln(pts1, options = "Tv FA")
qhull2 <- convhulln(pts2, options = "Tv FA")

qhull1$area
qhull1$vol

qhull2$area
qhull2$vol

elhull1 <- ellipsoidhull(pts1)
elhull2 <- ellipsoidhull(pts2)

cluster:::volume.ellipsoid(elhull1)
cluster:::volume.ellipsoid(elhull2)

print(elhull1)
print(elhull2)



sds <- seq(1, 10, by = 0.5)
eh_vols <- vector(mode="numeric", length = length(sds))
ch_vols <- vector(mode="numeric", length = length(sds))
hv_gauss_vols <- vector(mode="numeric", length = length(sds))
hv_svm_vols <- vector(mode="numeric", length = length(sds))
hv_box_vols <- vector(mode="numeric", length = length(sds))

i<-0
for(sdVal in sds){
  
  i <- i+1
  pts1 <- matrix(rnorm(600, mean = 0, sd = sdVal),  nrow = 200, ncol = 3)
  elhull1 <- ellipsoidhull(pts1)
  qhull1 <- convhulln(pts1, options = "Tv FA")
  hv_gauss <- hypervolume_gaussian(pts1, name="test", verbose = FALSE)
  hv_svm <- hypervolume_svm(pts1, name="test", verbose = FALSE)
  hv_box <- hypervolume_box(pts1, name="test", verbose = FALSE)
  
  eh_vols[i] <- cluster:::volume.ellipsoid(elhull1)
  ch_vols[i] <- qhull1$vol
  hv_gauss_vols[i] <- hv_gauss@Volume
  hv_svm_vols[i] <- hv_svm@Volume
  hv_box_vols[i] <- hv_box@Volume
}

all_vols <- data.frame(eh=eh_vols,
           ch=ch_vols,
           gs=hv_gauss_vols,
           sv=hv_svm_vols,
           bx=hv_box_vols)

round(cor(all_vols, method="spearman"), 2)

plot(sds, log10(eh_vols), xlab="Standard-deviation",ylab="Log10 Ellipsoid-hull volume")
plot(sds, log10(ch_vols))

plot(sds, log10(hv_vols), xlab="Standard-deviation",ylab="Log10 Hypervolume volume")
plot(sds, log10(hv_vols))

plot(eh_vols,hv_vols)
plot(ch_vols,hv_vols)

cor.test(eh_vols,hv_vols)
cor.test(ch_vols,hv_vols)

hypervolume_distance(hv_gauss,hv_box, num.points.max = 1E6)

