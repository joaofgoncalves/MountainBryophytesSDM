
library(raster)
library(biomod2)
library(ggplot2)
library(dplyr)


## ---------------------------------------------------------------------------------------------- ##
## Get predicted range size by species for current and future scenarios ----
## ---------------------------------------------------------------------------------------------- ##


fl <- list.files("/mnt/dados/datasets/jg/OUT/MODS", pattern=".RData$", full.names = TRUE)

spNames <- list.dirs("/mnt/dados/datasets/jg/OUT/MODS", full.names = FALSE, recursive = FALSE)

dirs <- list.dirs("/mnt/dados/datasets/jg/OUT/MODS", full.names = TRUE, recursive = FALSE)

projNames <- c("current","he45bi50","he85bi50","mp45bi50","mp85bi50")


for(projName in projNames){
  
  cat("\n\n-> Working with projection:",projName,"......\n\n")
  pb <- txtProgressBar(min=1,max=length(spNames),style=3)
  
  for(i in 1:length(spNames)){
    
    rstProj <- paste(dirs[i],"/proj_",projName,"/GeoTIFF/proj_",projName,"_",spNames[i],"_ensemble_TSSbin.tif", sep="")
    
    rstStackProj <- stack(rstProj)
    
    # Get raster data
    rstDF <- na.omit(values(rstStackProj))
    # Calculate the amount of suitable area (in terms of pixels)
    dfTmp <- data.frame(spNames[i], t(apply(rstDF[,1:3],2,sum)), nrow(rstDF))
    colnames(dfTmp) <- c("spNames", "suitArea_TSS", "suitArea_ROC", "suitArea_KAPPA", "totalArea_npix")
    
    if(i==1){
      DF_ProjArea <- dfTmp
    }else{
      DF_ProjArea <- rbind(DF_ProjArea, dfTmp)
    }
    
    setTxtProgressBar(pb,i)
  }
  
  write.csv(DF_ProjArea, file = paste("ProjectedAreaBySpecies_proj_",projName,"-v1.csv",sep=""), row.names = FALSE)
  rm(DF_ProjArea)
  
}


## ---------------------------------------------------------------------------------------------- ##
## Get binary predictions (TSS) for current period ----
##
## (this data is used to enable the calculation of Jaccard and Sorensen similarity indices as a 
## proxy to the overlap of species)
## ---------------------------------------------------------------------------------------------- ##


projName <- "current"
pb <- txtProgressBar(min=1,max=length(spNames),style=3)

for(i in 1:length(spNames)){
  
  rstProj <- paste(dirs[i],"/proj_",projName,"/GeoTIFF/proj_",projName,"_",spNames[i],"_ensemble_TSSbin.tif", sep="")
  
  tmp <- stack(rstProj)[[1]]
  
  if(i==1){
    rstStackProj_curr <- tmp
  }else{
    
    # A small hack to correct for different extents in the projections...
    if(!compareRaster(tmp,rstStackProj_curr,stopiffalse = FALSE)){
      tmp <- projectRaster(tmp, rstStackProj_curr, method = 'ngb')
    }
    
    rstStackProj_curr <- stack(rstStackProj_curr, tmp)
  }
  setTxtProgressBar(pb,i)
  
}

rstBinDF <- values(rstStackProj_curr) %>% na.omit

colnames(rstBinDF) <- spNames

saveRDS(rstBinDF,file = "./OUT/rstBinDF_TSS_v3.RData")


## ---------------------------------------------------------------------------------------------- ##
## Calculate niche overlap indices based on dismo::nicheOverlap function ----
## ---------------------------------------------------------------------------------------------- ##


library(dismo)


# Push data for the average ensemble and make a compatible RasterStack
projName <- "current"
pb <- txtProgressBar(min=1,max=length(spNames),style=3)

for(i in 1:length(spNames)){
  
  rstProj <- paste(dirs[i],"/proj_",projName,"/GeoTIFF/proj_",projName,"_",spNames[i],"_ensemble.tif", sep="")
  
  tmp <- stack(rstProj)[[1]]
  
  if(i==1){
    rstStackProj_curr <- tmp
  }else{
    
    # A small hack to correct for different extents in the projections...
    if(!compareRaster(tmp,rstStackProj_curr,stopiffalse = FALSE)){
      tmp <- projectRaster(tmp, rstStackProj_curr, method = 'ngb')
    }
    
    rstStackProj_curr <- stack(rstStackProj_curr, tmp)
  }
  setTxtProgressBar(pb,i)
  
}

# Calculate the indices based on dismo::nicheOverlap function ------- #

# Matrix size equal to the number of species
len <- nlayers(rstStackProj_curr)

# Start the matrix holding calculations
ltrimat <- matrix(1:len^2,len,len) %>% lower.tri
ovlp_indI <- matrix(NA,len,len,dimnames = list(spNames,spNames))
ovlp_indD <- matrix(NA,len,len,dimnames = list(spNames,spNames))

tot2run <- sum(ltrimat)
pb <- txtProgressBar(1,tot2run,style = 3)

# Iterate through the distance bu only in the lower triangular part
k <- 0
for(i in 1:len){
  for(j in 1:len){
    
    if(!ltrimat[i,j]){ # skip if not the lower-tri matrix
      next 
    }
    
    k <- k+1
    #sp1 <- as.character(spNames[i])
    #sp2 <- as.character(spNames[j])
    r1 <- rstStackProj_curr[[i]]
    r2 <- rstStackProj_curr[[j]]
    
    ovlpInd_I <- nicheOverlap(r1, r2, stat = "I")
    ovlpInd_D <- nicheOverlap(r1, r2, stat = "D")
    
    ovlp_indI[i,j] <- ovlpInd_I
    ovlp_indD[i,j] <- ovlpInd_D
    
    setTxtProgressBar(pb, k)
  }
}

saveRDS(ovlp_indI, "./OUT/dismo_nicheOverlap_indI_v3.rds")
saveRDS(ovlp_indD, "./OUT/dismo_nicheOverlap_indD_v3.rds")


## ---------------------------------------------------------------------------------------------- ##
## Calculate habitat suitability dynamics based on raster combos ----
## ---------------------------------------------------------------------------------------------- ##


allRstStacks <- list()

for(projName in projNames){
  
  cat("\n\n-> Working with projection:",projName,"......\n\n")
  pb <- txtProgressBar(min=1,max=length(spNames),style=3)
  

  for(i in 1:length(spNames)){
    
    rstProj <- paste(dirs[i],"/proj_",projName,"/GeoTIFF/proj_",projName,"_",spNames[i],"_ensemble_TSSbin.tif", sep="")
    
    tmp <- stack(rstProj)[[1]]
    
    if(i==1){
      rstStackProj_curr <- tmp
    }else{
      
      # A small hack to correct for different extents in the projections...
      if(!compareRaster(tmp, rstStackProj_curr, stopiffalse = FALSE)){
        tmp <- projectRaster(tmp, rstStackProj_curr, method = 'ngb')
      }
      
      rstStackProj_curr <- stack(rstStackProj_curr, tmp)
    }
    setTxtProgressBar(pb,i)
    
  }
  
  names(rstStackProj_curr) <- spNames
  
  # Make a complete RasterStack with all projections
  allRstStacks[[projName]] <- rstStackProj_curr
  
  # Write the raster data for all species
  writeRaster(rstStackProj_curr, filename = paste("./OUT/MODS/BIN_PREDS/binPreds_allSps_TSS_proj_",projName,".tif",sep=""))
}


# Start the matrices used to hold results from the overlap analyses
#
dynHab_he45 <- matrix(0,nrow=length(spNames),ncol=4,
                      dimnames=list(spNames,paste("p",1:4,sep="")))
dynHab_he85 <- matrix(0,nrow=length(spNames),ncol=4,
                      dimnames=list(spNames,paste("p",1:4,sep="")))
dynHab_mp45 <- matrix(0,nrow=length(spNames),ncol=4,
                      dimnames=list(spNames,paste("p",1:4,sep="")))
dynHab_mp85 <- matrix(0,nrow=length(spNames),ncol=4,
                      dimnames=list(spNames,paste("p",1:4,sep="")))

pb <- txtProgressBar(min=1,max=length(spNames),style=3)

i<-0
for(spName in spNames){
  
  i<-i+1 # Iterator count
  
  # Current projection for the species
  r1 <- allRstStacks[["current"]][[spName]]
  
  # Future projections for each GCM/RCP combination
  r2 <- allRstStacks[["he45bi50"]][[spName]]
  r3 <- allRstStacks[["he85bi50"]][[spName]]
  r4 <- allRstStacks[["mp45bi50"]][[spName]]
  r5 <- allRstStacks[["mp85bi50"]][[spName]]
  
  # Compare current vs he45bi50 projections ---------- #
  #
  rr1 <- r1
  rr1[(r1==0) & (r2==0)] <- 1 # Stable unsuitable
  rr1[(r1==0) & (r2==1)] <- 2 # Gain
  rr1[(r1==1) & (r2==1)] <- 3 # Stable suitable
  rr1[(r1==1) & (r2==0)] <- 4 # Loss 
  
  # Compare current vs he85bi50 projections ---------- #
  #
  rr2 <- r1
  rr2[(r1==0) & (r3==0)] <- 1 # Stable unsuitable
  rr2[(r1==0) & (r3==1)] <- 2 # Gain
  rr2[(r1==1) & (r3==1)] <- 3 # Stable suitable
  rr2[(r1==1) & (r3==0)] <- 4 # Loss 
  
  # Compare current vs mp45bi50 projections ---------- #
  #
  rr3 <- r1
  rr3[(r1==0) & (r4==0)] <- 1 # Stable unsuitable
  rr3[(r1==0) & (r4==1)] <- 2 # Gain
  rr3[(r1==1) & (r4==1)] <- 3 # Stable suitable
  rr3[(r1==1) & (r4==0)] <- 4 # Loss 
  
  # Compare current vs mp85bi50 projections ---------- #
  #
  rr4 <- r1
  rr4[(r1==0) & (r5==0)] <- 1 # Stable unsuitable
  rr4[(r1==0) & (r5==1)] <- 2 # Gain
  rr4[(r1==1) & (r5==1)] <- 3 # Stable suitable
  rr4[(r1==1) & (r5==0)] <- 4 # Loss 
  
  # Calculate cross-tabulations ---------------------- #
  # and save values to matrices
  #
  rr1_tb <- table(values(rr1)) %>% as.matrix %>% t
  colnames(rr1_tb) <- paste("p",colnames(rr1_tb),sep="")
  dynHab_he45[spName ,colnames(rr1_tb)] <- rr1_tb[1,]
    
  rr2_tb <- table(values(rr2)) %>% as.matrix %>% t
  colnames(rr2_tb) <- paste("p",colnames(rr2_tb),sep="")
  dynHab_he85[spName ,colnames(rr2_tb)] <- rr2_tb[1,]
  
  rr3_tb <- table(values(rr3)) %>% as.matrix %>% t
  colnames(rr3_tb) <- paste("p",colnames(rr3_tb),sep="")
  dynHab_mp45[spName ,colnames(rr3_tb)] <- rr3_tb[1,]
  
  rr4_tb <- table(values(rr4)) %>% as.matrix %>% t
  colnames(rr4_tb) <- paste("p",colnames(rr4_tb),sep="")
  dynHab_mp85[spName ,colnames(rr4_tb)] <- rr4_tb[1,]
  
  setTxtProgressBar(pb,i)
}

## Calculate net habitat gains in % ----------------- #
##

d1 <- apply(dynHab_he45[,3:4],1,sum)
dynHab_he45_rp <- cbind((dynHab_he45[,2:4] / d1) * 100 , net = (((dynHab_he45[,2] + dynHab_he45[,3]) - dynHab_he45[,4]) / d1) *100 )
colnames(dynHab_he45_rp) <- paste(c("pGain","pStable","pLoss","netHab"),"_he45",sep="")

d1 <- apply(dynHab_he85[,3:4],1,sum)
dynHab_he85_rp <- cbind((dynHab_he85[,2:4] / d1) * 100 , net = (((dynHab_he85[,2] + dynHab_he85[,3]) - dynHab_he85[,4]) / d1) *100 )
colnames(dynHab_he85_rp) <- paste(c("pGain","pStable","pLoss","netHab"),"_he85",sep="")

d1 <- apply(dynHab_mp45[,3:4],1,sum)
dynHab_mp45_rp <- cbind((dynHab_mp45[,2:4] / d1) * 100 , net = (((dynHab_mp45[,2] + dynHab_mp45[,3]) - dynHab_mp45[,4]) / d1) *100 )
colnames(dynHab_mp45_rp) <- paste(c("pGain","pStable","pLoss","netHab"),"_mp45",sep="")

d1 <- apply(dynHab_mp85[,3:4],1,sum)
dynHab_mp85_rp <- cbind((dynHab_mp85[,2:4] / d1) * 100 , net = (((dynHab_mp85[,2] + dynHab_mp85[,3]) - dynHab_mp85[,4]) / d1) *100 )
colnames(dynHab_mp85_rp) <- paste(c("pGain","pStable","pLoss","netHab"),"_mp85",sep="")

# Join all matrices
relHabShifts <- data.frame(spName = rownames(dynHab_he45),
                      dynHab_he45_rp,
                      dynHab_he85_rp,
                      dynHab_mp45_rp,
                      dynHab_mp85_rp,
                      stringsAsFactors = FALSE)

write.csv(relHabShifts, "./OUT/relHabitatChange_BySpScenario_v3.csv", row.names = FALSE)


