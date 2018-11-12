

dirs <- list.dirs("/home/cibio/myfiles/MountainBryophytesSDM/OUT/MODS", recursive = FALSE)

#i=1

for(i in 1:length(dirs)){
  
  cmd <- paste("mv",dirs[i],"/mnt/dados/datasets/jg/OUT/MODS")
  print(cmd)
  system(cmd)
}



fl <- list.files("/home/cibio/myfiles/MountainBryophytesSDM/OUT/MODS", 
                 pattern=".RData$",full.names = TRUE)

for(i in 1:length(fl)){
  
  cmd <- paste("mv",fl[i],"/mnt/dados/datasets/jg/OUT/MODS")
  print(cmd)
  system(cmd)
}



list.dirs("/mnt/dados/datasets/jg/OUT/MODS", recursive = FALSE)

fl <- list.files("/mnt/dados/datasets/jg/OUT/MODS", pattern=".RData$",full.names = TRUE)
print(fl)

modCount <- as.data.frame(matrix(nrow=length(fl), ncol=2))

for(i in 1:length(fl)){

  load(fl[i])
  
  modCount[i,1] <- myBiomodData@sp.name
  modCount[i,2] <- length(get_built_models(myBiomodModelOut))
  
  fl <- list.files("/mnt/dados/datasets/jg/OUT/MODS", pattern=".RData$",full.names = TRUE)
  print(fl[i])
}

write.csv(modCount,file = "./OUT/modCount-v1.csv")
