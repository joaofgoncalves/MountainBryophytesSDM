

library(raster)
library(dplyr)

r <- stack(list.files("./DATA/RASTER/TopoVars", full.names = TRUE, pattern = ".tif$"))

DF <- r %>% values %>% na.omit

cmSpear <- cor(DF, method = "spearman") %>% round(digits = 2)

View(cmSpear)

# aspect_beer / slope / tpi_3r / tri_3r

r1 <- stack(c(list.files("./DATA/RASTER/TopoVars", full.names = TRUE, pattern = ".tif$"),
            list.files("./DATA/RASTER/Worldclim_1km/a2000", full.names = TRUE, pattern = ".tif$")))


DF1 <- r1 %>% values %>% na.omit

cmSpear1 <- cor(DF1, method = "spearman") %>% round(digits = 2)

View(cmSpear1)


