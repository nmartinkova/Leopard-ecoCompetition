##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Prepares land mask, resamples SEDAC rasters to match     #
#   climatic rasters geometry                                #
#                                                            #
##############################################################


for(i in c("terra")){
  if(!require(i, character.only = TRUE)){
    install.packages(i, dependencies = TRUE)
    library(i, character.only = T)
  }
}


AfroEurasia <- vect("POLYGON((34.9 -36.9,  18.4 -40.6,  12.2 -33.5,   7.6 -16.9 ,  3.9  -0.8, -12.8   2.4, -23.2   9.2, -29.4  17.8, -26.2  28.1, -18.2  36.0, -12.6  38.7, -14.5  51.1, -12.1  59.6,  -9.4  65.4,   2.1  71.2 , 24.5  75.0,  41.6  74.7 , 55.3  77.6,  70.7  79.3 , 84.2  83.4 , 97.3  83.2 , 122.2  77.9, 130.2  77.6, 142.0  79.2, 154.6  78.4, 168.2  75.5, 175.8  75.3, 182.1  75.2, 181.9  58.8, 178.6  50.7, 172.6  48.8, 162.0  45.5, 150.6  39.6, 142.6  30.9, 135.2  20.2, 129.6  14.6, 130.8   4.6, 130.1   0.9, 128.9  -1.9, 126.7  -5.3, 129.2  -7.7, 128.4 -10.5, 124.3 -12.1, 115.1 -14.5,  97.9  -8.1,  88.7   3.8,  80.4   2.1,  76.2  -2.9 , 67.8 -17.1  ,57.4 -26.4,  44.7 -32.8,  34.9 -36.9))")
crs(AfroEurasia) = "epsg:4326"

# world mask
land_mask <- get_land_mask(time_ce = 1985, dataset = "WorldClim_2.1_5m")
# Afro-Eurasia mask
land_mask <- crop(land_mask, AfroEurasia)
land_mask <- mask(land_mask, AfroEurasia)


# resample SEDAC layers
vrstvy = dir("Data/SEDAC")
vrstvy = vrstvy[!grepl("env", vrstvy)]
for(i in 1:length(vrstvy)){
  if(i == 1){
    vrstva = resample(rast(paste0("Data/SEDAC/", vrstvy[i])), env, method = "bilinear")
  } else {
    add(vrstva) <- resample(rast(paste0("Data/SEDAC/", vrstvy[i])), env, method = "bilinear")
  }
  setMinMax(vrstva[[i]], force = TRUE)
  vrstva[[i]] = crop(vrstva[[i]], AfroEurasia)
  vrstva[[i]] = mask(vrstva[[i]], land_mask)
  if(i == 3) names(vrstva[[i]]) = "human_modification"
}

# find lakes as NAs is all three SEDAC layers
sucet = mosaic(vrstva[[1]], vrstva[[2]], vrstva[[3]], fun = "sum")
maska = countNA(sucet)
maska = subst(subst(maska, 1, NA), 0, 1)
names(maska) = "land_mask"
# save land mask without lakes
writeRaster(maska, "Data/land_mask.tiff", overwrite = TRUE)

# replace missing values in SEDAC with respective minima
for(i in 1:length(vrstvy)){
  vrstva[[i]] = subst(vrstva[[i]], NA, minmax(vrstva[[i]])[1])
}

# mask SEDAC layers to land
vrstva = mask(vrstva, maska)


# save SEDAC resampled and cropped layers
writeRaster(vrstva[[1]], filename = "Data/SEDAC/env-anthropogenic.tiff", overwrite = TRUE)
writeRaster(vrstva[[2]], filename = "Data/SEDAC/env-cropland.tiff", overwrite = TRUE)
writeRaster(vrstva[[3]], filename = "Data/SEDAC/env-human_modification.tiff", overwrite = TRUE)
