##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Calculates species distribution models (SDMs) for        #
#   carnivores that overlap with leopard                     #
#                                                            #
##############################################################



for(i in c("pastclim", "terra", "tidysdm", "sf", "EFA.dimensions", "xgboost")){
  if(!require(i, character.only = TRUE)){
    stop("Run script 0.5setup.r first")
  }
}

# debugging controls, TRUE recommended
preskocit.hotove = TRUE
verbose = TRUE

# create folders for results
if(!dir.exists("Results")) dir.create("Results")
if(!dir.exists("Data/SDMpredikcie")) dir.create("Data/SDMpredikcie")
if(verbose & !dir.exists("Results/log")) dir.create("Results/log")



# species that overlap with leopard range
prekryv = read.table("Results/historickePrekryvyALL.txt", header = T, sep = "\t")
selmy = c("Panthera_pardus", rownames(prekryv)[prekryv$velkost.prekryvu.m2 > 0])


if(preskocit.hotove){
	hotove = unique(sub("-.+", "", dir("Results/log")))
	selmy = selmy[!(selmy %in% hotove)]
	cat(selmy, sep = "\n")
}


# Afro-Eurasia for climatic data
AfroEurasia <- vect("POLYGON((34.9 -36.9,  18.4 -40.6,  12.2 -33.5,   7.6 -16.9 ,  3.9  -0.8, -12.8   2.4, -23.2   9.2, -29.4  17.8, -26.2  28.1, -18.2  36.0, -12.6  38.7, -14.5  51.1, -12.1  59.6,  -9.4  65.4,   2.1  71.2 , 24.5  75.0,  41.6  74.7 , 55.3  77.6,  70.7  79.3 , 84.2  83.4 , 97.3  83.2 , 122.2  77.9, 130.2  77.6, 142.0  79.2, 154.6  78.4, 168.2  75.5, 175.8  75.3, 182.1  75.2, 181.9  58.8, 178.6  50.7, 172.6  48.8, 162.0  45.5, 150.6  39.6, 142.6  30.9, 135.2  20.2, 129.6  14.6, 130.8   4.6, 130.1   0.9, 128.9  -1.9, 126.7  -5.3, 129.2  -7.7, 128.4 -10.5, 124.3 -12.1, 115.1 -14.5,  97.9  -8.1,  88.7   3.8,  80.4   2.1,  76.2  -2.9 , 67.8 -17.1  ,57.4 -26.4,  44.7 -32.8,  34.9 -36.9))")
crs(AfroEurasia) = "epsg:4326"


# land mask with lakes
land_mask <- rast("Data/land_mask.tiff")

# variables to consider
premenne = get_vars_for_dataset("WorldClim_2.1_5m")

# climatic variables in AfroEurasia
env = region_slice(time_ce = 1985, bio_variables = premenne,
	data = "WorldClim_2.1_5m", crop = AfroEurasia)


# add SEDAC layers
vrstvy = dir("Data/SEDAC")
vrstvy = vrstvy[grepl("env-", vrstvy)]
for(i in 1:length(vrstvy)){
	add(env) <- rast(paste0("Data/SEDAC/", vrstvy[i]))
}

# unify environmental layers for Afro-Eurasia and variable names
env = mask(env, land_mask)
premenne = names(env)


# SDM from climatic and landscape parameters
for(druh in selmy){
  
  message("Starting ", druh, ", ", Sys.time())
  # species occurrence records projection
  nalezy = read.table(paste0("Data/GBIF/", druh, ".csv"), 
    header = TRUE, sep = "\t", 
    quote = "\"", comment.char = "")
  nalezy = nalezy[,c("decimalLongitude", "decimalLatitude")]
  nalezy = nalezy[complete.cases(nalezy),]
  nalezy = st_as_sf(nalezy, coords = 1:2, crs = 4326)
  
  # filter occurrence records within historical species range
  areal = st_geometry(st_make_valid(read_sf(paste0("Data/IUCN/", druh, "/data_0.shp"))))
  if(!(druh %in% c("Melursus_ursinus", "Panthera_uncia", "Ursus_arctos", "Ursus_thibetanus"))){
    # ! Loop 0 is not valid: Edge 388 crosses edge 390
    # the exempt taxa need to be manually checked, records within IUCN range in August 2023
    nalezy = st_filter(nalezy, areal)
  }
  if(druh == "Ursus_thibetanus"){
    nalezy = st_filter(nalezy, st_polygon(list(matrix(c(40, 1, 146, -6, 160, 62, 46, 66, 40, 1), byrow = TRUE, ncol = 2))))
  }
  
  # records thinning (stochastic process; fix seed for debugging)
  nalezy = thin_by_cell(nalezy, raster = land_mask)
  nalezy = thin_by_dist(nalezy, dist_min = km2m(30))


  # coordinates for background points (stochastic process; fix seed for debugging)
  k = ifelse(nrow(nalezy) < 500, 4, 3)
  pozadie <- sample_pseudoabs(nalezy, n = k * nrow(nalezy), 
     raster = land_mask, method = c("random"))
  st_crs(pozadie) = 4326
  cat("Pseudoabsences chosen, ", nrow(pozadie), "rows for ", nrow(nalezy), " retained occurrence points\n")
  
  # environmental data for background points
  dat = pozadie %>% bind_cols(terra::extract(env, pozadie, ID = FALSE))
  cat("Background done\n")

  # PCA from background environmental variables
  pca = prcomp(as.data.frame(dat)[(nrow(nalezy) + 1):nrow(dat), premenne], scale = TRUE)
  # Velicer's minimum average partial (MAP) test
  # number of PCs to retain in SDM
  npcs = MAP(as.data.frame(dat)[(nrow(nalezy) + 1):nrow(dat), premenne], 
    "spearman", verbose = F)$NfactorsMAP4
  cat("PCA done\n")
  
  
  # predict PC loadings in dat (for work with models) and env (for work with predictions)
  dat = bind_cols(dat, predict(pca, as.data.frame(dat)))
  env.pca = predict(env, pca)
  cat("PCA prediction done\n")
  
  # keep only informative PCs
  dat = dat %>% select(all_of(c(paste0("PC", 1:npcs), "class")))
  cat("Informative PCs selected\n")

  # scheme for SDM models
  env_recipe = recipe(dat, formula = class ~ .)
  env_modely = workflow_set(preproc = list(default = env_recipe),
    models = list(glm = sdm_spec_glm(),
      maxent = sdm_spec_maxent(),
      rf = sdm_spec_rf(),
      gbm = sdm_spec_boost_tree()),
    cross = TRUE) %>% option_add(control = control_ensemble_grid())
  
  # scheme for spatial cross-validation
  env_cv = spatial_block_cv(dat, v = 5)
  cat("CV blocks done. Starting fitting SDMs\n")

  # run SDM
  env_modely = env_modely %>% workflow_map("tune_grid", resamples = env_cv, 
    grid = 5, metrics = sdm_metric_set(), verbose = TRUE)

  # select tuned models for each method based on AUC
  modely_spojene = simple_ensemble() %>% add_member(env_modely, metric = "roc_auc")
  cat("Best models selected\nStarting SDM prediction\n")

  # SDM prediction across Afro-Eurasia
  predikcia = predict_raster(modely_spojene, env.pca, metric_thresh = c("roc_auc", .75))
  cat("Prediction for ", druh, " done\n")
  
  # save SDM prediction for the species
  writeRaster(predikcia, filename = paste0("Data/SDMpredikcie/", druh, ".tiff"), overwrite = TRUE)

  
  # save PCA interpretation, model performance and pdf of SDM prediction
  if(verbose){
    pdf(paste0("Results/log/", druh, "-predikcia.pdf"), height = 5, width = 6)
      plot(predikcia)
      plot(nalezy, add = TRUE, pch = ".")
    dev.off()
    capture.output(modely_spojene$metrics, file = paste0("Results/log/", druh, "-modelPerformance.txt"))
    save(pca, file = paste0("Results/log/", druh, "-PCAmodel.RData"))
    pdf(paste0("Results/log/", druh, "-PCAinterpretacia.pdf"), height = 5, width = 12)
      par(mar = c(9,4,0,0)+.1)
      image(1:(nrow(pca$rotation)+1), 1:(npcs+1), abs(pca$rotation[,1:npcs]), 
        axes = F, xlab = "", ylab = "")
      axis(1, at = 1.5:(nrow(pca$rotation)+.5), labels = rownames(pca$rotation), las = 2)
      axis(2, at = 1.5:(npcs+.5), labels = colnames(pca$rotation)[1:npcs], las = 1)
    dev.off()
  }

}