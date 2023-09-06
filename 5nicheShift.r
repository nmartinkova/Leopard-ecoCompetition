##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Calculates niche shift, niche overlap and niche          #
#   similarity                                               #
#                                                            #
##############################################################



for(i in c("terra")){
  if(!require(i, character.only = TRUE)){
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}


if(!dir.exists("Results")) dir.create("Results")

# whether to recalculate the pairwise correlations
pocitat.korelaciu = FALSE

# colors for legend
res = 60
svetle = 31
farby = c(hcl.colors(svetle, "Viridis"),
  rev(hcl.colors(res-svetle, "Inferno")))

# species that overlap with leopard range
prekryv = read.table("Results/historickePrekryvyALL.txt", header = T, sep = "\t")
druhy = rownames(prekryv)[prekryv$velkost.prekryvu.m2 > 0]


# predicted SDM rasters
selmy = rast(paste0("Data/SDMpredikcie/", druhy[1], ".tiff"))

names(selmy) = druhy[1]

for(druh in druhy[2:length(druhy)]){
	druh.pred = rast(paste0("Data/SDMpredikcie/", druh, ".tiff"))
	names(druh.pred) = druh
	add(selmy) <- druh.pred
}

leopard = rast("Results/Panthera_pardus.tiff")

# leopard niche shift given the biotic variables
selmy$carnivore_impact <- leopard - selmy$Panthera_pardus
writeRaster(selmy$carnivore_impact, filename = "Results/leopard-carnivoreImpact.tiff")

pdf("Results/leopard-carnivoreImpact.pdf", width = 15, height = 8)
  layout(matrix(c(1,2,3,3,3,3), nrow = 2))
  plot(selmy$Panthera_pardus, main = "Climatic & landscape")
  plot(leopard, main = "Climatic, landscape & biotic")
  plot(selmy$carnivore_impact, 
   col = farby, 
   pax = list(retro = TRUE), las = 1,
   main = "Niche shift")
dev.off()



# calculate pairwise correlations
if(pocitat.korelaciu){
  selmy = selmy[[names(selmy) != "Panthera_pardus"]]
  korelacia.pearson = layerCor(selmy, \(x,y) cor(x,y,"pearson", use = "complete.obs"), na.rm= T)
  korelacia.spearman = layerCor(selmy, \(x,y) cor(x,y,"spearman", use = "complete.obs"), na.rm= T)
  save(korelacia.pearson, file = "Results/pairwise-corPearson.RData")
  save(korelacia.spearman, file = "Results/pairwise-corSpearman.RData")

}


# (un)comment which data to consider
# leopard-carnivoreImpact = niche shift
# leopard-interaction = likelihood model

# carnivoreImpact = rast("Results/leopard-carnivoreImpact.tiff")
 carnivoreImpact = rast("Results/leopard-interaction.tiff")

test.pearson = test.spearman = list()
selmy = selmy[[!(names(selmy) %in% c("Panthera_pardus", "carnivore_impact"))]]
# correlations between leopards and all other taxa
for(i in 1:dim(selmy)[3]){
 test.pearson[[i]] = cor.test(unlist(as.data.frame(carnivoreImpact)), 
   unlist(as.data.frame(selmy[[i]])), method = "pearson")
 test.spearman[[i]] = cor.test(unlist(as.data.frame(carnivoreImpact)), 
   unlist(as.data.frame(selmy[[i]])), method = "spearman")
}

names(test.pearson) <- names(test.spearman) <- names(selmy)
# save(test.pearson, file = "Results/carnivoreImpact-pearson.RData")
# save(test.spearman, file = "Results/carnivoreImpact-spearman.RData")
 save(test.pearson, file = "Results/interaction-pearson.RData")
 save(test.spearman, file = "Results/interaction-spearman.RData")



# niche similarity I, niche overlap D

leopard.selmy = rast("Results/Panthera_pardus.tiff")

# standardize rasters
leopard.klima = leopard / unlist(global(leopard, fun = "sum", na.rm = T))
leopard.selmy = leopard.selmy / unlist(global(leopard.selmy, fun = "sum", na.rm = T))

# niche shift in leopard
D = 1 - unlist(global(abs(leopard.klima - leopard.selmy), fun = "sum", na.rm = T)) / 2
I = 1 - unlist(global((sqrt(leopard.klima) - sqrt(leopard.selmy))^2, fun = "sum", na.rm = T)) / 2

cat("Leopard climate and landscape vs climate, landscape and biotic predictors\n\nNiche overlap statistic:\nSchoener's D = ", D, 
  "\n\nNiche similarity:\nI = ", I, file = "Results/leopard-nicheSimilarity.txt")
  
  
# pairwise niche overlap and nich similarity from climatic and landscape data only
selmy = sort(c("Panthera_pardus", selmy))

niche = as.data.frame(t(combn(selmy, m = 2)))
niche$D = niche$I = NA

for(i in 1:nrow(niche)){
  druh1 = rast(paste0("Data/SDMpredikcie/", niche[i,1], ".tiff"))
  druh1 = druh1 / unlist(global(druh1, fun = "sum", na.rm = T))
  druh2 = rast(paste0("Data/SDMpredikcie/", niche[i,2], ".tiff"))
  druh2 = druh2 / unlist(global(druh2, fun = "sum", na.rm = T))
  niche$D[i] = 1 - unlist(global(abs(druh1 - druh2), fun = "sum", na.rm = T)) / 2
  niche$I[i] = 1 - unlist(global((sqrt(druh1) - sqrt(druh2))^2, fun = "sum", na.rm = T)) / 2
}

D = I = matrix(NA, ncol = length(selmy), nrow = length(selmy), dimnames = list(selmy,selmy))
D[lower.tri(D)] = niche$D
I[lower.tri(I)] = niche$I

write.table(D, file = "Results/nicheOverlapD.txt", sep = "\t")
write.table(I, file = "Results/nicheSimilarityI.txt", sep = "\t")

