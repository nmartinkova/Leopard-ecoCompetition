##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Calculates likelihood of observing leopards given their  #
#   interactions with other carnivores                       #
#                                                            #
##############################################################

for(i in c("terra")){
  if(!require(i, character.only = TRUE)){
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}


# species that overlap with leopard range
# extinction rates (E_i, E_l_i)
prekryv = read.table("Results/historickePrekryvyALL.txt", header = TRUE, sep = "\t")
selmy = rownames(prekryv)[prekryv$velkost.prekryvu.m2 > 0]


# colors for legend for the likelihood
res = 60
svetle = 48
farby.interakcie = c(colorRampPalette(c(hcl.colors(res-svetle, "Inferno")[res-svetle], hcl.colors(res-svetle, "Inferno")[res-svetle-1]))(svetle),
  rev(hcl.colors(res-svetle, "Inferno")[1:(res-svetle-1)]))



# calculate likelihood of observing leopards given other carnivores (eq. 1)
leopard = rast("Data/SDMpredikcie/Panthera_pardus.tiff")

osud = log10(rast("Data/land_mask.tiff"))

for(druh in selmy){
	cat(druh, "\n")
  iny.druh = rast(paste0("Data/SDMpredikcie/", druh, ".tiff"))
  # raster for probability of meeting (a_i)
  stretnutie = leopard * iny.druh
  mrtvy.leopard = stretnutie * (prekryv[druh, "leopard.vyhynul"])
  zivy.leopard = stretnutie * ((1 - prekryv[druh, "druh.vyhynul"]))
  osud.druh = zivy.leopard - mrtvy.leopard
  osud.druh = ifel(osud.druh <= 0, 1e-10, osud.druh)
  osud = osud + log10(osud.druh)
}

writeRaster(osud, "Results/leopard-interaction.tiff", overwrite = TRUE)


pdf("Results/leopard-interaction.pdf", height = 5, width = 6)
plot(osud, 
 col = farby,
 pax = list(retro = TRUE), las = 1)
dev.off()



