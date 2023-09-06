##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Plots figures with main results                          #
#                                                            #
##############################################################




for(i in c("terra", "ape", "sf")){
  if(!require(i, character.only = TRUE)){
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}


# species that overlap with leopard range
prekryv = read.table("Results/historickePrekryvyALL.txt", header = TRUE, sep = "\t")
selmy = rownames(prekryv)[prekryv$velkost.prekryvu.m2 > 0]

selmy = sort(c("Panthera_pardus", selmy))


# should niche shift be recalculated?
pocitaj.niche = FALSE

# colors for the likelihood model
res = 60
svetle = 48
farby.interakcie = c(colorRampPalette(c(hcl.colors(res-svetle, "Inferno")[res-svetle], hcl.colors(res-svetle, "Inferno")[res-svetle-1]))(svetle),
  rev(hcl.colors(res-svetle, "Inferno")[1:(res-svetle-1)]))

# colors for the niche shift
res = 60
svetle = 31
farby.nika = c(hcl.colors(svetle, "Viridis"),
          rev(hcl.colors(res-svetle, "Inferno")))



# data
carnivoreImpact = rast("Results/leopard-carnivoreImpact.tiff")
osud = rast("Results/leopard-interaction.tiff")


if(pocitaj.niche){
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
  niche = I
} else {
  niche = read.table("Results/nicheSimilarityI.txt", header = T, sep = "\t")
}


pdf("Results/fig1.pdf", width = 6.2, height = 5)
layout(matrix(c(1,1,2,3,2,3), nrow=2))
# UPGMA from (1 - niche similarity I)
par(mar = c(.5,1.5,1,1) + .2, xpd = NA)
plot(as.phylo(hclust(as.dist(1 - niche), "average")), cex = 1)
mtext("A", side = 3, line = -.5, adj = 0, font = 2)
# niche shift
par(xpd = TRUE)
plot(carnivoreImpact, 
     col = farby.nika,
     pax = list(retro = TRUE), las = 1, plg = list(cex = .7))
mtext("B", side = 3, line = .5, adj = 0, font = 2)

# likelihood model 
plot(osud, 
     col = farby.interakcie,
     pax = list(retro = TRUE), las = 1, plg = list(cex = .7))
mtext("C", side = 3, line = .5, adj = 0, font = 2)
dev.off()



# SDM and leopard range
leopard.klima = rast("Data/SDMpredikcie/Panthera_pardus.tiff")
leopard.selmy = rast("Results/Panthera_pardus.tiff")
leopard.range = read_sf("Data/IUCN/Panthera_pardus/data_0.shp")
land_mask = rast("Data/land_mask.tiff")
# simplify the range polygons
leopard.range.all = as.polygons(rasterize(vect(leopard.range), land_mask, values = 1, background = NA))
# GBIF occurrence points for leopard - thinned
nalezy = read.table("Data/leopard-thinned.csv", sep = ",", header = TRUE)
nalezy = st_as_sf(nalezy, coords = 1:2, crs = 4326)

# correlations
load("Results/interaction-pearson.RData") 
load("Results/interaction-spearman.RData")

# sensitivity
sens = read.table(file = "Results/sensitivity.txt", sep = "\t", header = TRUE)
x = apply(sens, 2, mean, na.rm = T)
n.druhov = length(selmy[-1])
poradie = order(x[seq(n.druhov + 1, 2 * n.druhov, by = 1)], decreasing = FALSE)


pdf("Results/figS1.pdf", width = 8.4, height = 5.6)
layout(matrix(c(1,1,1,3,3,
				1,1,1,3,3,
				1,1,1,3,3,
				2,2,2,4,4,
				2,2,2,4,4,
				2,2,2,5,5), nrow = 5))
par(mar = c(.5,.5,.1,.5)+.2)
plot(leopard.klima, pax = list(retro = TRUE), las = 1, plg = list(cez = .7))
plot(leopard.range.all, add = T, lwd = .7)
mtext("A", side = 3, line = .5, adj = 0, font = 2)
plot(leopard.selmy, pax = list(retro = TRUE), las = 1, plg = list(cez = .7))
plot(nalezy, add = TRUE, pch = ".", cex = .7)
mtext("B", side = 3, line = .5, adj = 0, font = 2)

selmy = rownames(prekryv)[prekryv$velkost.prekryvu.m2 > 0]
par(mar = c(9, 4,1,0)+.2)
  plot(1, type = "n", xlim = c(1, 3 * 24), ylim = c(-1, 1), xlab = "", ylab = "Correlation", axes = FALSE)
  axis(2, las = 1)
  axis(1, las = 2, font = 3, labels = sub("_", " ", selmy), at = seq(1.5, 3 * 24 + .5, by = 3))
  box()
  rect(seq(0, 3 * 24 - 1, by = 6), 
    par("usr")[3], 
    seq(3, 3*24+1.5, by = 6), 
    par("usr")[4], col = "grey", border = NA)
  points(x = seq(1, 3 * 24, by = 3), y = lapply(test.pearson, FUN = \(x) x$estimate), pch = 15)
  points(x = seq(2, 3 * 24, by = 3), y = lapply(test.spearman, FUN = \(x) x$estimate), pch = 17)
   abline(h = 0, lty = 3)
mtext("C", side = 3, line = .5, adj = 0, font = 2)


par(mar = c(4,4,0,0)+.2)
plot(x[seq(1, 2 * n.druhov, by =1)], x[seq(2 * n.druhov + 1, 4 * n.druhov, by = 1)], 
	type = "n", xlab = expression(mu^"*"), ylab = expression(sigma), las = 1, cex.lab = 1)
points(x[seq(1, n.druhov, by = 1)][poradie], x[seq(2 * n.druhov + 1, 3 * n.druhov, by = 1)][poradie], pch = 3, col = hcl.colors(24, "Inferno", rev = TRUE), cex = 1.2, lwd = 1.2)
points(x[seq(n.druhov + 1, 2 * n.druhov, by = 1)][poradie], x[seq(3 * n.druhov + 1, 4 * n.druhov, by = 1)][poradie], pch = 19, col = hcl.colors(24, "Inferno", rev = TRUE), cex = 1.2)

mtext("D", side = 3, line = .5, adj = 0, font = 2)

par(mar = c(0,0,-.3,0) + .5)
plot.new()
legend("topleft", 
	legend = rev(sub("_", " ", selmy)[poradie]), 
	fill = hcl.colors(24, "Inferno"),
	text.font = 3, cex = .8, border = NA, y.intersp = .8)

dev.off()


#### PCA interpretation

load("Results/Panthera_pardus-PCAmodel.RData")
npcs = 15
pdf(paste0("Results/leopard-PCAinterpretation.pdf"), height = 4, width = 10)
layout(matrix(c(rep(1, 20), 2), nrow = 1))
par(mar = c(10,4,0,0) + .1)
image(1:(nrow(pca$rotation) + 1), 1:(npcs + 1), abs(pca$rotation[, 1:npcs]), 
      axes = FALSE, xlab = "", ylab = "")
for(i in 1:nrow(pca$rotation)){
  axis(1, at = i +.5, labels = rownames(pca$rotation)[i], las = 2,
       font = ifelse(i <= 24, 1, 3))
}
axis(2, at = 1.5:(npcs + .5), labels = colnames(pca$rotation)[1:npcs], las = 1)
for(i in 1:nrow(pca$rotation)){
  for(j in 1:npcs){
    if(pca$rotation[i,j] > 0){
      segments(i, j, i + 1, j + 1, col = "grey")
    } else {
      segments(i + 1, j, i, j + 1, col = "grey")
    }
  }
}

par(mar = c(10,0.3,0,2.5)+.1)
cisla = seq(0,max(abs(pca$rotation[, 1:npcs])),length.out = 100)
image(1, cisla, matrix(cisla, nrow = 1), 
      axes = FALSE, col = rev(hcl.colors(100, "YlOrBr")),
      xlab = "", ylab = "")
axis(4, las = 1)
box()

dev.off()
