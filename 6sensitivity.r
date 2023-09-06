##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Sensitivity analysis for the extinction rates in the     #
#   likelihood model (eq. 1)                                 #
#                                                            #
##############################################################



for(i in c("terra", "sensitivity")){
  if(!require(i, character.only = TRUE)){
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}


# likelihood model (eq. 1) for one species
# x[1:2] probability of occurrence of leopard and the other species
# x[3:4] reciprocal extinction rates, 3 = other species, 4 = leopard

prezitie = function(x){
	x = unlist(x)
	p = (x[1] * x[2] * (1 - x[3])) - (x[1] * x[2] * x[4]) 
	unlist(log10(ifelse(p <= 0, 1e-10, p)))
}



# species that overlap with leopard range
prekryv = read.table("Results/historickePrekryvyALL.txt", header = T, sep = "\t")
druhy = rownames(prekryv)[prekryv$velkost.prekryvu.m2 > 0]


# SDM predictions for all species
selmy = rast(paste0("Data/SDMpredikcie/Panthera_pardus.tiff"))

names(selmy) = "Panthera_pardus"

for(druh in druhy[1:length(druhy)]){
	druh.pred = rast(paste0("Data/SDMpredikcie/", druh, ".tiff"))
	names(druh.pred) = druh
	add(selmy) <- druh.pred
}

selmy = as.data.frame(selmy)
n.druhov = length(druhy)

# subsampling the Afro-Eurasian raster cells (stochastic process, set seed for debugging)
n = 10000
bunky = sample(1:nrow(selmy), size = n)

res = matrix(NA, nrow = n, ncol = n.druhov * 6, 
	dimnames = list(NULL, paste0(rep(colnames(selmy)[-1], 6), rep(c(".muStar", "L.muStar", ".sigma", "L.sigma", "mu", "L.mu"), each = n.druhov))))

# row counter for results	
riadok = 1	

for(k in bunky){
	# likelihood model (eq. 1) with fixed a_i to a values for a specific raster cell
	osud = function(X){
		n = length(X) / 2
		res = 0	
		for(i in seq(1, n, by = 1)){
			res = res + prezitie(c(selmy[k, 1], selmy[k, i], X[i], X[n + i]))
		}
		unname(unlist(res))
	}
	
	# design for the Morris's test
	morris.test = morris(osud, 
		factors = n.druhov * 2,
		r = 15,
		design = list(type = "oat", levels = 30, grid.jump = 15)
	)
	
	# supply model predictions
	y = apply(morris.test$X, 1, osud)
	tell(morris.test, y = y)
	
	# calculate Morris's mu* and sigma
	mu.star <- apply(morris.test$ee, 2, function(x) mean(abs(x)))
	sigma <- apply(morris.test$ee, 2, sd)
	mu = apply(morris.test$ee, 2, function(x) mean(x))
	
	# save results for each selected raster cell
	res[riadok, seq(1, n.druhov * 2, by = 1)] = mu.star
	res[riadok, seq(n.druhov * 2 + 1, 4 * n.druhov, by = 1)] = sigma
	res[riadok, seq(n.druhov * 4 + 1, ncol(res), by = 1)] = mu
 	riadok = riadok + 1
 	write.table(res, file = "Results/sensitivity.txt", sep = "\t")

	x = apply(res, 2, mean, na.rm = T)

	poradie = order(x[seq(n.druhov+1, 2*n.druhov, by = 1)], decreasing = F)

	pdf("Results/sensitivity.pdf", width = 7.3, height = 3.9)
	layout(matrix(c(1,1,1,2), nrow = 1))
	par(mar = c(4,4,0,0)+.2)
	plot(x[seq(1, 2* n.druhov, by =1)], x[seq(2*n.druhov+1, 4*n.druhov, by = 1)], 
		type = "n", xlab = expression(mu^"*"), ylab = expression(sigma), las = 1, cex.lab = 1.3)
	points(x[seq(1, n.druhov, by = 1)][poradie], x[seq(2*n.druhov+1, 3*n.druhov, by = 1)][poradie], pch = 1, col = hcl.colors(24, "Inferno", rev = T), cex = 1.5, lwd = 1.2)
	points(x[seq(n.druhov+1, 2*n.druhov, by = 1)][poradie], x[seq(3*n.druhov+1, 4*n.druhov, by = 1)][poradie], pch = 19, col = hcl.colors(24, "Inferno", rev = T), cex = 1.5)
	par(mar = c(0,0,-.3,0)+.5)
	plot.new()
	legend("topleft", 
		legend = rev(sub("_", " ", names(selmy)[-1])[poradie]), 
		fill = hcl.colors(24, "Inferno"),
		text.font  =3)
	dev.off()


	pdf("Results/sensitivityMU.pdf", width = 7.3, height = 3.9)
	layout(matrix(c(1,1,1,2), nrow = 1))
	par(mar = c(4,4,0,0)+.2)
	plot(x[seq(4*n.druhov+1, 6* n.druhov, by =1)], x[seq(2*n.druhov+1, 4*n.druhov, by = 1)], 
		type = "n", xlab = expression(mu), ylab = expression(sigma), las = 1, cex.lab = 1.3)
	points(x[seq(4*n.druhov+1, 5*n.druhov, by = 1)][poradie], x[seq(2*n.druhov+1, 3*n.druhov, by = 1)][poradie], pch = 1, col = hcl.colors(24, "Inferno", rev = T), cex = 1.5, lwd = 1.2)
	points(x[seq(5*n.druhov+1, 6*n.druhov, by = 1)][poradie], x[seq(3*n.druhov+1, 4*n.druhov, by = 1)][poradie], pch = 19, col = hcl.colors(24, "Inferno", rev = T), cex = 1.5)
	par(mar = c(0,0,-.3,0)+.5)
	plot.new()
	legend("topleft", 
		legend = rev(sub("_", " ", names(selmy)[-1])[poradie]), 
		fill = hcl.colors(24, "Inferno"),
		text.font  =3)
	dev.off()

}

# save final result for all n raster cells, row names correspond to non-NA cells in land_mask
rownames(res) = bunky
write.table(res, file = "Results/sensitivity.txt", sep = "\t")


