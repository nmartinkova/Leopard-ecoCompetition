##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Calculates reciprocal extinction rates from IUCN range   #
#   overlap                                                  #
#                                                            #
##############################################################


for(i in c("sf")){
  if(!require(i, character.only = TRUE)){
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}


# turn off spherical geometry to avoid invalid polygons after intersecting them
sf_use_s2(FALSE)

# species names
druhy = dir("Data/IUCN", include.dirs = TRUE)
druhy = druhy[!(druhy %in% c("Panthera_pardus"))]

leopard = read_sf("Data/IUCN/Panthera_pardus/data_0.shp")

# https://nc.iucnredlist.org/redlist/resources/files/1539614211-Mapping_attribute_codes_v1.16_2018.pdf
# 1 Extant The species is known or thought very likely to occur currently in the area, which encompasses localities with current or recent (last 20-30 years) records where suitable habitat at appropriate altitudes remains. Extant ranges should be considered in the calculation of EOO. When mapping an “assisted colonisation” it is important to note that this range should be treated as Extant.
# 2 Probably Extant This code value has been discontinued for reasons of ambiguity. It may exist in the spatial data but will gradually be phased out.
# 3 Possibly Extant There is no record of the species in the area, but the species may possibly occur, based on the distribution of potentially suitable habitat at appropriate altitudes, although the area is beyond where the species is Extant (i.e., beyond the limits of known or likely records), and the degree of probability of the species occurring is lower (e.g., because the area is beyond a geographic barrier, or because the area represents a considerable extension beyond areas of known or probable occurrence). Identifying Possibly Extant areas is useful to flag up areas where the taxon should be searched for. Possibly Extant ranges should not be considered in the calculation of EOO.
# 4 Possibly Extinct The species was formerly known or thought very likely to occur in the area (post 1500 AD), but it is most likely now extirpated from the area because habitat loss and/or other threats are thought likely to have extirpated the species, and there have been no confirmed recent records despite searches. Possibly Extinct ranges should not be considered in the calculation of EOO.
# 5 Extinct The species was formerly known or thought very likely to occur in the area (post 1500 AD), but it has been confirmed that the species no longer occurs because exhaustive searches have failed to produce recent records, and the intensity and timing of threats could plausibly have extirpated the taxon. Extinct ranges should not be considered in the calculation of EOO.
# 6 Presence Uncertain

# velkost.prekryvu.m2 = range overlap in m2
# oba.vyhynuli = range where both species are extinct or possibly extinct (4,5)
# oba.pritomni = range where both species are extant (1)
#   Note that probably extant and possibly extant (2,3) are ignored
# leopard.vyhynul = range where leopard is extinct (4,5) but the other species is extant (1)
# druh.vyhynul = range where the other species is extinct (4,5) but the leopard is extant (1)
historicke.prekryvy = data.frame(matrix(NA, 
  nrow = length(druhy), ncol = 5,
  dimnames = list(druhy, c("velkost.prekryvu.m2", "oba.vyhynuli", "oba.pritomni", "leopard.vyhynul", "druh.vyhynul"))))


# calculate area of range overlap in m2
vypocitaj.plochu = function(x,y){
	if(any(nrow(x) == 0, nrow(y) == 0)){
		return(0)
	} else {
		as.numeric(sum(st_area(st_intersection(st_make_valid(st_geometry(x)), st_make_valid(st_geometry(y))))))
	}
}


for(i in 1:length(druhy)){
  druh = druhy[i]
  cat(druh, "\n")
  areal = st_make_valid(read_sf(paste0("Data/IUCN/", druh, "/data_0.shp")))
  prekryv.arealov = st_intersection(leopard, areal)
  vyhynuly.leopard = st_collection_extract(
    st_intersection(prekryv.arealov, leopard[leopard$PRESENCE %in% c(4:5),]),
    "POLYGON")
  pritomny.leopard = st_collection_extract(
    st_intersection(prekryv.arealov, leopard[leopard$PRESENCE == 1,]),
    "POLYGON")
  vyhynuly.druh = st_collection_extract(
    st_intersection(prekryv.arealov, areal[areal$PRESENCE %in% c(4:5),]),
    "POLYGON")
  pritomny.druh = st_collection_extract(
    st_intersection(prekryv.arealov, st_make_valid(areal)[areal$PRESENCE == 1,]),
    "POLYGON")
  historicke.prekryvy[i,] = c(
    # velkost.preryvy.m2
    as.numeric(sum(st_area(st_geometry(prekryv.arealov)))),
    # oba.vyhynuli
    vypocitaj.plochu(vyhynuly.leopard, vyhynuly.druh),
    # oba.pritomni
    vypocitaj.plochu(pritomny.leopard, pritomny.druh),
    # leopard.vyhynul
    vypocitaj.plochu(vyhynuly.leopard, pritomny.druh),
    # druh vyhynul
    vypocitaj.plochu(pritomny.leopard, vyhynuly.druh)
  )
  # save preliminary table after each species is calculated
  write.table(historicke.prekryvy, "Results/historickePrekryvy.txt", sep = "\t")
}

# from m2 to extinction rate proxy
historicke.prekryvy[,2:5] = historicke.prekryvy[,2:5] / historicke.prekryvy[,1]
write.table(historicke.prekryvy, "Results/historickePrekryvy.txt", sep = "\t")