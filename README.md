# Leopard-ecoCompetition

**Natália Martínková and Michal Škrobánek**


R code to replicate analyses from the Leopard-ecoCompetition paper. Required data must be supplied in a `Data/` folder in the working directory.
* `Data/GBIF` - tab-delimited files from [GBIF](https://www.gbif.org). Note that files downloaded from GBIF will have a csv suffix, despite being tab-delimited. File names are expected to be species names. Geographic coordinates are expected to be in columns `decimalLongitude` and `decimalLatitude`.
* `Data/IUCN` - shapefiles from the [IUCN Red List](https://www.iucnredlist.org), one folder per species. Folder names must correspond to those used in `Data/GBIF`.
* `Data/SEDAC` - rasters in TIFF format with landscape environmental layers from the [SEDAC](https://sedac.ciesin.columbia.edu) database.

The code should be run in sequence. The results might slightly differ from those included in the paper as the code includes multiple stochastic processes that could influence the results.

![leopard-PCAinterpretacia.pdf](https://github.com/nmartinkova/Leopard-ecoCompetition/files/12537277/leopard-PCAinterpretacia.pdf)
