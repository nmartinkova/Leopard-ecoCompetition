##############################################################
#                                                            #
#                          Code for                          #
#         Carnivore interactions shape leopard presence      #
#            Natalia Martinkova, Michal Skrobanek            #
#                                                            #
#   Code author: Natalia Martinkova                          #
#                martinkova@ivb.cz                           #
#                                                            #
#   Installs packages and downloads climatic data            #
#                                                            #
##############################################################


for(i in c( "devtools", "pastclim", "terra", "tidysdm", "tools", "sf", "EFA.dimensions", "xgboost")){
  if(!require(i, character.only = TRUE)){
    if(i %in% c("tidysdm", "terra", "pastclim")){
    switch(i,  
      # dev versions required 
      # https://evolecolgroup.github.io/tidysdm/index.html
      tidysdm = devtools::install_github("EvolEcolGroup/tidysdm"),
      terra = install.packages('terra', repos='https://rspatial.r-universe.dev'),
      pastclim = devtools::install_github("EvolEcolGroup/pastclim", ref="dev")
    )
    } else {
    install.packages(i, dependencies = TRUE)
    }
    library(i, character.only = TRUE)
  }
}


# set up path for climatic data
set_data_path(tools::R_user_dir("pastclim","data"))



# download climatic data
download_dataset("WorldClim_2.1_5m")