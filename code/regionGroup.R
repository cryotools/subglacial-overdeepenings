# ###############################
# ## Subglacial Overdeepenings ##
# ###############################
# 
# This is an additional code file for the first work package of the project 'Future glacial lakes
# in High Mountain Asia - Modeling and Risk Analysis' (GLAMoR).
# This script was written by Wilhelm Furian.
# 
# By running this script, a subsequent python routine is able to uniquely identify each single overdeepening. 
# Therefore, it is necessary for all file names to match the names in the python script!
# This script can be used as a faster alternative to the region group function of arcpy, 
# which apparently sometimes breaks when handling large datasets.
# 
# The code is available on github. https://github.com/cryotools/subglacial-overdeepenings
# For more information on this work package see the README.
# For more information on the project as a whole see https://hu-berlin.de/glamor.
# 
# You are allowed to use and modify this code in a noncommercial manner and by
# appropriately citing the above mentioned developer. If you would like to share your own improvements,
# please fork the repository on GitHub, commit your changes, and file a merge request.
# 
# Correspondence: furiawil@geo.hu-berlin.de

library(raster)
input <- list.files(path = "~", 
                    pattern = "sink", recursive = TRUE, full.names = T)
input <- input[!grepl('aux|vat|tfw|ovr|xml', input)] 

fun_regionGroup <- function(input){
  print(input)
  name <- as.character(input)
  r <- raster(name)
  regions <- clump(r) # region group alternative
  out_path = paste0(substr(input[1],1,29),"sinksRegion.tif")
  writeRaster(regions, out_path,format = "GTiff", overwrite = TRUE)
}

for (i in 1:length(input)){
  fun_regionGroup(input[i])
}
