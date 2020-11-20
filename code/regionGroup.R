# Fast alternative to region group algorithm in arcpy

library(raster)
# load every single overdeepening in HMA
input <- list.files(path = "~", # this path has to be the same as the work_path from the python script
                    pattern = "sink", recursive = TRUE, full.names = T)
input <- input[!grepl('aux|vat|tfw|ovr|xml', input)] # grep only the .shp file of the overdeepenings

fun_regionGroup <- function(input){
  print(input)
  name <- as.character(input)
  test <- raster(name)
  regions <- clump(test) # region group
  out_path = paste0(substr(input[1],1,29),"sinksRegion.tif")
  writeRaster(regions, out_path,format = "GTiff", overwrite = TRUE)
}

for (i in 1:length(input)){
  fun_regionGroup(input[i])
}