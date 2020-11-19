input <- list.files(path = "~", # work_path from python script
                    pattern = "sink", recursive = TRUE, full.names = T)
input <- input[!grepl('aux|vat|tfw|ovr|xml', input)]

fun_regionGroup <- function(input){
  print(input)
  name <- as.character(input)
  test <- raster(name)
  regions <- clump(test) # alternative to region group routine from arcpy
  out_path = paste(substr(input[1],1,29),"sinksRegion.tif", sep="")
  writeRaster(regions, out_path,format = "GTiff", overwrite = TRUE)
}

for (i in 1:length(input)){
  fun_regionGroup(input[i])
}


  
  
  