# open and format birdlife gdb file
# note: this is just formatting the version Jacob previously accessed. I have 
# just done this conversion myself directly so I can be sure of how it got done. 
# The version I accessed was downloaded 10/08/2021. Can switch to using this once
# I've got the pipeline off the ground
library(sf)

ogrListLayers("birdlife_maps/BOTW.gdb")

ranges <- st_read("birdlife_maps/BOTW.gdb", "All_Species")
saveRDS(ranges, "birdlife_maps/BOTW_extracted.rds")
