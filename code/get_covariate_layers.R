# Export elevation and forest cover rasters for Eastern Cordillera for plotting
# maps. 

## Packages ----
library(rgee); library(sf); library(stars); library(dplyr); 
library(rnaturalearth)
cols <- function(x) colorRampPalette(c('#bbe029', '#0a9501', '#074b03'))(x-1)

## Set up rgee session ----
ee_Initialize()

## EE datasets ----
ALOS <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")
tc_modis <- ee$ImageCollection("MODIS/006/MOD44B")
worldclim <- ee$Image("WORLDCLIM/V1/BIO")

# rectangle designating neotropics ----
geometry <- ee$Geometry$Rectangle(
    coords = c(-113,-23.5,-34,23.5),
    proj = "EPSG:4326",
    geodesic = FALSE
)

# get treecover in 2010
pct_tc <- tc_modis$select('Percent_Tree_Cover')
listOfImages <- pct_tc$toList(pct_tc$size())
tc_2010 <- ee$Image(listOfImages$get(10))

# export for neotropics 10k resolution, WGS84 projection
tc_stars <- ee_as_stars(tc_2010, region=geometry, scale = 50000, crs="epsg:4326")
write_stars(tc_stars, "data/covariate_layers/treecover_modis_2010_50km.tif")


# mask out non-land
tc_stars_local <- read_stars("data/covariate_layers/treecover_modis_2010.tif")
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_union(world)

tc_masked <- tc_stars_local[world]

# plot and save
png("figures/treecover_modis_2010.png", width=180, height=100, units="mm", res=200)
plot(tc_masked, nbreaks=10, breaks="equal", col=cols(10), reset=F)
plot(st_geometry(world), add=T)
dev.off()

write_stars(tc_masked, "data/covariate_layers/treecover_modis_2010_masked.tif")

# Worldclim data ----
worldclim_sub <- worldclim$select(opt_selectors = c("bio01", "bio04", "bio13", 
                                                    "bio14", "bio15"), 
                                  opt_names = c("mat", "temp_seasonality", 
                                                "precip_wmonth", "precip_dmonth", 
                                                "precip_seasonality"))
Map$addLayer(worldclim_sub, list(min=0, max=320), name = "mat")
Map$addLayer(worldclim_sub, list(min=60, max=22721), name="temp_seasonality")

worldclim_export <- ee_as_stars(worldclim_sub, region=geometry, scale = 50000, 
                                crs="epsg:4326")
write_stars(worldclim_export, "data/covariate_layers/worldclim_50km.tif")


worldclim_local <- read_stars("data/covariate_layers/worldclim.tif")

# plot worldclim data ----
worldclim_spl <- split(worldclim_local, "band")
layer_names <- names(worldclim_spl)
plot_worldclim <- function(x) plot(worldclim_spl[x], nbreaks=10, 
                                   col=viridis::cividis(9), reset=F, key.pos=NULL)

png("figures/worldclim.png", width=300, height=200, units="mm", res=300)
par(mfrow=c(2, 3))
lapply(layer_names, plot_worldclim)
dev.off()

# check for identical CRS ----
identical(st_crs(tc_stars_local), st_crs(worldclim_local))
dim(tc_stars_local)
dim(worldclim_local)
