## rasterise range polygons
# notes: ultimately generates a dataframe that only contains presences for 
# land (i.e. drops ocean cells). Also 
# need to create a buffer raster that contains cells within some large radius 
# that can be treated as absences in modelling. 

# packages ----
library(stars); library(sf); library(dplyr)

# climate layers for neotropics
clim <- read_stars("data/covariate_layers/worldclim_50km.tif")

# polygons for neotropical species
range_maps <- readRDS("birdlife_maps/range_maps_neotropical_endemics.rds") %>%
    select(-is_valid)

grd <- clim[,,,1]
grd[[1]] <- NA

rast <- st_rasterize(sf = range_maps[1:5,], template = grd)

rasts <- rep(NA, nrow(range_maps)) %>% as.list
names(rasts) <- range_maps$SCINAME
for(i in 1:nrow(range_maps)) {
    rasts[[i]] <- st_rasterize(range_maps[i,], grd, options=c("ALL_TOUCHED=T"))
}

saveRDS(rasts, "birdlife_maps/range_maps_neotropical_endemics_rast50k.rds")


clim_df <- as.data.frame(clim) %>%
    reshape2::dcast(., x+y~band) %>%
    filter(complete.cases(.))

as.data.frame(rasts[[1]]) %>%
reshape2::dcast(., x+y~.) 

r_df <- lapply(rasts, function(x) as.data.frame(x) %>%
           filter(complete.cases(.)))

# get land ----
range_maps_df <- lapply(r_df, function(x) {
    left_join(clim_df[c("x", "y")], x, by=c("x", "y")) %>% 
        rename(presence = ID) %>%
        mutate(presence = ifelse(is.na(presence), 0, presence))
}) %>%
    bind_rows(., .id="species") %>% 
    as_tibble 

ggplot(range_maps_df %>% filter(species == "Grallaria quitensis"), aes(x, y)) +
    geom_tile(aes(fill=factor(presence))) +
    coord_equal()


saveRDS(range_maps_df, "birdlife_maps/neotropical_endemics_df_50km.rds")
saveRDS(clim_df, "data/climate_df_50km.rds")
