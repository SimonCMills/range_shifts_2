## rasterise range polygons
# generates a final dataframe with presence/absence for all land cells in 
# neotropics

# packages ----
library(stars); library(sf); library(dplyr)

# climate layers for neotropics
clim <- read_stars("data/covariate_layers/worldclim_50km.tif")

# polygons for neotropical species
range_maps <- readRDS("birdlife_maps/range_maps_neotropical_endemics.rds") %>%
    select(-is_valid)

grd <- clim[,,,1]
grd[[1]] <- NA

rasts <- rep(NA, nrow(range_maps)) %>% as.list
names(rasts) <- range_maps$SCINAME
for(i in 1:nrow(range_maps)) {
    rasts[[i]] <- st_rasterize(range_maps[i,], grd, options=c("ALL_TOUCHED=T"))
}

saveRDS(rasts, "birdlife_maps/range_maps_neotropical_endemics_rast50k.rds")


clim_df <- as.data.frame(clim) %>%
    reshape2::dcast(., x+y~band) %>%
    filter(complete.cases(.))

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

# get buffered cells ----
# updates range_maps_df sequentially
ranges_buffered <- range_maps %>%
    mutate(Shape = st_buffer(Shape, 5e5))

rasts_buffered <- rep(NA, nrow(range_maps)) %>% as.list
names(rasts_buffered) <- ranges_buffered$SCINAME
for(i in 1:nrow(ranges_buffered)) {
    rasts_buffered[[i]] <- st_rasterize(ranges_buffered[i,], grd, options=c("ALL_TOUCHED=T"))
}

r_buffered_df <- lapply(rasts_buffered, function(x) as.data.frame(x) %>%
                   filter(complete.cases(.))) %>% 
    bind_rows(., .id="species")

range_maps_df <- left_join(range_maps_df, r_buffered_df) %>%
    rename(within_500km = ID) %>%
    mutate(within_500km = ifelse(is.na(within_500km), 0, 1))

range_maps_df %>% 
    filter(species == "Abeillia abeillei") %>%
    ggplot(aes(x, y, fill=factor(presence + within_500km))) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_d()

## repeat but for 1000km 
ranges_buffered <- range_maps %>%
    mutate(Shape = st_buffer(Shape, 1e6))

rasts_buffered <- rep(NA, nrow(range_maps)) %>% as.list
names(rasts_buffered) <- ranges_buffered$SCINAME
for(i in 1:nrow(ranges_buffered)) {
    rasts_buffered[[i]] <- st_rasterize(ranges_buffered[i,], grd, options=c("ALL_TOUCHED=T"))
}

r_buffered_df <- lapply(rasts_buffered, function(x) as.data.frame(x) %>%
                            filter(complete.cases(.))) %>% 
    bind_rows(., .id="species")

range_maps_df <- left_join(range_maps_df, r_buffered_df) %>%
    rename(within_1000km = ID) %>%
    mutate(within_1000km = ifelse(is.na(within_1000km), 0, 1))

range_maps_df %>% 
    filter(species == "Abeillia abeillei") %>%
    ggplot(aes(x, y, fill=factor(presence + within_500km + within_1000km))) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_d()

#
saveRDS(range_maps_df, "birdlife_maps/neotropical_endemics_df_50km_with_buffers.rds")
