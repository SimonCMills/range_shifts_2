# clean range polygons
# a number of polygons are not valid, and need some rebuilding to get them to be
# Relies on having previous extracted the range maps (in extract_birdlife_ranges.R)

# note: there's a really bizarre issue that when I select by within(tropics) rather
# than within(neotropics) it returns no species, despite the fact that the neotrop
# polygon is a subset of the tropics polygon. For now, just subset for neotropical
# species, but it needs resolving. 

# packages ----
library(dplyr); library(sf); library(stars); library(rnaturalearth)
library(rnaturalearthdata)


# datasets & polygons ----
# coarse habitat associations (for removing any marine birds)
habitat_assoc <- readRDS("data/habitat_associations.rds")

species_subset <- habitat_assoc %>% 
    group_by(species) %>%
    filter(!any(grepl("Marine", habitat_1))) %>%
    filter(any(grepl("Tropical|Subtropical", habitat_2))) %>%
    pull(species) %>%
    unique

# polygons for tropics and neotropics
neotrop <- matrix(c(-113, 23.5, -113, -23.5, -34, -23.5, -34, 23.5, -113, 23.5),
                  byrow=T, ncol=2) %>%
    list() %>%
    st_polygon() %>%
    st_sfc()
st_crs(neotrop) <- st_crs("WGS84")

# tropics <- matrix(c(-179, 23.45, -179, -23.45, 1, -23.45, 1, 23.45, -179, 23.45),
#                   byrow=T, ncol=2) %>%
#     list() %>%
#     st_polygon() %>%
#     st_sfc()
# st_crs(tropics) <- st_crs("WGS84")

world <- ne_countries(scale = "medium", returnclass = "sf")
map_neotrop <- st_intersection(world, neotrop) %>% summarise()

# Subset ranges ----
## Remove non-tropical species
# read in range maps and drop species that are at least partially marine
range_maps <- readRDS("birdlife_maps/BOTW_extracted.rds") %>% 
    filter(SCINAME %in% species_subset) 
class(range_maps) <- c("sf", "tbl_df", "tbl", "data.frame")

## Fix invalid polygons ----
check_valid <- range_maps %>%
    mutate(is_valid = st_is_valid(Shape))

# 2 species have invalid geometry
check_valid %>% 
    filter(!is_valid)

# rebuild 2 species with invalid geometries
fixed_geoms <- check_valid %>%
    filter(!is_valid) %>%
    st_as_sf() %>%
    st_as_s2(., rebuild=T) %>%
    s2::s2_buffer_cells(., 0) %>%
    st_as_sf 

clean_polys <- bind_rows(check_valid %>% filter(is_valid==TRUE), 
          fixed_geoms)

## Get species in neotropics ----
# now want to know: is it entirely within the tropics. 
# union all polygons associated with a species
union_polys <- clean_polys %>%
    group_by(SCINAME) %>%
    summarise(Shape = st_union(Shape)) %>%
    st_make_valid()

saveRDS(union_polys, "birdlife_maps/range_maps_unioned_intsct_tropics.rds")

## only want species fully within the neotropics. 
within <- union_polys %>%
    filter(st_within(Shape, neotrop, sparse = FALSE))

saveRDS(within, "birdlife_maps/range_maps_neotropical_endemics.rds")

# Okay, for future: ----
# - write code for assigning continental occupancy (e.g. neotropics, Africa, 
# East Asia?). Need this for selecting climate layers (doesn't make sense to 
# include climate from a continent a species doesn't exist in)
# 
# neotropics <- ne_countries(scale="large", returnclass = "sf", continent = c("North America", "Oceania")) %>%
#     summarise(geometry = st_union(geometry))
# 
# neotropics
# test
# as_tibble(test)
# library(ggplot2)
# ggplot() + geom_sf(data=union_polys[1,]) +
#     geom_sf(data=tropics, fill="red", alpha=.1) +
#     geom_sf(data=neotrop, fill="blue", alpha=.1)
# 
# # strange buffering issue.. 
# test <- st_buffer(union_polys[1,], 9000000)# %>%
# test2 <- st_buffer(union_polys[1,], 1000000)
# test3 <- st_buffer(union_polys[1,], 5000000)
# test4 <- st_buffer(union_polys[1,], 8000000)
# 
# 
# plot(test3)
# plot(test4)
# 
# ggplot(test) + geom_sf(fill="red")     +
#     geom_sf(data=test2) +
#     geom_sf(data=test3, alpha=.1) +
#     geom_sf(data=test4, alpha=.1, fill="blue")
# 
# ggplot(union_polys[1,]) + geom_sf()
# summarise(st_within(Shape, tropics, sparse=FALSE))
#
## JUNK: old code, delete at some point
# within_1 <- ranges_1 %>%
#     filter(st_within(Shape, neotrop, sparse = FALSE))
# 
# within_2 <- ranges_2 %>%
#     filter(st_within(Shape, neotrop, sparse = FALSE))
# 
# intsct_1 <- ranges_1 %>%
#     filter(st_intersects(Shape, neotrop, sparse = FALSE))
# 
# intsct_2 <- ranges_2 %>%
#     filter(st_intersects(Shape, neotrop, sparse = FALSE))
# 
# 
# within_both <- bind_rows(within_1, within_2)
# intsct_both <- bind_rows(intsct_1, intsct_2)
# 
# # save polygons ----
# saveRDS(within_both, "polygons_neotropical_endemics.rds")
# saveRDS(intsct_both, "polygons_neotropical_overlap.rds")
# 
# # this takes a *very* long time
# intsct_valid <- st_make_valid(intsct_both)
# saveRDS(intsct_valid, "polygons_valid_neotropical_overlap.rds")
# st_write(intsct_valid, "polygons_valid_intersect_neotropics.shp")
# 
# # produce clipped set of polygons
# clipped <- st_intersection(intsct_valid, neotrop)
# saveRDS(clipped, "polygons_valid_clipped_neotropical_overlap.rds")
# 
# # create grid for calculating species richness
# grd <- st_as_stars(st_bbox(neotrop), dx = .1, dy=.1)#nx = 200, ny = 200)
# 
# test <- clipped %>% mutate(sp = 1)
# rast <- st_rasterize(test["sp"], grd, options = c("MERGE_ALG=ADD"))
# rast_clipped <- rast[map_neotrop]
# 
# saveRDS(rast, "raster_neotropical_SpRic.rds")
# saveRDS(rast_clipped, "raster_neotropical_SpRic_clipped.rds")
# 
# 
# # plot
# png("Species_richness.png", width=300, height=150, units="mm", res=200)
# n_breaks <- 20
# plot(rast_clipped, nbreaks = n_breaks, breaks="equal", col = hcl.colors(n_breaks-1, "Spectral"), 
#      main="", key.pos = 4, reset = FALSE)
# title(main="Species richness in .1 deg cell", adj=0)
# plot(map_neotrop["geometry"], add=T)
# dev.off()