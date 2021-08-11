# create dataframe of habitat associations for crude subsetting (i.e. species 
# that do not occupy marine environments and have ranges that intersect with the
# tropics)

library(dplyr)

# scraped bird traits
load("data/birdlife_traits.Rdata")

hab_1 <- lapply(birdlife_traits$habitats, function(x) unique(sub("(^.*?);.*", "\\1", x))) 
hab_2 <- lapply(birdlife_traits$habitats, function(x) unique(sub("^.*?; (.*?);.*", "\\1", x)))

hab_1_df <- lapply(hab_1, as_tibble) %>%
    bind_rows(., .id="species") %>%
    rename(habitat_1 = value) %>%
    group_by(species) %>%
    mutate(habitat_id = paste0("habitat_", 1:n())) 

hab_2_df <- lapply(hab_2, as_tibble) %>%
    bind_rows(., .id="species") %>%
    rename(habitat_2 = value) %>%
    group_by(species) %>%
    mutate(habitat_id = paste0("habitat_", 1:n())) 


both <- left_join(hab_1_df, hab_2_df) %>%
    select(species, habitat_id, everything())

# 3463 species in the tropics
both %>% 
    group_by(species) %>%
    summarise(is_tropical = any(grepl("Tropical|Subtropical", habitat_2))) %>%
    summarise(sum(is_tropical))

# 3391 in the tropics & not marine
both %>% 
    group_by(species) %>%
    filter(!any(grepl("Marine", habitat_1))) %>%
    summarise(is_tropical = any(grepl("Tropical|Subtropical", habitat_2))) %>%
    summarise(sum(is_tropical))

saveRDS(both, "data/habitat_associations.rds")
