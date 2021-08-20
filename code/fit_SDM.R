# initial exploration

# packages
library(dplyr); library(sf); library(mgcv); library(ggplot2); library(stars)
library(raster); library(spdep); library(rgdal); library(sp); library(INLA)

# functions ----
sortBase <- function(vec, n.knots = 2, prev_base) {
    ## Function to calculate bases for regression splines. Modified from code
    ## provided in Crainiceanu, C., Ruppert, D. & Wand, M.P. Bayesian analysis for
    ## penalized spline regression using WinBUGS. J. Stat. Soft. 14, 1?24(2005).
    ## Parameter vec is a vector defining the raw data vector, n.knots defines the
    ## number of knots in the GAM.
    N         <- length(vec)
    x.time    <- c(vec)
    ZFE       <- cbind(rep(1,N), x.time)
    if(missing(prev_base)) {
        x.knots   <- quantile(unique(x.time), seq(0, 1, length = (n.knots+2))[-c(1,
                                                                                 (n.knots+2))], na.rm = TRUE)
    } else {
        x.knots <- attr(prev_base, "knot_pos")
        n.knots <- length(x.knots)
    }
    Z_K       <- (abs(outer(x.time,x.knots,"-")))^3
    OMEGA.all <- (abs(outer(x.knots,x.knots,"-")))^3
    svd.OMEGA.all  <- svd(OMEGA.all)
    sqrt.OMEGA.all <- t(svd.OMEGA.all$v %*% (t(svd.OMEGA.all$u) *
                                                 sqrt(svd.OMEGA.all$d)))
    Z.out     <- t(solve(sqrt.OMEGA.all, t(Z_K)))
    attr(Z.out, "knot_pos") <- x.knots
    return(Z.out)
}
scaling_fun <- function(x, center, scale) {
    y <- x
    for(i in 1:ncol(x)) {
        y[,i] <- (x[,i] - center[i])/scale[i]
    }
    y
}

# sort layers for plotting ----
# neotropics
neotrop <- matrix(c(-113, 23.5, -113, -23.5, -34, -23.5, -34, 23.5, -113, 23.5),
                  byrow=T, ncol=2) %>%
    list() %>%
    st_polygon() %>%
    st_sfc()
st_crs(neotrop) <- st_crs("WGS84")

countries <- rnaturalearth::ne_countries(continent=c("South America", "North America"), 
                                         returnclass = "sf")  %>% 
    st_intersection(., neotrop)

# datasets ----
df <- readRDS("birdlife_maps/neotropical_endemics_df_50km_with_buffers.rds")
ranges <- readRDS("birdlife_maps/range_maps_neotropical_endemics.rds")
ranges_rast <- readRDS("birdlife_maps/range_maps_neotropical_endemics_rast50k.rds")
clim <- readRDS("data/climate_df_50km.rds")

# species names
species_list <- unique(df$species)

# species_i_name <- "Psophia crepitans" # Grey-winged trumpeter # pale-winged south west
species_i_name <- "Pipra filicauda" # Wire-tailed manakin
# species_i_name <- "Jacamerops aureus"
#species_i_name <- "Microrhopias quixensis" # Dot-winged antwren
# species_i_name <- "Henicorhina leucophrys" 
# species_i_name <- "Andigena nigrirostris"
# species_i_name <- "Adelomyia melanogenys" 
species_i_name <- "Rhegmatorhina gymnops"
species_i <- which(species_list == species_i_name)
print(species_i)
# species_list[grepl("Campylor", species_list)]

# select single species ----
subset_rast <- ranges_rast[[species_i]]
st_rangemap_pixels <- st_union(st_as_sf(subset_rast))
plot(st_rangemap_pixels)
subset <- df %>% 
    filter(species == species_list[species_i], #[999], 
           within_1000km == 1) %>%
    left_join(., clim) %>%
    mutate(mat_std = scale(mat), 
           temp_seasonality_std = scale(temp_seasonality), 
           precip_wmonth_std = scale(precip_wmonth)) %>%
    arrange(desc(y), x) %>% 
    as.data.frame

# create spatial matrix ----
spg <- subset
# spg <- st_as_sf(subset[c("presence", "x", "y")], coords=c("x", "y"))
# st_crs(spg) <- "epsg:4326"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
r <- raster(spg)
r2 <- stars::st_as_stars(r)
poly_subset <- st_as_sf(r2)

nb_subset <- poly2nb(poly_subset, queen=F)
adj <- nb2mat(nb_subset, style="B", zero.policy = T)
adj <- as(adj, "dgTMatrix")

# sort out climate variables ----
var_names <- c("Mean annual temp", "Temp seasonality", "Precip: driest month", 
               "Precip: wettest month")
subset$idx <- 1:nrow(subset)
temp1 <- scale(sortBase(subset$mat, n.knots = 2))
temp2 <- scale(sortBase(subset$temp_seasonality, n.knots = 2))
precip1 <- scale(sortBase(subset$precip_dmonth, n.knots = 2))
precip2 <- scale(sortBase(subset$precip_wmonth, n.knots = 2))
# long <- scale(sortBase(scale(subset$x)))
# lat <- scale(sortBase(scale(subset$y)))

# sort out variables for predicting margin
temp1_pred <- with(subset, seq(min(mat), max(mat), len=30))
temp2_pred <- with(subset, seq(min(temp_seasonality), max(temp_seasonality), len=30))
precip1_pred <- with(subset, seq(min(precip_dmonth), max(precip_dmonth), len=30))
precip2_pred <- with(subset, seq(min(precip_wmonth), max(precip_wmonth), len=30))

temp1_pred_sc <- scaling_fun(sortBase(temp1_pred, prev_base = temp1), 
                             attr(temp1, "scaled:center"), 
                             attr(temp1, "scaled:scale"))

temp2_pred_sc <- scaling_fun(sortBase(temp2_pred, prev_base = temp2), 
                             attr(temp2, "scaled:center"), 
                             attr(temp2, "scaled:scale"))

precip1_pred_sc <- scaling_fun(sortBase(precip1_pred, prev_base = precip1), 
                             attr(precip1, "scaled:center"), 
                             attr(precip1, "scaled:scale"))

precip2_pred_sc <- scaling_fun(sortBase(temp1_pred, prev_base = precip2), 
                             attr(precip2, "scaled:center"), 
                             attr(precip2, "scaled:scale"))


# make linear combinations ----
lc1 = inla.make.lincombs(`(Intercept)` = rep(1, nrow(temp1_pred_sc)), 
                         temp11 = temp1_pred_sc[,1], 
                         temp12 = temp1_pred_sc[,2])
names(lc1) <- paste0(names(lc1), "_1")

lc2 = inla.make.lincombs(`(Intercept)` = rep(1, nrow(temp2_pred_sc)), 
                         temp21 = temp2_pred_sc[,1], 
                         temp22 = temp2_pred_sc[,2])
names(lc2) <- paste0(names(lc2), "_2")

lc3 = inla.make.lincombs(`(Intercept)` = rep(1, nrow(precip1_pred_sc)), 
                         precip11 = precip1_pred_sc[,1], 
                         precip12 = precip1_pred_sc[,2])
names(lc3) <- paste0(names(lc3), "_3")

lc4 = inla.make.lincombs(`(Intercept)` = rep(1, nrow(precip2_pred_sc)), 
                         precip21 = precip2_pred_sc[,1], 
                         precip22 = precip2_pred_sc[,2])
names(lc4) <- paste0(names(lc4), "_4")

lc5 = inla.make.lincombs(`(Intercept)` = rep(1, nrow(temp1)), 
                         temp11 = temp1[,1], 
                         temp12 = temp1[,2], 
                         temp21 = temp2[,1], 
                         temp22 = temp2[,2], 
                         precip11 = precip1[,1], 
                         precip12 = precip1[,2],
                         precip21 = precip2[,1], 
                         precip22 = precip2[,2])
names(lc5) <- paste0(names(lc5), "_5")

lc6 <-  inla.make.lincombs(`(Intercept)` = rep(1, nrow(temp1)), 
                           idx = subset$idx, 
                           temp11 = temp1[,1], 
                           temp12 = temp1[,2])
names(lc6) <- paste0(names(lc6), "_6")

lc_all <- c(lc1, lc2, lc3, lc4, lc5)

# fit ----
m1 <- inla(presence ~ 1 + temp1 + temp2 + precip1 + precip2, 
           data=subset, family="binomial", control.predictor = list(compute=T), lincomb=lc_all)

m2 <- inla(presence ~ 1 +  temp1 + temp2 + precip1 + precip2 + f(idx, model="besag", graph=adj), 
           data=subset, family="binomial", control.predictor = list(compute=T), 
           control.compute = list(waic = TRUE), 
           lincomb=lc_all)

m2_null <- inla(presence ~ 1 + f(idx, model="besag", graph=adj), 
           data=subset, family="binomial", control.predictor = list(compute=T), 
           control.compute = list(waic = TRUE))

# waics
m2$waic$waic - m2_null$waic$waic

# covariate effects ----
lincomb_summ_1 <- m1$summary.lincomb.derived %>%
    mutate(lc_id = row.names(.)) %>%
    as_tibble %>%
    filter(!grepl("_5$|_6$", lc_id)) %>%
    rename(q2.5 = `0.025quant`,
           q50 = `0.5quant`, 
           q97.5 = `0.975quant`) %>%
    mutate(lc_id2 = as.integer(gsub(".*_(.*)", "\\1", lc_id)),
           var_name = factor(var_names[lc_id2], levels=var_names),
           var = c(temp1_pred/10, temp2_pred/100, precip1_pred, precip2_pred),
           mean = mean - m1$summary.fixed$mean[1], 
           q2.5 = q2.5 - m1$summary.fixed$mean[1], 
           q97.5 = q97.5 - m1$summary.fixed$mean[1]) #%>%
    # mutate_at(c("mean", "q2.5", "q97.5"), boot::inv.logit)

lincomb_summ_2 <- m2$summary.lincomb.derived %>%
    mutate(lc_id = row.names(.)) %>%
    as_tibble %>%
    filter(!grepl("_5$|_6$", lc_id)) %>%
    rename(q2.5 = `0.025quant`,
           q50 = `0.5quant`, 
           q97.5 = `0.975quant`) %>%
    mutate(lc_id2 = as.integer(gsub(".*_(.*)", "\\1", lc_id)),
           var_name = factor(var_names[lc_id2], levels=var_names),
           var = c(temp1_pred/10, temp2_pred/100, precip1_pred, precip2_pred),
           mean = mean - m2$summary.fixed$mean[1], 
           q2.5 = q2.5 - m2$summary.fixed$mean[1], 
           q97.5 = q97.5 - m2$summary.fixed$mean[1]) #%>%
    # mutate_at(c("mean", "q2.5", "q97.5"), boot::inv.logit)

lincomb_summ_both <- bind_rows(lincomb_summ_1, lincomb_summ_2, .id="id") %>%
    mutate(model_type = c("non-spatial", "spatial")[as.integer(id)])
    
plot_effs <- ggplot(lincomb_summ_both, aes(x=var, y=mean, group=model_type)) + 
    geom_line(aes(col=model_type)) +
    geom_ribbon(aes(ymin=q2.5, ymax=q97.5), alpha=.2) + 
    facet_wrap(~var_name, scales="free_x") +
    theme(strip.text = element_text(face="bold", hjust=0), 
          strip.background = element_blank())

ggsave(paste0("figures/effPlot_", species_i_name, ".png"), plot=plot_effs, units="mm", 
       width=221, height=146) 

# plot map ----
lincomb_spatial <- m2$summary.lincomb.derived %>%
    mutate(lc_id = row.names(.)) %>%
    as_tibble %>%
    filter(grepl("_5$", lc_id)) %>%
    bind_cols(subset[c("x", "y")], .) %>%
    rename(q2.5 = `0.025quant`,
           q50 = `0.5quant`, 
           q97.5 = `0.975quant`) %>%
    mutate(#lc_id2 = as.integer(gsub(".*_(.*)", "\\1", lc_id)), 
           #var_name = factor(var_names[lc_id2], levels=var_names),
           #var = c(temp1_pred, temp2_pred, precip1_pred, precip2_pred), 
           mean = mean - m2$summary.fixed$mean[1], 
           q2.5 = q2.5 - m2$summary.fixed$mean[1], 
           q97.5 = q97.5 - m2$summary.fixed$mean[1], 
           overlap_0 = sign(q2.5) == sign(q97.5)) #%>%
    # mutate_at(c("mean", "q2.5", "q97.5"), boot::inv.logit)


# lincomb_spatial <- m2$summary.lincomb.derived %>%
#     mutate(lc_id = row.names(.)) %>%
#     as_tibble %>%
#     filter(grepl("_6$", lc_id)) %>%
#     bind_cols(subset[c("x", "y")], .) %>%
#     rename(q2.5 = `0.025quant`,
#            q50 = `0.5quant`, 
#            q97.5 = `0.975quant`) %>%
#     mutate(#lc_id2 = as.integer(gsub(".*_(.*)", "\\1", lc_id)), 
#         #var_name = factor(var_names[lc_id2], levels=var_names),
#         #var = c(temp1_pred, temp2_pred, precip1_pred, precip2_pred), 
#         mean = mean - m2$summary.fixed$mean[1], 
#         q2.5 = q2.5 - m2$summary.fixed$mean[1], 
#         q97.5 = q97.5 - m2$summary.fixed$mean[1], 
#         overlap_0 = sign(q2.5) == sign(q97.5)) #%>%
#     # mutate_at(c("mean", "q2.5", "q97.5"), boot::inv.logit)


# ggplot(lincomb_spatial) + 
#     geom_tile(aes(x, y, fill=(inv.logit(mean) - inv.logit(m2$summary.fitted.values$mean)))) + #, alpha=as.numeric(overlap_0))) + 
#     #scale_fill_viridis_c() +
#     scale_fill_gradient2(midpoint = 0, low="indianred", high="darkgreen") + 
#     scale_alpha_continuous(range = c(0.2, 1)) #+
#     geom_sf(data=st_rangemap_pixels, fill=NA, col="red") 

p1 <- ggplot(subset) + geom_tile(aes(x, y, fill=m1$summary.fitted.values$mean)) +
    geom_sf(data=st_rangemap_pixels, fill=NA, col="red") +
    scale_fill_viridis_b(limits=c(0, 1), breaks=seq(0, 1, .2)) +
    # guides(fill="none") +
    labs(x="", y="", fill="", title="(a) ignoring spatial autocorrelation") + 
    geom_sf(data=cnt, fill=NA)

p2 <- ggplot(subset) + geom_tile(aes(x, y, fill=m2$summary.fitted.values$mean)) +
    geom_sf(data=st_rangemap_pixels, fill=NA, col="red") +
    scale_fill_viridis_b(limits=c(0, 1), breaks=seq(0, 1, .2)) +
    # guides(fill="none") +
    labs(x="", y="", fill="", title="(b) climate + spatial autocorrelation")  + 
    geom_sf(data=cnt, fill=NA)

p3 <- ggplot(lincomb_spatial) + 
    geom_tile( aes(x, y, fill=mean, alpha=as.numeric(overlap_0))) +
    scale_alpha_continuous(range=c(0.3, 1)) +
    scale_fill_gradient2() +
    geom_sf(data=st_rangemap_pixels, fill=NA, col="red") +
    labs(x="", y="", fill="", title="(b) climate suitability (spatial)", 
         caption = paste0("dWAIC = ", round(m2$waic$waic - m2_null$waic$waic, 2))) +
    guides(alpha = "none", fill="none")  + 
    geom_sf(data=cnt, fill=NA)

p_all3 <- egg::ggarrange(p1, p2, p3, ncol=2)  
ggsave(paste0("figures/3maps_", species_i_name, ".png"), plot=p_all3, units="mm", width=290, height=210)
