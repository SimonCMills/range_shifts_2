# demonstrate effectiveness of tensor spline over space. 
library(gstat); library(ggplot2); library(dplyr); library(mgcv); library(INLA)

# generate some spatially autocorrelated predictors
# Variogram model, with defined sill and range
varioMod <- vgm(psill=0.05, range=10, model='Exp')
# Set up an additional variable from simple kriging
zDummy <- gstat(formula=z ~ 1, locations = ~x+y, dummy=TRUE, 
                beta=1, model=varioMod, nmax=20)

# simulate data over grid
gridDim <- 60
xy <- expand.grid(x=1:gridDim, y=1:gridDim)
set.seed(3)
xyz <- predict(zDummy, newdata=xy, nsim=2)

# simulate some data from these predictions
df <- xyz %>%
    mutate(psi = boot::inv.logit(sim1*4 - 3.5), 
           det = rbinom(n(), 1, psi))

p1 <- ggplot(df, aes(x, y, fill=sim1)) + geom_tile() + coord_equal()
p2 <- ggplot(df, aes(x, y, fill=det)) + geom_tile() + coord_equal()
egg::ggarrange(p1, p2)


# mod 1 ----
# the model reclaims correct parameters for the sim1 variable
mod <- gam(det ~ sim1, data=df, family="binomial")
summary(mod)

# dataset with boundary ----
# any cell in the upper triangle has a probability of 0
df2 <- xyz %>%
    mutate(psi = ifelse(x+y < 60, boot::inv.logit(sim1*4 - 3.5), 0), 
           det = rbinom(n(), 1, psi))

p1 <- ggplot(df2, aes(x, y, fill=psi)) + geom_tile() + coord_equal()
p2 <- ggplot(df2, aes(x, y, fill=det)) + geom_tile() + coord_equal()
egg::ggarrange(p1, p2)

# mod 2 ----
# model that excludes spatial autocorrelation terms 
mod2 <- gam(det ~ sim1, data=df2, family="binomial")
summary(mod2)

# mod 3 ----
# model that includes spatial terms (and therefore is hopefully able to reclaim 
# correct parameters and posterior inference)
mod3 <- gam(det ~ sim1 + te(x, y), data=df2, family="binomial")
summary(mod3)

p1 <- ggplot(df2, aes(x, y, fill=psi)) + geom_tile() + coord_equal() +
    scale_fill_viridis_c(limits=c(0,1))
p2 <- ggplot(df2, aes(x, y, fill=boot::inv.logit(predict(mod3)))) + geom_tile() + coord_equal() +
    scale_fill_viridis_c(limits=c(0,1))
p3 <- ggplot(df2, aes(x, y, fill=boot::inv.logit(predict(mod2)))) + geom_tile() + coord_equal() +
    scale_fill_viridis_c(limits=c(0,1))
egg::ggarrange(p1, p2, p3)

# mod 4 ---
# INLA
df3 <- df2 %>%
    arrange(desc(y), x)
df_3_df <- df3 %>% mutate(idx = 1:n())
coordinates(df3) <- ~ x + y
gridded(df3) <- TRUE
r <- raster(df3)
r2 <- stars::st_as_stars(r)
poly_df3 <- st_as_sf(r2)

nbs <- poly2nb(poly_df3, queen=F)
adj <- nb2mat(nbs, style="B", zero.policy = T)
adj <- as(adj, "dgTMatrix")

mod4 <- inla(det ~ sim1 + f(idx, model="besag", graph=adj), 
           data=df_3_df, family="binomial", control.predictor = list(compute=T))

# CAR version (brms)
# set up distance and neighbourhood matrices
distance <- as.matrix(dist(df2[c("x", "y")]))
K <- nrow(distance)
W <- array(0, c(K, K))
W[distance == 1] <- 1 	

mod4 <- brm(det ~ sim1 + car(W, type = "icar"), 
           data = df2, data2 = list(W = W),
           family = bernoulli(), 
           chains = 4, cores = 4) 

summary(mod4)
pred4 <- fitted(mod4) %>%
    as_tibble

p1 <- ggplot(df2, aes(x, y, fill=psi)) + geom_tile() + coord_equal() +
    scale_fill_viridis_c(limits=c(0,1)) +
    labs(title = "Simulated P(Occupancy)", x="", y="") +
    guides(fill="none")
p2 <- ggplot(df2, aes(x, y, fill=boot::inv.logit(predict(mod3)))) + geom_tile() + coord_equal() +
    scale_fill_viridis_c(limits=c(0,1)) +
    labs(title = "Spatial GAM prediction", x="", y="") +
    guides(fill="none")
p3 <- ggplot(df2, aes(x, y, fill=boot::inv.logit(predict(mod2)))) + geom_tile() + coord_equal() +
    scale_fill_viridis_c(limits=c(0,1)) +
    labs(title = "Non-spatial GAM prediction", x="", y="") +
    guides(fill="none")

# p4 <- ggplot(df2, aes(x, y, fill=pred4$Estimate)) + geom_tile() + coord_equal() +
#     scale_fill_viridis_c(limits=c(0,1))

p4 <- ggplot(df_3_df, aes(x, y, fill = mod4$summary.fitted.values$mean)) + geom_tile() +
    coord_equal() + scale_fill_viridis_c(limits=c(0,1)) +
    labs(title = "CAR prediction (INLA)", x="", y="") +
    guides(fill="none")

p_all4 <- egg::ggarrange(p1, p3, p2, p4, ncol=2)
ggsave("figures/example_spatial_models.png", plot=p_all4)

# p1 <- ggplot(df2, aes(x, y, fill=boot::inv.logit(predict(mod3)) - psi)) + geom_tile() + coord_equal() +
#     scale_fill_viridis_c()
# p2 <- ggplot(df2, aes(x, y, fill=pred4$Estimate - psi)) + geom_tile() + coord_equal() +
#     scale_fill_viridis_c()
# egg::ggarrange(p1, p2)
# 
# ggplot(df2, aes(boot::inv.logit(predict(mod3)) - psi)) + 
#     geom_histogram(boundary=0, alpha=.5) +
#     geom_histogram(aes(pred4$Estimate - psi), boundary=0, 
#                    alpha=.5, fill="indianred")

## what does coverage look like?
# Note: currently only doing for spatial GAM
reps <- replicate(100, {
    df_sim <- predict(zDummy, newdata=xy, nsim=1) %>%
        mutate(psi = ifelse(x+y < 60, boot::inv.logit(sim1*4 - 3.5), 0), 
               det = rbinom(n(), 1, psi))
    
    mod1 <- gam(det ~ sim1, data=df_sim, family="binomial")
    mod2 <- gam(det ~ sim1 + te(x, y), data=df_sim, family="binomial")
    summ1 <- summary(mod1)
    summ2 <- summary(mod2)
    
    bind_rows(tibble(coef = summ1$p.coeff, p = summ1$p.pv, type="mod_naive"), 
              tibble(coef = summ2$p.coeff, p = summ2$p.pv, type="mod_autocorr"))
}, simplify = FALSE) %>% 
    bind_rows(., .id="id") %>%
    group_by(type, id) %>%
    mutate(term=c("intcpt", "slope")) 

df_sim <- predict(zDummy, newdata=xy, nsim=2) %>%
    mutate(psi = ifelse(x+y < 60, boot::inv.logit(sim2*4 - 3.5), 0), 
           det = rbinom(n(), 1, psi))

reps <- replicate(100, {
    # note using the second autocorrelated simulation for generating psi, but 
    # then modelling as a function of the first
    df_sim <- predict(zDummy, newdata=xy, nsim=2) %>%
        mutate(psi = ifelse(x+y < 60, boot::inv.logit(sim2*4 - 3.5), 0), 
               det = rbinom(n(), 1, psi))
    
    mod1 <- gam(det ~ sim1, data=df_sim, family="binomial")
    mod2 <- gam(det ~ sim1 + te(x, y), data=df_sim, family="binomial")
    summ1 <- summary(mod1)
    summ2 <- summary(mod2)
    
    bind_rows(tibble(coef = summ1$p.coeff, p = summ1$p.pv, type="mod_naive"), 
              tibble(coef = summ2$p.coeff, p = summ2$p.pv, type="mod_autocorr"))
}, simplify = FALSE) %>% 
    bind_rows(., .id="id") %>%
    group_by(type, id) %>%
    mutate(term=c("intcpt", "slope")) 

reps %>%
    filter(term == "slope") %>%
    ggplot(aes(coef, fill=type)) +
    geom_histogram(alpha=.5, position = "identity", boundary=0, binwidth=.1)

# There is still type 1 error in the model with spatial terms, but it is far 
# less pronounced than the naive model, which will identify significant effects
# of sim1 80% of the time, even at an alpha of 0.01!
# note: this is probably quite an extreme case because the spatial autocorrelation
# is a mixture of the gstat simulation spatial autocorrelation and the discrete
# boundary in the upper triangle. Certainly, the Type 1 error in Beale 2010 is 
# better than here
reps %>%
    filter(term == "slope") %>%
    group_by(type) %>%
    mutate(n = n()) %>%
    ggplot(aes(p, y = (..count../sum(..count..))*2, fill=type)) +
    # note: *2 is because sum(..count..) will divide by 200 rather than 100 
    geom_histogram(alpha=.5, position = "identity", boundary=0, binwidth=.05) +
    geom_hline(yintercept = .05, lty="longdash") + facet_wrap(~type)

reps %>% 
    filter(term=="slope") %>%
    group_by(type) %>%
    summarise(prop_05 = sum(p < 0.05)/n(), 
              prop_01 = sum(p < 0.01)/n())
