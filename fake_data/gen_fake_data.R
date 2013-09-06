# r script to generate fake rasters, occurrence and evaluation data
rm(list = ls())

setwd('C:/Users/lina1864/Dropbox/Github/seegSDM/data')

library(geoR)
library(seegSDM)

genRaster <- function(mask = NULL, pars = c(1, 15), side = 200) {
  require(geoR)
  # generate a raster form a matern GRF
  f <- matrix(grf(side ^ 2, grid = 'reg', cov.pars = pars)$data, side)
  # make it a raster
  r <- raster(f)
  if (!is.null(mask)) {
    # set the same extent and projection as the mask
    extent(r) <- extent(mask)
    projection(r) <- projection(mask)
    # mask the resulting raster
    r <- mask(r, mask)
  }
  return (r)
}

# ~~~~~~~~~~~

# set the RNG seed for reproducible results
set.seed(321)

# generate a first GRF as a template
template <- genRaster()

# mask out the lowest 20% of values
template[template[] < quantile(template, 0.2)] <- NA

# set non-NA values to 0
template <- template * 0

# set the extent and projections
extent(template) <- c(-10, 10, -5, 15)
projection(template) <- wgs84(TRUE)

# generate covariate rasters
n_covs <- 3
# make n_covs more random rasters, masked by template
covariates <- brick(lapply(1:n_covs, function(i) genRaster(mask = template)))
# make the last one discrete
covariates[[n_covs]][] <- round(covariates[[n_covs]][] * 5)

# generate artificial admin units
admin0 <- aggregate(template, 100)
admin1 <- aggregate(template, 50)
admin2 <- aggregate(template, 25)
admin3 <- aggregate(template, 5)

# assign unique IDs at each level
admin0[] <- order(rnorm(ncell(admin0)))
admin1[] <- order(rnorm(ncell(admin1)))
admin2[] <- order(rnorm(ncell(admin2)))
admin3[] <- order(rnorm(ncell(admin3)))

# get back to the same resolution
admin0 <- disaggregate(admin0, 100)
admin1 <- disaggregate(admin1, 50)
admin2 <- disaggregate(admin2, 25)
admin3 <- disaggregate(admin3, 5)

# remask
admin0 <- mask(admin0, template)
admin1 <- mask(admin1, template)
admin2 <- mask(admin2, template)
admin3 <- mask(admin3, template)

# and combine into a brick
admin <- brick(admin0,
               admin1,
               admin2,
               admin3)


# fix their names
names(covariates) <- paste0('cov_', letters[1:n_covs])
names(admin) <- paste0('admin', 0:3)
names(template) <- 'template'

# save as RData objects in the data folder
save(covariates, file = 'covariates.RData')
save(admin, file = 'admin.RData')
save(template, file = 'template.RData')

# ~~~~~~~~~~~~~~~~
# simulating a fake species / disease

# create a *nasty* function giving probability of presence
# as a function of covariates
f <- function (covs) {
  # covs is a dataframe
  # z is a continuous latent variable
  # an intercept
  z <- 1 +
    # a polynomial effect of cov_a 
    1.5 * covs$cov_a - 15 * covs$cov_a ^ 2 +
    # a double sin curve effect of cov_b
    sin(covs$cov_b) + sin(covs$cov_b * 2) +
    # a threshold * log effect of cov_c
    0.1 * (ifelse(covs$cov_c > -1, 1, -1) * exp(covs$cov_c / 6)) +
    # some additional random noise (on top of likelihood
    # just to be extra mean)
    rnorm(nrow(covs))
  # use the probit link to return probability scale (flattened a little)
  return(pnorm(0.1 * z))
}

# generate a surface of probability of presence
# create a blank raster
prob <- template
# fill it with probabilities
prob[] <- f(as.data.frame(getValues(covariates)))

# create a fake evidence consensus layer from probabilities
# at admin1 level, 4th root to remove spike in one cell
consensus <- aggregate(prob ^ (1 / 4), 50)

# scale to (0, 1)
consensus <- consensus - min(consensus[], na.rm = TRUE)
consensus <- consensus / max(consensus[], na.rm = TRUE)

# then  to (-100, 100)
consensus <- consensus * 200 - 100

# disaggregate back to same resolution and remask
consensus <- disaggregate(consensus, 50)
consensus <- mask(consensus, template)

# save
save(consensus, file = 'consensus.RData')

n_eval <- 500

# sample randomly to get an evaluation set of n_eval
eval_pts <- bgSample(prob, n_eval, spatial = FALSE)
PA <- rbinom(n_eval, 1, extract(prob, eval_pts))
evaluation <- cbind(PA = PA, eval_pts)

# save this as an RData file
save(evaluation, file = 'evaluation.RData')

# generate an observation bias grid, sharpened on probability level
bias <- genRaster(mask = template)
bias[] <- pnorm(scale(bias[]) ^ 2)

# probability that the species would be observed
observation <- prob * bias

# number of occurrence points
n_occ <- 100

# sample them from the probability of observing
# (prob(observing|present) * prob(present))
occurrence <- bgSample(observation, n_occ, prob = TRUE, spatial = FALSE)

# add additional columns to mimic a real occurrence dataset
occurrence <- data.frame(UniqueID = 1:n_occ,
                         Admin = sort(sample(c(-9999, 3:0),
                                             n_occ,
                                             replace = TRUE,
                                             prob = c(0.5, 0.2, 0.15, 0.1, 0.05))),
                         occurrence)

occurrence$Admin <- as.integer(occurrence$Admin)
# add Area column, based on admin unit size

# calculate areas for different admin levels
area0 <- zonal(template + 0.01, admin0, fun = sum)
area1 <- zonal(template + 0.01, admin1, fun = sum)
area2 <- zonal(template + 0.01, admin2, fun = sum)
area3 <- zonal(template + 0.01, admin3, fun = sum)

# create rasters of the same
area0 <- reclassify(admin0, rcl = area0)
area1 <- reclassify(admin1, rcl = area1)
area2 <- reclassify(admin2, rcl = area2)
area3 <- reclassify(admin3, rcl = area3)

# add a column with NAs
occurrence$Area <- rep(NA, nrow(occurrence))

# fill in with the areas
occurrence[occurrence$Admin == 0, 'Area'] <-
  extract(area0, occurrence[occurrence$Admin == 0, c('x', 'y')])
occurrence[occurrence$Admin == 1, 'Area'] <-
  extract(area1, occurrence[occurrence$Admin == 1, c('x', 'y')])
occurrence[occurrence$Admin == 2, 'Area'] <-
  extract(area2, occurrence[occurrence$Admin == 2, c('x', 'y')])
occurrence[occurrence$Admin == 3, 'Area'] <-
  extract(area3, occurrence[occurrence$Admin == 3, c('x', 'y')])

# check they're all covered
stopifnot(all.equal(is.na(occurrence$Area), occurrence$Admin == -9999))

# save as an RData file
save(occurrence, file = 'occurrence.RData')