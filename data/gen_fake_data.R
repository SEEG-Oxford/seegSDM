# r script to generate fake rasters, occurrence and evaluation data
rm(list = ls())

setwd('C:/Users/lina1864/Dropbox/Github/seegSDM/data')

library(geoR)
library(seegSDM)

genRaster <- function(mask = NULL, pars = c(1, 15), side = 100) {
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
projection(template) <- wgs84()

n_covs <- 3
# make n_covs more random rasters, masked by template
rasters <- lapply(1:n_covs, function(i) genRaster(mask = template))

# save these as tiffs
lapply(1:n_covs, function(i) {
  writeRaster(rasters[[i]],
              filename = paste0('cov_', letters[i]),
              format = 'GTiff',
              overwrite = TRUE)
})

# save the template layer too
writeRaster(template,
            filename = 'template',
            format = 'GTiff',
            overwrite = TRUE)


# ~~~~~~~~~~~~~~~~
# simulating a fake species / disease

# create a *nasty* function giving probability of presence
# as a function of covariates
f <- function (covs) {
  # covs is a dataframe
  # z is a continuous latent variable
  # an intercept
  z <- -1 +
    # a polynomial effect of cov_a 
    1.5 * covs$cov_a - 15 * covs$cov_a ^ 2 +
    # a double sin curve effect of cov_b
    sin(covs$cov_b) + sin(covs$cov_b * 2) +
    # a threshold * log effect of cov_c
    ifelse(covs$cov_c > 0.2, 1, -1) * exp(2 * covs$cov_c) +
    # some additional random noise (on top of likelihood
    # just to be extra mean)
    rnorm(nrow(covs))
  # use the probit link to return probability scale (flattened a little)
  return(pnorm(0.3 * z))
}

# make the rasters a brick with correct names
b <- brick(rasters)
names(b) <- paste0('cov_', letters[1:n_covs])

# generate a surface of probability of presence
# create a blank raster
prob <- template
# fill it with probabilities
prob[] <- f(as.data.frame(getValues(b)))

# plot(prob)

n_eval <- 500


# sample randomly to get an evaluation set of n_eval
eval_pts <- bgSample(prob, n_eval, spatial = FALSE)
PA <- rbinom(n_eval, 1, extract(prob, eval_pts))
eval_pts <- cbind(PA = PA, eval_pts)

# write these out
write.csv(eval_pts,
          file = 'evaluation.csv',
          row.names= FALSE)

# plot these
# points(eval_pts[, 2:3], pch = ifelse(eval_pts[, 1], 16, 21))

# generate an observation bias grid, sharpened on probability level
bias <- genRaster(mask = template)
bias[] <- pnorm(scale(bias[]) ^ 2)

# probability that the species would be observed
observation <- prob * bias

# number of occurrence points
n_occ <- 100
occ_pts <- bgSample(observation, n_occ, prob = TRUE, spatial = FALSE)

write.csv(occ_pts,
          file = 'occurrence.csv',
          row.names= FALSE)

