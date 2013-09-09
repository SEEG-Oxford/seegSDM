# function file for seegSDM package

# checking inputs
checkOccurrence <- function(occurrence, consensus, admin, consensus_threshold = -25,
                             area_threshold = 1, max_distance = 0.05, spatial = TRUE,
                             verbose = TRUE) {
  
  # if occurrence
  
  
  # check that the consensus layer is projected
  if (projection(consensus) != wgs84(TRUE)@projargs) {
    stop ('consensus is not projected, please project it (or correct the coordinate system)')
  }
  
  # expected column names and data types
  expected_names <- c('UniqueID',
                      'Admin',
                      'Year',
                      'x',
                      'y',
                      'Area')
  
  expected_classes <- c('integer',
                        'integer',
                        'integer',
                        'numeric',
                        'numeric',
                        'numeric')
  
  # ~~~~~~~~~~~~~~
  # check if any expected column names are missing
  missing <- expected_names[!expected_names %in% names(occurrence)]
  
  # if they are, throw an error 
  if (length(missing) > 0) {
    stop(paste('missing columns:', paste(missing, collapse = ', ')))
  }
  
  # ~~~~~~~~~~~~~~
  # check column data types
  
  # get data types of columns in occurrence
  occurrence_classes <- sapply(occurrence, class)
  
  # get the expected data types in the same order
  order <- match(expected_names, names(occurrence))
  expected_classes_ordered <- expected_classes[order]
  
  # do they match up?
  classes_match <- occurrence_classes == expected_classes_ordered
  
  # if not, stop with a helpful error
  if (!all(classes_match)) {

    # list problems
    message_vector <- sapply(which(!classes_match),
                             function(i) paste(names(occurrence)[i],
                                               'is',
                                               occurrence_classes[i],
                                               'but should be',
                                               expected_classes_ordered[i]))
    
    stop(paste("data types don't match\n",
               paste(message_vector, collapse = '\n')))
  }
  
  
  # ~~~~~~~~~~~~~~
  # check that Admin contains a reasonable number (either admin level 0:3, or -9999 for points)
  bad_admin <- !(occurrence$Admin %in% c(0:3, -9999))
  
  if (any(bad_admin)) {
    stop (paste('bad Admin codes for records with these UniqueIDs:',
                 occurrence$UniqueID[bad_admin]))
  }  
  
  
  # ~~~~~~~~~~~~~~
  # check that there are no duplicated polygon/year combinations
  
  # if there are any polygons
  if (any(occurrence$Admin != -9999)) {
    
    # find and add GAUL codes, maintaining occurrence as a dataframe
    occurrence$GAUL <- getGAUL(occurrence2SPDF(occurrence), admin)$GAUL
    
    # get an index for which records are polygons
    poly_idx <- occurrence$Admin != -9999
    
    # check for duplicates (Admin, GAUL and Year all the same)
    poly_dup <- duplicated(occurrence[poly_idx, c('Admin', 'GAUL', 'Year')])
  
    # if any are duplicated give a warning but proceed
    if (any(poly_dup)) {
      stop (paste('there are duplicated polygon / year combinations
                  at records with UniqueId:',
                     paste(occurrence$UniqueId[poly_idx][poly_dup], collapse = ', ')))
    }
  }
  
  # ~~~~~~~~~~~~~~
  # remove polygons over area limit
  big_polygons <- occurrence$Admin != -9999 &
    occurrence$Area > area_threshold
  
  # notify the user
  if (verbose) {
    cat(paste(sum(big_polygons),
              "polygons had areas greater than the area threshold of",
              area_threshold,
              "and will be removed.\n"))
  }
  
  # remove these from occurrence
  occurrence <- occurrence[!big_polygons, , drop = FALSE]
  
  
  # ~~~~~~~~~~~~~~
  # see if any coordinates have values below (or equal to) the evidence consensus threshold
  consensus_scores <- extract(consensus, occurrence[, c('x', 'y')])
  low_consensus <- consensus_scores <= consensus_threshold 
  
  # if there were any
  if (any(low_consensus)) {
    
    # let the user know
    if (verbose) {
      cat(paste('removing',
                sum(low_consensus),
                "points which were in areas with evidence consensus value below the threshold of",
                consensus_threshold,
                '\n'))
    }
    
    # then remove them
    occurrence <- occurrence[!low_consensus, , drop = FALSE]
  }
  
  # ~~~~~~~~~~~~~~
  # find points (and centroids) outside mask
  vals <- extract(consensus, occurrence[, c('x', 'y')], drop = FALSE)
  outside_mask <- is.na(vals)
  
  if (any(outside_mask)) {
    # if there are any outside the mask...
    
    # notify the user
    if (verbose) {
      cat(paste())
    }
    
    # try to find new coordinates for these
    new_coords <- nearestLand(xyFromCell(which(outside_mask)),
                              consensus,
                              max_distance)

    # replace those coordinates in occurrence
    occurrence[outside_mask, c('x', 'y')] <- new_coords
    
    # how many of these are still not on land
    still_out <- is.na(new_coords[, 1])
    
    # tell the user
    if (verbose) {
      cat(paste(sum(outside_mask),
                'points were outside consensus.\n',
                sum(!still_out),
                'points were moved to dry land, ',
                sum(still_out),
                'points were further than',
                max_distance,
                'decimal degrees out and have been removed'))
    }
    
    # update outside_mask
    outside_mask[outside_mask] <- still_out

    # and remove still-missigng points from occurrence and vals
    occurrence <- occurrence[!outside_mask, , drop = FALSE]
    values <- values[!outside_mask]
  }
  
  # ~~~~~~~~~~~~~~
  # if spatial = TRUE, return an SPDF, otherwise the dataframe
  if (spatial) {
    occurrence <- SpatialPointsDataFrame(cbind(occurrence$x, occurrence$y),
                                          occurrence,
                                          proj4string = wgs84(TRUE))
  }
  
  # return corrected occurrence dataframe/SPDF
  return (occurrence)
}

checkRasters <- function (rasters, template, cellbycell = FALSE) {
  
  # check whether the raster* object 'rasters' conforms to the 'template' rasterLayer
  # By default the extent and projection are compared.
  # If 'cellbycell = TRUE' the pixel values are individually compared with the template
  # though this can be very slow & RAM-hungry with large rasters.
  # the function throws an error if anything is wrong and returns 'rasters' otherwise.
  
  if (extent(rasters) != extent(template)) {
    stop('extents do not match, see ?extent')
  }
  
  if (projection(rasters) != projection(template)) {
    stop('projections do not match, see ?projection')
  }
  
  if (ncell(rasters) != ncell(template)) {
    stop('number of cells do not match, see ?ncell')
  }
  
  if (cellbycell) {
    # get the template pixels and find the nas
    template_pixels <- getValues(template[[1]])
    template_nas <- is.na(template_pixels)
    
    # same for rasters
    rasters_pixels <- getValues(rasters)
    rasters_nas <- is.na(rasters_pixels)
    
    # how many layers in rasters
    n <- nlayers(rasters)
    
    # if rasters is multi-band
    if (n > 1) {
      
      # create an empty vector to store results
      pixel_mismatch <- vector(length = n)
      
      # loop through, checking NAs are in the same place
      for (i in 1:n) {
        pixel_mismatch[i] <- any(rasters_nas[, i, drop = TRUE] != template_nas)
      }
      
      if (any(pixel_mismatch)) {
        stop(paste0('mismatches between layers ',
                   names(rasters)[pixel_mismatch],
                   ', see ?resample'))
      }
      
    }
    
    if (n == 1 & any(rasters_nas != template_nas)) {
      stop('mismatch between layers, see ?resample')
    }
    
  }
  
  return (rasters)
}

nearestLand <- function (points, raster, max_distance) {
  # get nearest non_na cells (within a maximum distance) to a set of points
  # points can be anything extract accepts as the y argument
  # max_distance is in the map units if raster is projected
  # or metres otherwise
  
  # function to find nearest of a set of neighbours or return NA
  nearest <- function (lis, raster) {
    neighbours <- matrix(lis[[1]], ncol = 2)
    point <- lis[[2]]
    # neighbours is a two column matrix giving cell numbers and values
    land <- !is.na(neighbours[, 2])
    if (!any(land)) {
      # if there is no land, give up and return NA
      return (c(NA, NA))
    } else{
      # otherwise get the land cell coordinates
      coords <- xyFromCell(raster, neighbours[land, 1])
      
      if (nrow(coords) == 1) {
        # if there's only one, return it
        return (coords[1, ])
      }
      
      # otherwise calculate distances
      dists <- sqrt((coords[, 1] - point[1]) ^ 2 +
                      (coords[, 2] - point[2]) ^ 2)
      
      # and return the coordinates of the closest
      return (coords[which.min(dists), ])
    }
  }
  
  # extract cell values within max_distance of the points
  neighbour_list <- extract(raster, points,
                            buffer = max_distance,
                            cellnumbers = TRUE)
  
  # add the original point in there too
  neighbour_list <- lapply(1:nrow(points),
                           function(i) {
                             list(neighbours = neighbour_list[[i]],
                                  point = points[i, ])
                           })
  
  return (t(sapply(neighbour_list, nearest, raster)))
}

occurrence2SPDF <- function (occurrence) {
  # helper function to convert an occurrence dataframe
  # i.e. one which passes checkOccurrence into a SpatialPointsDataFrame object
  
  # get column numbers for coordinates
  coord_cols <- which(names(occurrence) %in% c('x', 'y'))
  
  # convert to SPDF
  occurrence <- SpatialPointsDataFrame(occurrence[, coord_cols],
                                occurrence,
                                coords.nrs  = coord_cols,
                                proj4string = wgs84(TRUE))
  return (occurrence)
}

getGAUL <- function (occurrence, admin) {
  # given a SpatialPointsDataFrame (occurrence) and a brick of GAUL codes
  # for admin levels 0 to 3, extract codes for polygon records
  
  # check that the coordinate references match
  if (projection(occurrence) != projection(admin)) {
    stop ('projection arguments for occurrence and admin do not match')
  }
  
  # start with GAUL column filled with NAs
  occurrence$GAUL <- rep(NA, nrow(occurrence))
  
  # loop through the 4 admin levels
  for (level in 0:3) {
    
    # find relevant records
    idx <- occurrence$Admin == level
    
    # if there are any polygons at that level
    if (any(idx)) {
      
      # extract codes into GAUL column
      occurrence$GAUL[idx] <- extract(admin[[level + 1]], occurrence[idx, ])
    }
  }
  
  return (occurrence)
}


extractAdmin <- function (occurrence, covariates, admin, fun = mean) {
  
  # get indices for polygons in occurrence
  poly_idx <- which(occurrence$Admin != -9999)
  
  # create an empty matrix to store the results
  ex <- matrix(NA,
               nrow = length(poly_idx),
               ncol = nlayers(covariates))
  
  # give them useful names
  colnames(ex) <- names(covariates)
  
  # pull out relevant vectors for these
  all_admins <- occurrence$Admin[poly_idx]
  all_GAULs <- occurrence$GAUL[poly_idx]
  
  for (level in 0:3) {
    
    # are each of the *polygon records* of this level?
    level_idx <- all_admins == level
    
    # if there are any polygons at that level
    if (any(level_idx)) {
      
      # get GAUL_codes levels of interest
      level_GAULs <- all_GAULs[level_idx]
      
      # clone the admin layer we want
      ad <- admin[[level + 1]]
      
      # get all (and possibly some excess) unique GAUL codes
      level_all_codes <- ad@data@min : ad@data@max
      
      # remove the ones we *do* want from this list, by index
      level_all_codes <- level_all_codes[-(level_GAULs - ad@data@min + 1)]
      
      # then reclassify ad to mask out the ones we don't want
      # this should speed up the zonal operation
      ad <- reclassify(ad,
                       cbind(level_all_codes,
                             rep(NA, length(level_all_codes))))
      
      # extract values for each zone, aggregating by 'fun'
      zones <- zonal(covariates, ad, fun = fun)
      
      # match them to polygons (accounts for change in order and for duplicates)
      which_zones <- match(level_GAULs, zones[, 1])
      
      # add them to the results matrix
      ex[level_idx, ] <- zones[which_zones, -1]
      
    }
    
  }
  
  return (ex)
  
}


runBRT <- function (data, gbm.x, gbm.y, pred.raster,
                    wt.fun = function(PA) rep(1, length(PA)),
                    max_tries = 5, verbose = FALSE,
                    tree.complexity = 4, learning.rate = 0.005,
                    bag.fraction = 0.75, n.trees = 10,
                    n.folds = 10, max.trees = 40000,
                    step.size = 10, ...)
  # wrapper to run a BRT model with Sam's defaults & a default 50:50 PA weighting
  # and return covariate effects, relative influence,
  # and a prediction map (on the probability scale).
  # background points are weighted at 4 * presence points,
  # mimicking prevalence of 0.2
{
  # keep running until model isn't null (can happen with wrong leaning rate etc.)

  # set up for the while loop
  m <- NULL
  tries <- 0
  
  # try 5 times
  while (is.null(m) & tries < max_tries) {
    # fit the model, if it fails m will be NULL and the loop will continue
    m <- gbm.step(data = data,
                  gbm.x = gbm.x,
                  gbm.y = gbm.y,
                  step.size = 10,
                  tree.complexity = tree.complexity,
                  verbose = verbose,
                  learning.rate = learning.rate,
                  bag.fraction = bag.fraction,
                  n.trees = n.trees,
                  max.trees = max.trees,
                  plot.main = FALSE,
                  site.weights = wt.fun(data[, gbm.y]),
                  ...)
    
    # add one to the tries counter
    tries <- tries + 1
  }
  
  # throw an error if it still hasn't worked after max_tries
  if (tries >= max_tries) {
    stop (paste('Unexpected item in the bagging area!\nStopped after',
                max_tries,
                'attempts, try reducing the step.size or learning rate'))
  }
    
  
  # otherwise return the list of model objects
  # (predict grids, relative influence stats and prediction map)
  list(model = m,
       effects = lapply(1:length(gbm.x), function(i) plot(m, i, return.grid = TRUE)),
       relinf = summary(m, plotit = FALSE),
       pred = predict(pred.raster, m, type = 'response', n.trees = m$n.trees))
}

getRelInf <- function (models, plot = FALSE, quantiles = c(0.025, 0.975), ...)
  # given a list of BRT model bootstraps (each an output from runBRT)
  # get the mean and quantiles (95% by default) of relative influence
  # for covariates. Optionally plot a boxplot of the relative influences.
  # The dots argument is passed to 'boxplot'.
{
  rel.infs <- t(sapply(models, function(x) x$relinf[, 2]))
  colnames(rel.infs) <- models[[1]]$relinf[, 1]
  if (plot) {
    boxplot(rel.infs, ylab = 'relative influence', col = 'light grey', ...)
  }
  relinf <- cbind(mean = colMeans(rel.infs),
                  t(apply(rel.infs, 2, quantile, quantiles)))
  return (relinf)
}

getEffectPlots <- function (models, plot = FALSE, quantiles = c(0.025, 0.975), ...) {
  # given a list of BRT model bootstraps (each an output from runBRT)
  # get the mean and quantiles (95% by default) lines for effect plots
  # returns a list of matrices for each covariate and optionally plots the results
  # dots argument allows some customisation of the plotting outputs
  
  getLevels <- function (cov, models) {
    # get all the possible levels of covariate 'cov'
    # from the BRT ensemble 'models'
    if (models[[1]]$model$var.type[cov] != 0) {
      # if gbm thinks it's a factor, calculate all possible levels
      levels <- lapply(models, function (m) m$effects[[cov]][, 1])
      levels <- unique(unlist(levels))
      return (levels)
    } else {
      # if it isn't a factor return NULL
      return (NULL)
    }
  }
  
  matchLevels <- function (m, cov, levels) {
    # get effects of covariate for this model
    eff <- m$effects[[cov]]
    # blank vector of response for ALL levels
    y <- rep(NA, length(levels))
    # match those levels present here
    y[match(eff[, 1], levels)] <- eff[, 2]
    # return the marginal prediction
    return (y)
  }
  
  getEffect <- function (cov, models) {
    # get levels
    levels <- getLevels(cov, models)
    if (is.null(levels)) {
      # if it's continuous
      # get x
      x <- models[[1]]$effects[[cov]][, 1]
      # get matrix of y
      y <- sapply(models, function (x) x$effects[[cov]]$y)
    } else {
      # it it's a factor
      x <- levels
      y <- sapply(models, matchLevels, cov, levels)
    }
    
    # get means and quantiles
    mn <- rowMeans(y, na.rm = TRUE)
    qs <- t(apply(y, 1, quantile, quantiles, na.rm = TRUE))
    # and return these
    return (cbind(covariate = x,
                  mean = mn,
                  qs))
  }
  
  ncovs <- length(models[[1]]$effects)
  
  names <- sapply(models[[1]]$effects, function (x) names(x)[1])
  
  effects <- lapply(1:ncovs, getEffect, models)
  
  names(effects) <- names
  
  if (plot) {
#     op <- par(mfrow = n2mfrow(ncovs), mar = c(5, 4, 4, 2) + 0.1)
    for (i in 1:ncovs) {
      eff <- effects[[i]]
      
      if (is.null(getLevels(i, models))) {
        # if it's a continuous covariate do a line and CI region
        
        # blank plot
        plot(eff[, 2] ~ eff[, 1], type = 'n', ylim = range(eff[, 3:4]),
             xlab = names[i], ylab = paste('f(', names[i], ')', sep = ''), ...)
        # uncertainty region
        polygon(c(eff[, 1], rev(eff[, 1])),
                c(eff[, 3], rev(eff[, 4])),
                col = 'light grey', lty = 0)
        # mean line
        lines(eff[, 2] ~ eff[, 1], lwd = 2)
        
      } else {
        # if it's a factor, do dots and bars
        
        # blank plot
        plot(eff[, 2] ~ eff[, 1], type = 'n', ylim = range(eff[, 3:4]),
             xlab = names[i], ylab = paste('f(', names[i], ')', sep = ''), ...)
        # uncertainty lines
        segments(eff[, 1], eff[, 3], y1 = eff[, 4],
                 col = 'light grey', lwd = 3)
        # points
        points(eff[, 2] ~ eff[, 1], pch = 16, cex = 1.5)
      }
      
    }
#     par(op)
  }
  
  return(effects)
}

combinePreds <- function (preds, quantiles = c(0.025, 0.975),
                          parallel = FALSE, maxn = NULL)
  # function to calculate means and quantiles for each cell, given a rasterBrick
  # or rasterStack where each layer is a single ensemble prediction.
  # If a snowfall cluster is running, 'parallel = TRUE' sets the function to run
  # in parallel. 'maxn' gives the maximum number of cells in each batch. If NULL
  # this is set to fill all the cpus once.
{
  # function to get mean and quantiles
  getStats <- function (x, quants) {
    if (nrow(x) == 1) c(mean(x), quantile(x, quants, na.rm = TRUE))
    else cbind(mean = rowMeans(x),
               t(apply(x, 1, quantile, quants, na.rm = TRUE)))
  }
  
  stopifnot(nlayers(preds) > 1)
  
  # extract data to RAM and find non-NA values
  alldata <- getValues(preds)
  nas <- is.na(alldata[, 1])
  dat <- alldata[!nas, ]
  
  if (parallel) {
    
    n <- nrow(dat)
    cpu <- sfCpus()
    stopifnot(cpu > 1)
    
    if (is.null(maxn)) {
      maxn <- n %/% cpu + ifelse(n %% cpu > 0, 1, 0)
    }
    
    # split up dat
    idxs <- splitIdx(n, maxn)
    # get stats in parallel
    stats <- sfLapply(idxs,
                      function(idx, dat, quants) getStats(dat[idx[1]:idx[2], ], quants),
                      dat, quantiles)
    
    stats <- do.call(rbind, stats)
    
  } else {
    
    stats <- getStats(dat, quants = quantiles)
    
  }
  
  if (nlayers(preds) < 3) ans <- brick(preds[[1]], preds[[1]], preds[[1]])
  else ans <- preds[[1:3]]
  ans[!nas, ] <- stats
  names(ans) <- c('mean', 'lower', 'upper')  
  ans
}

rmse <- function(truth, prediction)
  # root mean squared error of prediction from true probability
{
  sqrt(mean(abs(prediction - truth) ^ 2))
}

devBern <- function (truth, prediction)
  # predictive deviance from a bernoulli distribution
{
  -2 * sum(dbinom(truth, 1, prediction, log = TRUE))
}

subsample <- function (data, n, minimum = c(5, 5), prescol = 1, replace = FALSE) {
  # get a random subset of 'n' records from 'data', ensuring that there
  # are at least 'minimum[1]' presences and 'minimum[2]' absences.
  # assumes by default that presence/absence code is in column 1 ('prescol')
  OK <- FALSE
  # until criteria are met
  while (!OK) {
    # take a subsample
    sub <- data[sample(1:nrow(data), n, replace = replace), ]
    # count presences and absences
    npres <- sum(sub[, prescol]) 
    nabs <- sum(1 - sub[, prescol])
    # if they are greater than or equal to the minimum, criteria are met
    if (npres >= minimum[1] & nabs >= minimum[2]) {
      OK <- TRUE
    }
  }
  return (sub)
}

bgSample <- function(raster, n = 1000, prob = FALSE, replace = FALSE, spatial = TRUE)
  # sample N random pixels (not NA) from raster. If 'prob = TRUE' raster is assumed
  # to be a bias grid and points sampled accordingly. If 'sp = TRUE' a
  # SpatialPoints* object is returned, else a matrix of coordinates
{
  pixels <- which(!is.na(raster[]))
  if (prob) {
    prob <- raster[pixels]
  } else {
    prob <- NULL
  }
  points <- sample(x = pixels, size = n, replace = replace, prob = prob)
  points <- xyFromCell(raster, points)
  if (spatial) {
    return (SpatialPoints(points, proj4string = raster@crs))
  } else {
    return (points)
  }
}

bgDistance <- function (n, points, raster, distance, ...) {
  # sample (nbg) background points from within a region within a given distance
  # (distance) of occurrence points (an sp object given by sp). The dots argument
  # can be used to pass options to bgSample. This is waaaay more efficient than
  # anything using gBuffer.
  
  r <- rasterize(points@coords, raster)
  buff <- buffer(r, width = distance) * raster
  bgSample(buff, n = n, ...)
}

extractBhatt <- function (pars,
                          occurrence,
                          covariates,
                          consensus,
                          admin,
                          threshold = -25,
                          factor = rep(FALSE, nlayers(covariates)),
                          return_points = FALSE,
                          ...) {
  # generate and extract pseudo-absence/pseudo-presence and occurrence
  # covariates for a single BRT run using the procedure in Bhatt paper
  # na = number of pseudo-absences per occurrence
  # np = number of pseudo-presences per occurrence
  # mu = distance (in decimal degrees) from which to select these
  # a dataframe is created with the presence/absence code in the first column
  # and extracted values in the other columns. factor can be used to
  # coerce covariates into factors in the dataframe.
  # dots is passed to bgDistance
  
  # coerce pars into a numeric vector (lapply can pass it as a dataframe)
  pars <- as.numeric(pars)
  
  
  # make sure occurrence is the correct type
  if (class(occurrence) != "SpatialPointsDataFrame") {
    stop ('occurrence must be a SpatialPointsDataFrame object')
  }
  
  # make sure the template raster and occurrence are projected
  if (projection(consensus) != wgs84(TRUE)@projargs) {
    stop ('consensus must be projected WGS84, try ?projectExtent and ?wgs84.')
  }
  if (projection(covariates) != wgs84(TRUE)@projargs) {
    stop ('covariates must be projected WGS84, try ?projectExtent and ?wgs84.')
  }
  if (projection(admin) != wgs84(TRUE)@projargs) {
    stop ('admin must be projected WGS84, try ?projectExtent and ?wgs84.')
  }
  if (!is.projected(occurrence)) {
    stop ('occurrence must be projected, try ?spTransform and ?wgs84.')
  }
  
  # get number of occurrence records
  no <- length(occurrence)
  
  # calculate
  na <- ceiling(pars[1] * no)
  np <- ceiling(pars[2] * no)
  mu <- pars[3]
  
  # if any are required generate pseudo-absences and pseudo presences,
  # extract and add labels
  
  # pseudo-absences
  if (na > 0) {

    # modify the consensus layer (-100:100) to probability scale (0, 1)
    abs_consensus <- 1 - (consensus + 100) / 200
    
    # sample from it, weighted by consensus (more likely in -100, impossible in +100)
    p_abs <- bgDistance(na, occurrence, abs_consensus, mu, prob = TRUE, ...)
    p_abs_covs <- extract(covariates, p_abs)
    p_abs_data <- cbind(PA = rep(0, nrow(p_abs_covs)), p_abs_covs)
    
  } else {
    
    p_abs <- p_abs_data <- NULL

  }
  
  # pseudo-presences
  if (np > 0) {

    if (!exists('abs_consensus')) {
      
      # if a pseudo-absence consensus layer already exists then save some computation
      pres_consensus <- 1 - abs_consensus
      
    } else {
      
      # otherwise create it from consensus
      pres_consensus <- (consensus + 100) / 200
      
    }
    
    # threshold it (set anything below threshold to 0)
    threshold <- (threshold + 100) / 200
    pres_consensus[pres_consensus <= threshold] <- 0
    
    # sample from it, weighted by consensus (more likely in 100, impossible
    # below threshold)
    p_pres <- bgDistance(np, occurrence, pres_consensus, mu, prob = TRUE, ...)
    p_pres_covs <- extract(covariates, p_pres)
    p_pres_data <- cbind(PA = rep(1, nrow(p_pres_covs)), p_pres_covs)
    
  } else {
    
    p_pres <- p_pres_data <- NULL
    
  }
  
  # extract covariates for occurrence and add presence label

  # create an empty matrix for th occurrence covariate records
  occ_covs <- matrix(NA,
                     nrow = nrow(occurrence),
                     ncol = nlayers(covariates))
  
  # give them their names
  colnames(occ_covs) <- names(covariates)
  
  # find points
  points <- occurrence$Admin == -9999

  # if there are points
  if (any(points)) {
    # extract and add to the results
    occ_covs[points, ] <- extract(covariates, occurrence[points, ])
  }
  
  # if there are any polygons
  if (any(!points)) {
    # extract them, but treat factors and non-factors differently

    if (any(factor)) {
      # if there are any factors, get mode of polygon
      occ_covs[!points, factor] <- extractAdmin(occurrence,
                                                 covariates[[which(factor)]],
                                                 admin,
                                                 fun = modal)
    }
    if (any(!factor)) {
      # if there are any continuous, get mean of polygon
      occ_covs[!points, !factor] <- extractAdmin(occurrence,
                                                 covariates[[which(!factor)]],
                                                 admin,
                                                 fun = mean)
      
    }
  }
  
  # add a vector of ones
  occ_data <- cbind(PA = rep(1, nrow(occ_covs)), occ_covs)
  
  # combine the different datasets and convert to a dataframe
  all_data <- rbind(occ_data, p_abs_data, p_pres_data)
  all_data <- as.data.frame(all_data)
  
  # if there are any factor covariates, convert the columns in the dataframe
  if (any(factor)) {
    facts <- which(factor)
    for (i in facts) {
      # plus one to avoid the PA column
      all_data[, i + 1] <- factor(all_data[, i + 1])
    }
  }

  # return list with the dataframe and possibly the SpatialPoints objects
  
  if (return_points) {
    # if the points are required, return a list
    return (list(data = all_data,
                 pseudo_absence = p_abs,
                 pseudo_presence = p_pres))
  } else {
    # otherwise just the dataframe
    return (all_data)
  }
}


biasGrid <- function(polygons, raster, sigma = 30)
  # create a bias grid from polygons using a gaussian moving window smoother
  # sigma is given in the units of the projection
{
  ras <- raster
  ras[] <- 0
  # weighted occurrence raster (areas for polys roughly sum to 1)
  areas <- sapply(polygons@polygons, function(x) {
    x@Polygons[[1]]@area
  }) / prod(res(raster))
  ras <- rasterize(polygons, ras, field = 1/areas,
                   fun = sum, update = TRUE)
  # get sigma in pixels
  sigma <- ceiling(sigma / mean(res(raster)))
  # use a window of 5 * sigma
  w <- gaussWindow(sigma * 5, sigma)
  print(system.time(ras <- focal(ras, w = w / sum(w), pad = TRUE, na.rm = TRUE)))
  # replace mask
  ras[is.na(ras[])] <- 0
  mask(ras, raster)
}

# bufferMask <- function(feature, raster, km = 100, val = NULL, maxn = 1000, parallel = FALSE, dissolve = FALSE)
#   # Returns a mask giving the non-NA areas of 'raster' which are within
#   # 'km' kilometres of 'feature'. If val is specified, non-NA values are
#   # assigned val, if null the values in (the first layer of) raster are
#   # used. gBuffer does poorly with very complex feature, so buffering is
#   # carried out in batches of size 'maxn'. If a snowfall cluster is
#   # running, setting 'parallel = TRUE', runs the buffering in parallel. 
#   # Uses rgeos::gBuffer and rgdal::spTransform
# 
#   # NOTE: bgDistance is MUCH FASTER & more memory safe for generating pseuco-absences
#   
#   {
#   # if raster is a stack or brick, get the first layer
#   if (nlayers(raster) > 1) raster <- raster[[1]]
#   
#   # convert points to projected latlong
#   if (proj4string(feature) != wgs84(TRUE)@projargs) {
#     feature <- spTransform(feature, CRSobj = wgs84(TRUE))
#   }
#   
#   # buffer by 'km' kilometres, in batches of maxn since gBuffer is slow
#   # for complex features
#   split <- splitIdx(length(feature), maxn)
#   
#   # if a snowfall cluster is running
#   if (parallel) {
#     
#     # export necessary things
#     sfLibrary(rgeos)
#     #    sfExport('feature', 'km')
#     
#     # run in parallel
#     buffers <- sfLapply(split, function(idx) gBuffer(feature[idx[1]:idx[2], ], width = km * 1000))
#     #   buffers <- sfLapply(split, function(idx) gBuffer(feature[idx[1]:idx[2]], width = km * 1000))
#     
#   } else {
#     
#     # otherwise, run sequentially
#     buffers <- lapply(split, function(idx) gBuffer(feature[idx[1]:idx[2], ], width = km * 1000))
#     #   buffers <- lapply(split, function(idx) gBuffer(feature[idx[1]:idx[2]], width = km * 1000))
#     
#   }
#   
#   # recombine these with gUnion if necessary
#   #  if (dissolve) {
#   if (length(buffers) > 1) {
#     buffer <- buffers[[1]]
#     for (i in 2:length(buffers)) {
#       buffer <- gUnion(buffer, buffers[[i]])
#     }
#   } else {
#     buffer <- buffers[[1]]
#   }
#   # }
#   
#   # reproject to CRS of rasterlayer
#   buffer <- spTransform(buffer, CRS(proj4string(raster)))
#   
#   # mask the raster layer (waaay faster than the mask function)
#   buffmask <- raster * rasterize(buffer, raster)
#   
#   # if a val is given, overwrite values
#   if (!is.null(val)) {
#     buffmask[which(!is.na(buffmask[]))] <- val
#   }
#   
#   buffmask
# }

buildSP <- function (files)
  # given a character vector of ESRI shapefiles (.shp)
  # build a SpatialPolygons object combining them.
{
  shps <- lapply(files, shapefile)
  polylis <- lapply(shps, function(x) x@polygons[[1]])
  for(i in 1:length(polylis)) polylis[[i]]@ID <- paste(i)
  SpatialPolygons(polylis, proj4string = shps[[1]]@proj4string)
}

# centre2poly <- function(centre, res)
#   # given pixels centres and raster resolution,
#   # return SpatialPolygons for the pixels
# {
#   ldru <- rep(centre, 2) + c(-res, res)/2
#   Polygon(rbind(ldru[1:2], # ld
#                 ldru[c(1, 4)], # lu
#                 ldru[c(3, 4)], # ru
#                 ldru[c(3, 2)], # rd
#                 ldru[1:2])) # ld
# }
# 
# pixels2polys <- function(raster)
#   # turn all non-na pixels into polys - don't run on a big file!
# {
#   centre <- xyFromCell(raster, which(!is.na(raster[])))
#   polys <- apply(centre, 1, centre2poly, res(raster))
#   polys <- lapply(polys, function(x) Polygons(list(x), 1))
#   for(i in 1:length(polys)) polys[[i]]@ID <- as.character(i)
#   SpatialPolygons(polys, proj4string = raster@crs)
# }
# 


en2os <- function (en)
{
  # convert 2-column matrix of eastings/northings to alphanumeric OX grid references
  lookup <- data.frame(N = c('H','H','H','H','H','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','T','T','T','T','T','T','T','T'),
                       E = c('P','T','U','Y','Z','A','B','C','D','F','G','H','J','K','L','M','N','O','R','S','T','U','W','X','Y','Z','C','D','E','H','J','K','M','N','O','P','R','S','T','U','V','W','X','Y','Z','A','F','G','L','M','Q','R','V'),
                       X = c(4,3,4,3,4,0,1,2,3,0,1,2,3,4,0,1,2,3,1,2,3,4,1,2,3,4,2,3,4,2,3,4,1,2,3,4,1,2,3,4,0,1,2,3,4,5,5,6,5,6,5,6,5),
                       Y = c(12,11,11,10,10,9,9,9,9,8,8,8,8,8,7,7,7,7,6,6,6,6,5,5,5,5,4,4,4,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0,0,4,3,3,2,2,1,1,0))
  # get numbers for squares
  squarenums <- floor(en / 100000)
  # convert to letters
  squares <- apply(squarenums, 1, function(x) {
    idx <- which(lookup[, 3] == x[1] & lookup[, 4] == x[2])
    paste(lookup[idx, 1], lookup[idx, 2], sep = '')
  })
  # get Eastings/ Northings w/in square, to nearest meter
  nums <- round(en - squarenums * 100000) # to nearest metre
  # remove trailing 0s
  nums <- t(apply(nums, 1, function(x) {
    x / 10 ^ (which(sapply(1:6, function(i, x) any(x %% 10 ^ i != 0), x))[1] - 1)
  }))
  
  gridrefs <- cbind(squares, nums)
  apply(gridrefs, 1, paste, collapse = '')
}

os2en <- function (os)
{
  # convert vector of grid references to eastings/northings. Accounts for tetrads at any resolution (not just 2km)
  # lookup for grid square
  lookup <- data.frame(N = c('H','H','H','H','H','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','T','T','T','T','T','T','T','T'),
                       E = c('P','T','U','Y','Z','A','B','C','D','F','G','H','J','K','L','M','N','O','R','S','T','U','W','X','Y','Z','C','D','E','H','J','K','M','N','O','P','R','S','T','U','V','W','X','Y','Z','A','F','G','L','M','Q','R','V'),
                       X = c(4,3,4,3,4,0,1,2,3,0,1,2,3,4,0,1,2,3,1,2,3,4,1,2,3,4,2,3,4,2,3,4,1,2,3,4,1,2,3,4,0,1,2,3,4,5,5,6,5,6,5,6,5),
                       Y = c(12,11,11,10,10,9,9,9,9,8,8,8,8,8,7,7,7,7,6,6,6,6,5,5,5,5,4,4,4,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0,0,4,3,3,2,2,1,1,0))
  # lookup for tetrad
  lookup_tet <- data.frame(tetsquare = LETTERS[c(1:14, 16:26)],
                           X = 0.2 * rep(0:4, each = 5),
                           Y = 0.2 * rep(0:4, 5))
  
  squares <- substr(os, 1, 2)
  nums <- substr(os, 3, nchar(os))
  
  # are there DINTY tetrads in the gridrefs? If so, get them and clean from nums
  tet <- nchar(nums) %% 2 != 0
  tet_squares <- substr(nums[tet], nchar(nums[tet]), nchar(nums[tet]))
  if(!all(as.character(tet_squares) %in% lookup_tet[, 1])) stop('uneven grid reference and tetrad not recognised')
  nums[tet] <- substr(nums[tet], 1, nchar(nums[tet]) - 1)  
  
  # work out resolution
  len <- nchar(nums) / 2
  res <- 10 ^ (5 - len)
  # get (xy coordinates within 100km grid square)
  xy <- data.frame(E = as.numeric(substr(nums, 1, len)),
                   N = as.numeric(substr(nums, len + 1, 2 * len))) * res
  # add tetrad offset if needed
  tet_offset <- sapply(tet_squares, function(x) as.numeric(lookup_tet[which(lookup_tet[, 1] == x), 2:3]))
  tet_offset <- ifelse(length(tet_offset) == 0, 0, tet_offset)
  
  xy[tet, ] <- xy[tet, ] + tet_offset  * res[tet]
  # calculate 100km grid square offset
  offset <- t(sapply(squares, function(x) {
    as.numeric(lookup[lookup[, 1] == substr(x, 1, 1) & lookup[, 2] == substr(x, 2, 2), 3:4])
  }))
  # add to get lower-left coordinate
  ll <- offset * 100000 + xy
  # get upper right, accounting for tetrads
  ur <- ll + res * ifelse(tet, 0.2, 1)
  data.frame(LL = ll, UR = ur, centre = ll + (ur - ll) / 2, resolution = res * ifelse(tet, 0.2, 1))
}

en2poly <- function (coords)
{
  # matrix with columns giving ll.e, ll.n, ur.e, ur.n to SpatialPolygons
  getpoly <- function (x, ID) {
    x <- as.numeric(x)
    coords <- rbind(x[1:2], x[c(1, 4)], x[3:4], x[c(3, 2)], x[1:2])
    Polygons(list(Polygon(coords)), ID)
  }
  lis <- list()
  for(i in 1:nrow(coords)) lis <- c(lis, getpoly(coords[i, ], i))
  SpatialPolygons(lis, proj4string = osgb36())
}

extent2poly <- function(extent, proj4string)
  # convert an sp or raster extent object to a polygon,
  # given the correct projection
{
  coords <- rbind(c(extent@xmin, extent@ymin),
                  c(extent@xmin, extent@ymax),
                  c(extent@xmax, extent@ymax),
                  c(extent@xmax, extent@ymin),
                  c(extent@xmin, extent@ymin))
  poly <- Polygons(list(Polygon(coords)), 1)
  SpatialPolygons(list(poly), proj4string = proj4string)
}

# areaDensity <- function(polygons, raster)
#   # get polygon area weighting (in terms of number of pixels covered) on
#   # 0 (max pixels) to 1 (one pixel)  scale
# {
#   areas <- sapply(polygons@polygons, function(x) {
#     x@Polygons[[1]]@area
#   }) / prod(res(raster))
#   # set minimum to one pixel
#   areas[areas < 1] <-  1
#   # and scale from 0 to 1
#   areas <- areas - min(areas)
#   areas <- areas / max(areas)
#   1 - (areas / 3)
# }
# 
# 

featureDensity <- function(feature, raster, weights = 1, ...) {
  # Get density of features, masked by raster
  # (points/polys). 'feature' is passed to 'x' and 'raster' to 'y',
  # so these can take any value those can. 'weights' can be a single value,
  # a vector of the correct length, or the name of a field if 'feature'
  # is a Spatial*DataFrame.

  # set raster values to 0
  raster <- raster * 0
  #   raster[!is.na(raster[])] <- 0
  
  density <- rasterize(feature, raster, field = weights, fun = sum,
                       update = TRUE, ...)
  return (density)
}

gaussWindow <- function (n, sigma = n / 5)
  # generate a gaussian window for smoothing
  # n = width in cells of window (if even, 1 is added)
{
  if (n %% 2 == 0) n <- n + 1
  mat <- matrix(0, n, n)
  centre <- n / 2  + 0.5
  dist <- (col(mat) - centre) ^ 2 + (row(mat) - centre) ^ 2
  exp(-dist / (2 * sigma ^ 2))
}

importRasters <- function(path, as = brick, ext = '.grd')
  # import all rasters form the filepath (of a given type) into a stack or brick
{
  files <- mixedsort(list.files(path, pattern = ext))
  as(lapply(paste(path, files, sep = '/'), raster))
}

maxExtent <- function(list, margin = c(0, 0, 0, 0)) {
  # Given a list of matrices with latitudes (column 1) and longitudes (column 2),
  # return an extent object containing all the coordinates.
  # optionally add a margin to these, in the units of the coordinates
  exts <- sapply(list, function(x) c(range(x[, 2]), range(x[, 1])))
  ext <- c(min(exts[1, ]), max(exts[2, ]), min(exts[3, ]), max(exts[4, ]))
  extent(ext + margin * c(-1, 1, -1, 1))
}

osgb36 <- function()
{
  CRS('+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.999601272 +x_0=400000
      +y_0=-100000 +datum=OSGB36 +units=m +no_defs +ellps=airy
      +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894')
}

wgs84 <- function (projected = FALSE)
  # latlong, either projected or uprojected
{
  if (projected) {
    CRS("+init=epsg:3395") 
  } else {
    CRS('+init=epsg:4326')
  }
}

percCover <- function(raster, template, points, codes)
  # given discrete high res (raster1) and low res (raster2) images and points
  # calculate the % cover of each class of raster1 in the cell of raster2
  # identified by points. dropna drops raster2 cells with all na values
{
  pointras <- rasterize(points, template, mask = TRUE)
#   polys <- pixels2polys(pointras)
  polys <- rasterToPolygons(pointras)
  extr <- extract(raster, polys)
  cover <- function(x) sapply(codes, function(i, x) mean(x == i, na.rm = TRUE), x)
  perc <- t(sapply(extr, cover))
  colnames(perc) <- codes
  perc
}

# weighted standard deviation
sdWeighted <- function (x, weights, weighted_mean)
{
  n <- length(x)
  if (length(weights) != n) stop('x and weights must have the same length')
  num <- sum(weights * (x - weighted_mean) ^ 2)
  denom <- sum(weights) * (n - 1) / n
  sqrt(num / denom)
}

splitIdx <- function (n, maxn = 1000) {
  # get start and end indices to split a vector into bins of maximum length 'maxn'.
  names(n) <- NULL
  bins <- n %/% maxn + ifelse(n %% maxn > 0, 1, 0)
  start <- 1 + (1:bins - 1) * maxn
  end <- c(start[-1] - 1, n)
  lapply(1:bins, function(i) c(start[i], end[i]))
}
