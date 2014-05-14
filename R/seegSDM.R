# function file for seegSDM package


## TO DOCUMENT

notMissingIdx <- function(raster) {
  # return an index for the non-missing cells in raster
  which(!is.na(getValues(raster)))
}

missingIdx <- function(raster) {
  # return an index for the missing cells in raster
  which(is.na(getValues(raster)))
}

# clone of the auc function in PresenceAbsence
# but without the shocking
auc2 <- function (DATA,
                  st.dev = TRUE,
                  which.model = 1,
                  na.rm = FALSE) {
  if (is.logical(st.dev) == FALSE) {
    stop ("'st.dev' must be of logical type")
  }
  if (is.logical(na.rm) == FALSE) {
    stop ("'na.rm' must be of logical type")
  }
  if (sum(is.na(DATA)) > 0) {
    if (na.rm == TRUE) {
      NA.rows <- apply(is.na(DATA), 1, sum)
      warning (length(NA.rows[NA.rows > 0]), " rows ignored due to NA values")
      DATA <- DATA[NA.rows == 0, ]
    } else {
      return (NA)
    }
  }
  if (length(which.model) != 1) {
    stop ("this function will only work for a single model, 'which.model' must be of length one")
  }
  if (which.model < 1 || round(which.model) != which.model) {
    stop ("'which.model' must be a positive integer")
  }
  if (which.model + 2 > ncol(DATA)) {
    stop ("'which.model' must not be greater than number of models in DATA")
  }
  DATA <- DATA[, c(1, 2, which.model + 2)]
  DATA[DATA[, 2] > 0, 2] <- 1
  OBS <- DATA[, 2]
  PRED <- DATA[, 3]
  if (length(OBS[OBS == 1]) == 0 || length(OBS[OBS == 1]) == 
        nrow(DATA)) {
    if (st.dev == FALSE) {
      return (NaN)
    } else {
      return (data.frame(AUC = NaN, AUC.sd = NaN))
    }
  }
  rm(DATA)
  PRED.0 <- PRED[OBS == 0]
  PRED.1 <- PRED[OBS == 1]
  N <- length(PRED)
  n0 <- as.double(length(PRED.0))
  n1 <- as.double(length(PRED.1))
  R <- rank(PRED, ties.method = "average")
  R0 <- R[OBS == 0]
  R1 <- R[OBS == 1]
  U <- n0 * n1 + (n0 * (n0 + 1)) / 2 - sum(R0)
  AUC <- U / (n0 * n1)
  
  # This line turns terrible predictions into excellent ones.
  # Deleted because JESUS. CHRIST.
  #   AUC[AUC < 0.5] <- 1 - AUC
  rm(PRED)
  rm(OBS)
  if (st.dev == FALSE) {
    return (AUC = AUC)
  } else {
    RR0 <- rank(PRED.0, ties.method = "average")
    RR1 <- rank(PRED.1, ties.method = "average")
    pless.0 <- (R0 - RR0) / n1
    pless.1 <- (R1 - RR1) / n0
    var.0 <- var(pless.0)
    var.1 <- var(pless.1)
    var.AUC <- (var.0 / n0) + (var.1 / n1)
    st.dev.AUC <- var.AUC ^ 0.5
    return (data.frame(AUC = AUC, AUC.sd = st.dev.AUC))
  }
}

# Calculate range of validation statistics for fitted models.
# df is a dat.frame or matrix containing observed 0 or 1 data in
# the first column and (0, 1] predictions in the second
calcStats <- function(df) {
  
  # ad an id column (needed for PresenceAbsence functions)
  df <- data.frame(id = 1:nrow(df), df)
  
  # ~~~~~~~~~~~
  # copntinuous probability metrics
  
  # bernoulli deviance
  dev <- devBern(df[, 2], df[, 3])
  
  # root mean squared error
  rmse <- rmse(df[, 2], df[, 3])

  # auc (using my safe version - see above)
  auc <- auc2(df, st.dev = TRUE)
  
  # ~~~~~~~~~~~  
  # discrete classification metrics
  
  # calculate the 'optimum' threshold - one at which sensitivity == specificity
  opt <- optimal.thresholds(df, threshold = 101, which.model = 1, 
                            opt.methods = 3)
  
  # create confusiuon matrix at this threshold
  confusion <- cmx(df, threshold = opt[1, 2])

  # kappa (using threshold)
  kappa <- Kappa(confusion, st.dev = TRUE)
  
  # sensitivity and specificity using threshold
  sens <- sensitivity(confusion, st.dev = TRUE)
  spec <- specificity(confusion, st.dev = TRUE)
  
  # proportion correctly classified using threshold
  pcc <- pcc(confusion, st.dev = TRUE)
  
  # create results vector
  results <- c(deviance = dev,
               rmse = rmse,
               kappa = kappa[, 1],
               auc = auc[, 1],
               sens = sens[, 1],
               spec = spec[, 1],
               pcc = pcc[, 1],
               kappa_sd = kappa[, 2],
               auc_sd = auc[, 2],
               sens_sd = sens[, 2],
               spec_sd = spec[, 2],
               pcc_sd = pcc[, 2],
               thresh = opt[1, 2])
  
  # and return it
  return (results)
}


## DOCUMENTED

# get the mean stats for one model run either for the training data
# (if 'cv = FALSE') or as the mean of the statistics for the k-fold
# cross-validation statistics. if pwd  = TRUE use dismo::pwdSmaple
# to select evaluation points. dots is passed to pwdSample
getStats <- 
  function (object,
            cv = TRUE,
            pwd = TRUE,
            threshold = 1,
            ...) {

    # get observed data
    y.data <- object$model$data$y
    
    # squish covariates back into a matrix
    x.data <- matrix(object$model$data$x,
                     nrow = length(y.data))
    
    # then coerce to a dataframe
    x.data <- data.frame(x.data)
    
    # and give them back their names
    names(x.data) <- colnames(object$model$data$x.order)
    

    # if pair-wise distance sampling is required
    if (pwd) {
      
      # error handling
      
      # the user wants to do non-cross-validated pwd
      if (!cv) {
        # throw an error since this is undefined
        stop ('Can only run the pwd procedure id cv = TRUE. See ?getStats for details.')
      }
      
      # if coordinates weren't returned by runBRT
      if (is.null(object$coords)) {
        # throw an error asking for them nicely
        stop ('Coordinates were not available to run the pwd procedure.
              To include coordinates provide them to runBRT via the gbm.coords argument, otherwise rerun getStats setting pwd = FALSE')
      }
            
      # get fold models
      models <- object$model$fold.models
      n.folds <- length(models)
      
      # get fold vector
      fold_vec <- object$model$fold.vector
      
      # blank list to populate
      preds <- list()
      
      # loop through the folds
      for (i in 1:n.folds) {
        
        # fold-specific mask
        mask <- fold_vec == i
        
        # predicted probabilities to test set
        pred <- predict.gbm(models[[i]],
                            x.data,
                            n.trees = models[[i]]$n.trees,
                            type = 'response')
        
        # training presence point index
        train_p <- which(!mask & y.data == 1)
        
        # test presence point index
        test_p <- which(mask & y.data == 1)
        
        # test absence point index
        test_a <- which(mask & y.data == 0)
        
        x <- pwdSample(object$coords[test_p, ],
                       object$coords[test_a, ],
                       object$coords[train_p, ],
                       n = 1, tr = threshold, ...)

        
        keep_p <- which(!is.na(x[, 1]))
        keep_a <- na.omit(x[, 1])
        
        keep <- c(test_p[keep_p], test_a[keep_a])
        
        # add an evaluation dataframe to list
        preds[[i]] <- data.frame(PA = y.data[keep],
                                 pred = pred[keep])
        
      }
      
      # calculate cv statistics for all folds
      stats <- t(sapply(preds, calcStats))
      
      # return the mean of these
      return (colMeans(stats))

      # with return statement
      
    } else {  # close pwd if
      
      # if pair-wise sampling isn't required
      
      # if the overall training validation stats are required
      if (!cv) {
        
        model <- object$model
  
        # predicted probabilities 
        pred <- predict(model,
                        x.data,
                        n.trees = model$n.trees,
                        type = 'response')
        
        stats <- calcStats(data.frame(y.data,
                                      pred))
        
        return (stats)
        
      } else {  # close cv if
        # otherwise, for cv stats
        
        models <- object$model$fold.models
        
        n.folds <- length(models)
        
        # get fold vector
        fold_vec <- object$model$fold.vector
        
        # blank list to populate
        preds <- list()
        
        # loop through the folds
        for (i in 1:n.folds) {
          
          # fold-specific mask
          mask <- fold_vec == i
          
          # predicted probabilities for that model 
          pred <- predict.gbm(models[[i]],
                              x.data[mask, ],
                              n.trees = models[[i]]$n.trees,
                              type = 'response')
          
          # add an evaluation dataframe to list
          preds[[i]] <- data.frame(y.data[mask],
                                   pred)
          
        }
        
        # calculate cv statistics for all folds
        stats <- t(sapply(preds, calcStats))
        
        # return the mean of these
        return (colMeans(stats))
      } # close pwd = FALSE, cv else
    }# close pwd else
  }

# checking inputs
checkOccurrence <- function(occurrence,
                            consensus,
                            admin,
                            consensus_threshold = -25,
                            area_threshold = 1,
                            max_distance = 0.05,
                            spatial = TRUE,
                            verbose = TRUE) {
  
  # if occurrence
  
  
  # check that the consensus layer is projected
  if (CRSargs(projection(consensus, asText = FALSE)) !=
        CRSargs(wgs84(TRUE))) {
    stop ('consensus is not projected,
          please project it (or correct the coordinate system)')
  }
  
  # expected column names and data types
  expected_names <- c('UniqueID',
                      'Admin',
                      'Year',
                      'Longitude',
                      'Latitude',
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
  
  # get locations of expected names
  occurrence_name_idx <- which(names(occurrence) %in% expected_names)
  
  # get data types of columns in occurrence
  occurrence_classes <- sapply(occurrence[, occurrence_name_idx], class)
  
  # do they match up?
  classes_match <- occurrence_classes == expected_classes#_ordered
  
  # if not, stop with a helpful error
  if (!all(classes_match)) {
    
    # list problems
    message_vector <- sapply(which(!classes_match),
                             function(i) paste(names(occurrence)[i],
                                               'is',
                                               occurrence_classes[i],
                                               'but should be',
                                               expected_classes[i]))
    
    stop(paste("data types don't match,\n",
               paste(message_vector, collapse = '\n')))
  }
  
  
  # ~~~~~~~~~~~~~~
  # check that Admin contains a reasonable number
  # (either admin level 0:3, or -999 for points)
  bad_admin <- !(occurrence$Admin %in% c(0:3, -999))
  
  if (any(bad_admin)) {
    stop (paste('bad Admin codes for records with these UniqueIDs:',
                paste(occurrence$UniqueID[bad_admin], collapse = ', ')))
  }  
  
  # ~~~~~~~~~~~~~~
  # remove polygons over area limit
  big_polygons <- occurrence$Admin != -999 &
    occurrence$Area > area_threshold
  
  # notify the user
  if (verbose) {
    cat(paste(sum(big_polygons),
              "polygons had areas greater than the area threshold of",
              area_threshold,
              "and will be removed.\n\n"))
  }
  
  # remove these from occurrence
  occurrence <- occurrence[!big_polygons, , drop = FALSE]
  
  # ~~~~~~~~~~~~~~
  # find points (and centroids) outside mask
  vals <- extract(consensus,
                  occurrence[, c('Longitude', 'Latitude')],
                  drop = FALSE)
  
  outside_mask <- is.na(vals)
  
  if (any(outside_mask)) {
    # if there are any outside the mask...
    
    # try to find new coordinates for these
    new_coords <- nearestLand(xyFromCell(consensus,
                                         which(outside_mask)),
                              consensus,
                              max_distance)
    
    # replace those coordinates in occurrence
    occurrence[outside_mask, c('Longitude', 'Latitude')] <- new_coords
    
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
                'decimal degrees out and have been removed.\n\n'))
    }
    
    # update outside_mask
    outside_mask[outside_mask] <- still_out
    
    # and remove still-missigng points from occurrence and vals
    occurrence <- occurrence[!outside_mask, , drop = FALSE]
  }
  
  # ~~~~~~~~~~~~~~
  # see if any coordinates have values below (or equal to)
  # the evidence consensus threshold
  consensus_scores <- extract(consensus,
                              occurrence[, c('Longitude', 'Latitude')])
  
  low_consensus <- consensus_scores <= consensus_threshold 
  
  # if there were any
  if (any(low_consensus)) {
    
    # let the user know
    if (verbose) {
      cat(paste('removing',
                sum(low_consensus),
                "points which were in areas with evidence consensus value below
                the threshold of",
                consensus_threshold,
                '\n\n'))
    }
    
    # then remove them
    occurrence <- occurrence[!low_consensus, , drop = FALSE]
  }
  
  # ~~~~~~~~~~~~~~
  # check that there are no duplicated polygon/year combinations
  
  # if there are any polygons
  if (any(occurrence$Admin != -999)) {
    
    # find and add GAUL codes, maintaining occurrence as a dataframe
    occurrence$GAUL <- getGAUL(occurrence2SPDF(occurrence), admin)$GAUL
    
    # get an index for which records are polygons
    poly_idx <- occurrence$Admin != -999
    
    # check for duplicates (Admin, GAUL and Year all the same)
    poly_dup <- duplicated(occurrence[poly_idx, c('Admin', 'GAUL', 'Year')])
    
    # if any are duplicated give a warning but proceed
    if (any(poly_dup)) {
      stop (paste('there are',
                  length(poly_dup),
                  'duplicated polygon / year combinations
                  at records with UniqueId:',
                  paste(occurrence$UniqueID[poly_idx][poly_dup],
                        collapse = ', ')))
    }
  }
  
  
  # ~~~~~~~~~~~~~~
  # if spatial = TRUE, return an SPDF, otherwise the dataframe
  if (spatial) {
    # get column numbers for long and lat
    coord_cols <- match(c("Longitude", "Latitude"), names(occurrence))
    # build SPDF, retaining coordinate column numbers
    occurrence <- SpatialPointsDataFrame(occurrence[, coord_cols],
                                         occurrence,
                                         coords.nrs = coord_cols,
                                         proj4string = wgs84(TRUE))
  }
  
  # return corrected occurrence dataframe/SPDF
  return (occurrence)
  }

checkRasters <- function (rasters, template, cellbycell = FALSE) {
  
  # check whether the raster* object 'rasters' conforms to the 'template'
  # rasterLayer. By default the extent and projection are compared.
  # If 'cellbycell = TRUE' the pixel values are individually compared with
  # the template though this can be very slow & RAM-hungry with large rasters.
  # the function throws an error if anything is wrong and returns 'rasters'
  # otherwise.
  
  if (extent(rasters) != extent(template)) {
    stop('extents do not match, see ?extent')
  }
  
  if (CRSargs(projection(rasters, asText = FALSE)) !=
        CRSargs(projection(template, asText = FALSE))) {
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

tempStand <- function (occurrence, admin, verbose = TRUE) {
  
  # temporal standardisation
  # Expand the dataframe according to date ranges and check for
  # duplicate location/year combinations
  # occurrence should have a column named 'SourceID' giving unique
  # IDs for the multi-year source records, if this isn't presence,
  # one will be created. Columns named 'Start_Year' and 'End_Year'
  # must also be present. If a column named 'GAUL' is present it
  # is assumed to contain correct GAUL codes for the polygons (or
  # NA for points) otherwise a rasterbrick of admin units must be
  # supplied via the 'admin' argument ang getGAUL will be used to
  # add a 'GAUL' column. An 'Admin' column must also be present.
  #
  # Returns a list of the expanded dataframe and vector of duplicate
  # records (or NULL if there are none).
  
  # get number of rows  
  n <- nrow(occurrence)
  
  # if no source ID column, add one
  if (!('SourceID' %in% names(occurrence))) {
    occurrence$SourceID <- 1:n
    
    # tell the user
    if (verbose) {
      cat('SourceID column was missing, one has been added.\n\n')
    }
  }
  
  # if no GAUL column
  if (!('GAUL' %in% names(occurrence))) {
    # add one
    occurrence$GAUL <- getGAUL(occurrence2SPDF(occurrence), admin)$GAUL
    
    # and tell the user
    if (verbose) {
      cat('GAUL column was missing, one has been added using getGAUL.\n\n')
    }
    
    # check if any were missed
    failed_GAUL <- is.na(occurrence$GAUL[occurrence$Admin != -999])
    
    # and throw an error if so
    if (any(failed_GAUL)) {
      stop (paste(sum(failed_GAUL),
                  'polygon centroids fell outside the masked area and their
                  GAUL codes could not be determined. Try using nearestLand
                  to correct these points.\n\n'))
    }
    }
  
  # expand the dataframe  
  
  # get date ranges (1 if same year)
  range <- occurrence$End_Year - occurrence$Start_Year + 1
  
  # make index for repetitions (rep each index 'range' times)
  rep_idx <- rep(1:n, times = range)
  
  # calculate years (start_year + 0:range)
  years <- occurrence$Start_Year[rep_idx] + unlist(lapply(range, seq_len)) - 1
  
  # subset dataframe
  df <- occurrence[rep_idx, ]
  
  # remove start and end year columns
  df <- df[, !(names(df) %in% c('Start_Year', 'End_Year'))]
  
  # add years
  df$Year <- years
  
  # and a unique ID
  df$UniqueID <- 1:nrow(df)
  
  
  # check for duplicates
  
  # look for polygons
  polys <- df$Admin != -999
  
  # if there are any polygons
  if (any(polys)) {
    # get duplicates (Admin, GAUL and Year all the same)
    poly_dup <- duplicated(df[polys, c('Admin', 'GAUL', 'Year')])
    duplicated_polygons <- which(polys)[poly_dup]
  }
  
  # look for points
  pnts <- df$Admin == -999
  
  if (any(pnts)) {
    
    # get cell numbers in admin
    nums <- cellFromXY(admin, df[pnts, c('Longitude', 'Latitude')])
    
    # get duplicates (pixel and Year all the same)
    pnt_dup <- duplicated(cbind(nums, df$Year[pnts]))
    duplicated_points <- which(pnts)[pnt_dup]
  }
  
  return (list(occurrence = df,
               duplicated_polygons = duplicated_polygons,
               duplicated_points = duplicated_points))
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
                                  point = as.numeric(points[i, ]))
                           })
  
  return (t(sapply(neighbour_list, nearest, raster)))
}

occurrence2SPDF <- function (occurrence) {
  # helper function to convert an occurrence dataframe
  # i.e. one which passes checkOccurrence into a SpatialPointsDataFrame object
  
  # get column numbers for coordinates
  coord_cols <- match(c('Longitude', 'Latitude'), names(occurrence))
  
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
  if (CRSargs(projection(occurrence, asText = FALSE)) !=
        CRSargs(projection(admin, asText = FALSE))) {
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


extractAdmin <- function (occurrence, covariates, admin, fun = 'mean') {
  
  # get indices for polygons in occurrence
  poly_idx <- which(occurrence$Admin != -999)
  
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
      # this *should* speed up the zonal operation
      ad <- reclassify(ad,
                       cbind(level_all_codes,
                             rep(NA, length(level_all_codes))))
      
      # extract values for each zone, aggregating by 'fun'
      zones <- zonal(covariates, ad, fun = fun)
      
      # match them to polygons
      # (accounts for change in order and for duplicates)
      which_zones <- match(level_GAULs, zones[, 1])
      
      # add them to the results matrix
      ex[level_idx, ] <- zones[which_zones, -1]
      
    }
    
  }
  
  return (ex)
  
}


runBRT <- function (data, gbm.x, gbm.y, pred.raster,
                    gbm.coords = NULL,
                    wt = NULL,
                    max_tries = 5, verbose = FALSE,
                    tree.complexity = 4, learning.rate = 0.005,
                    bag.fraction = 0.75, n.trees = 10,
                    n.folds = 10, max.trees = 40000,
                    step.size = 10, step = TRUE, ...)
  
  # wrapper to run a BRT model with Sam's defaults
  # and return covariate effects, relative influence,
  # and a prediction map (on the probability scale).
  # background points are weighted at 4 * presence points,
  # mimicking prevalence of 0.2
{
  # keep running until model isn't null
  # (can happen with wrong leaning rate etc.)
  
  # calculate weights
  if (is.null(wt)) {
    # if it isn't given, assume full weight for all observations
    wt <- rep(1, nrow(data))
  } else if (class(wt) == 'function') {
    # otherwise, if it's a function of the presene-absence column
    # calculate the weights
    wt <- wt(data[, gbm.y])
  } else if (length(wt) == nrow(data)) {
    # otherwise use them directly as weights
    wt <- wt
  } else if (length(wt) == 1) {
    # if it's one long, use it as a column index
    wt <- data[, wt]
  } else {
    stop('wt must either be NULL, a function, a vector of weights or a column index')
  }
  
  if (step) {
    # if using gbm.step
    
    # set up for the while loop
    m <- NULL
    tries <- 0
    
    # try 'tries' times
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
                    site.weights = wt,
                    keep.fold.models = TRUE, 
                    keep.fold.vector = TRUE,
                    keep.fold.fit = TRUE,
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
    
  } else {
    # if not stepping, run a single BRT with n.trees trees
    
    m <- gbm(data[, gbm.y] ~ .,
             distribution = 'bernoulli',
             data = data[, gbm.x],
             n.trees = n.trees,
             interaction.depth = tree.complexity,
             verbose = verbose,
             shrinkage = learning.rate,
             bag.fraction = bag.fraction,
             weights = wt,
             ...)
  }
  
  # get effect plots
  effects <- lapply(1:length(gbm.x),
                   function(i) {
                     plot(m,
                          i,
                          return.grid = TRUE)
                   })
  
  # get relative influence
  relinf <- summary(m,
                   plotit = FALSE)
  
  # get prediction raster
  pred = predict(pred.raster,
                 m,
                 type = 'response',
                 n.trees = m$n.trees)
  
  # get coordinates
  if (is.null(gbm.coords)) {
    coords <- NULL
  } else {
    coords <- data[, gbm.coords]
  }
  
  # otherwise return the list of model objects
  # (predict grids, relative influence stats and prediction map)
  ans <- list(model = m,
              effects = effects,
              relinf = relinf,
              pred = pred,
              coords = coords)
  
  return (ans)
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

getEffectPlots <- function (models,
                            plot = FALSE,
                            quantiles = c(0.025, 0.975),
                            ...) {
  
  # given a list of BRT model bootstraps (each an output from runBRT)
  # get the mean and quantiles (95% by default) effect curves and the
  # effect curves for each submodel and return a list giving a matrix
  # for each covariate and optionally plot the results. The dots argument
  # allows some customisation of the plotting outputs
  
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
  
  getEffect <- function (cov, models, summarise) {
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
    
    # give y names
    colnames(y) <- paste0('model', 1:ncol(y))
    
    # calculate the mean
    mn <- rowMeans(y, na.rm = TRUE)
    # and quantiles
    qs <- t(apply(y, 1, quantile, quantiles, na.rm = TRUE))

    # and return these alond with the raw data
    return (cbind(covariate = x,
                  mean = mn,
                  qs,
                  y))
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

combinePreds <- function (preds, quantiles = c(0.025, 0.975))
  # function to calculate mean, median and quantile raster layers for
  # ensemble predictions given a rasterBrick or rasterStack where each
  # layer is a single ensemble prediction.
{
  
  # function to get mean, median and quantiles
  combine <- function (x) {
    ans <- c(mean = mean(x),
             median = median(x),
             quantile(x, quantiles, na.rm = TRUE))
    return (ans)
  }
  
  stopifnot(nlayers(preds) > 1)
  
  ans <- calc(preds, fun = combine)
  
  names(ans)[3:nlayers(ans)] <- paste0('quantile_', quantiles)
  
  return(ans)
  
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

subsample <- function (data,
                       n,
                       minimum = c(5, 5),
                       prescol = 1,
                       replace = FALSE,
                       max_tries = 10) {
  # get a random subset of 'n' records from 'data', ensuring that there
  # are at least 'minimum[1]' presences and 'minimum[2]' absences.
  # assumes by default that presence/absence code is in column 1 ('prescol')
  OK <- FALSE
  tries <- 0
  
  # until criteria are met or tried enough times
  while (!OK & tries < max_tries) {
    # take a subsample
    sub <- data[sample(1:nrow(data), n, replace = replace), ]
    # count presences and absences
    npres <- sum(sub[, prescol]) 
    nabs <- sum(1 - sub[, prescol])
    # if they are greater than or equal to the minimum, criteria are met
    if (npres >= minimum[1] & nabs >= minimum[2]) {
      OK <- TRUE
    }
    tries <- tries + 1
  }
  
  # if the number of tr=ies maxed out, throuw an error
  if (tries >= max_tries) {
    stop (paste('Stopped after',
                max_tries,
                'attempts,
                try changing the minimum number of presence and absence points'))
  }
  
  return (sub)
  }

bgSample <- function(raster,
                     n = 1000,
                     prob = FALSE,
                     replace = FALSE,
                     spatial = TRUE)
  # sample N random pixels (not NA) from raster. If 'prob = TRUE' raster
  # is assumed to be a bias grid and points sampled accordingly. If 'sp = TRUE'
  # a SpatialPoints* object is returned, else a matrix of coordinates
{
  pixels <- notMissingIdx(raster)  # which(!is.na(getValues(raster)))
  if (prob) {
    prob <- getValues(raster)[pixels]
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
  # (distance) of occurrence points (an sp object given by sp). The dots
  # argument can be used to pass options to bgSample. This is waaaay more
  # efficient than anything using gBuffer.
  
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
  if (CRSargs(projection(consensus, asText = FALSE)) !=
        CRSargs(wgs84(TRUE))) {
    stop ('consensus must be projected WGS84, try ?projectExtent and ?wgs84.')
  }
  if (CRSargs(projection(covariates, asText = FALSE)) !=
        CRSargs(wgs84(TRUE))) {
    stop ('covariates must be projected WGS84, try ?projectExtent and ?wgs84.')
  }
  if (CRSargs(projection(admin, asText = FALSE)) !=
        CRSargs(wgs84(TRUE))) {
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
    
    # sample from it, weighted by consensus
    # (more likely in -100, impossible in +100)
    p_abs <- bgDistance(na, occurrence, abs_consensus, mu, prob = TRUE, ...)
    p_abs_covs <- extract(covariates, p_abs)
    p_abs_data <- cbind(PA = rep(0, nrow(p_abs_covs)),
                        p_abs@coords,
                        p_abs_covs)
    
  } else {
    
    p_abs <- p_abs_data <- NULL
    
  }
  
  # pseudo-presences
  if (np > 0) {
    
    #     if (!exists('abs_consensus')) {
    #       
    #       # if a pseudo-absence consensus layer already exists
    #       # then save some computation
    #       pres_consensus <- 1 - abs_consensus
    #       
    #     } else {
    #       
    #       # otherwise create it from consensus
    #       pres_consensus <- (consensus + 100) / 200
    #       
    #     }
    
    #     # threshold it (set anything below threshold to 0)
    #     threshold <- (threshold + 100) / 200
    #     
    #     under_threshold_idx <- which(getValues(pres_consensus) <= threshold)
    #     pres_consensus <- setValuesIdx(pres_consensus, 0, index = under_threshold_idx)
    
    # convert presence consensus to the 0, 1 scale and threshold at 'threshold'
    pres_consensus <- calc(consensus,
                           fun = function(x) {
                             ifelse(x <= threshold,
                                    0,
                                    (x + 100) / 200)
                           })
    
    # sample from it, weighted by consensus (more likely in 100, impossible
    # below threshold)
    p_pres <- bgDistance(np, occurrence, pres_consensus, mu, prob = TRUE, ...)
    p_pres_covs <- extract(covariates, p_pres)
    p_pres_data <- cbind(PA = rep(1, nrow(p_pres_covs)),
                         p_pres@coords,
                         p_pres_covs)
    
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
  points <- occurrence$Admin == -999
  
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
                                                fun = 'modal')
    }
    if (any(!factor)) {
      # if there are any continuous, get mean of polygon
      occ_covs[!points, !factor] <- extractAdmin(occurrence,
                                                 covariates[[which(!factor)]],
                                                 admin,
                                                 fun = 'mean')
      
    }
  }
  
  # add a vector of ones and the coordinates
  occ_data <- cbind(PA = rep(1, nrow(occ_covs)),
                    occurrence@coords,
                    occ_covs)
  
  # combine the different datasets and convert to a dataframe
  all_data <- rbind(occ_data, p_abs_data, p_pres_data)
  all_data <- as.data.frame(all_data)
  
  # if there are any factor covariates, convert the columns in the dataframe
  if (any(factor)) {
    facts <- which(factor)
    for (i in facts) {
      # plus one to avoid the PA column
      all_data[, i + 3] <- factor(all_data[, i + 3])
    }
  }
  
  # if occurrence gave the names of the xy coordinate columns... use them,
  if (length(occurrence@coords.nrs) == 2) {
    # use them in all data
    colnames(all_data)[2:3] <- colnames(occurrence@data)[occurrence@coords.nrs]
  } else {
    # otherwise call them coord_x and coord_y
    colnames(all_data)[2:3] <- c('coord_x', 'coord_y')
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
  # duplicate a blank raster
  ras <- raster
  ras <- setValues(ras, 0)
  
  # weighted occurrence raster (areas for polys roughly sum to 1)
  areas <- sapply(polygons@polygons, function(x) {
    x@Polygons[[1]]@area
  }) / prod(res(raster))
  ras <- rasterize(polygons, ras, field = 1/areas,
                   fun = 'sum', update = TRUE)
  # get sigma in pixels
  sigma <- ceiling(sigma / mean(res(raster)))
  # use a window of 5 * sigma
  w <- gaussWindow(sigma * 5, sigma)
  print(system.time(ras <- focal(ras,
                                 w = w / sum(w),
                                 pad = TRUE,
                                 na.rm = TRUE)))
  
  
  # replace mask
  #   ras <- setValuesIdx(ras, 0, missingIdx(ras))
  ras <- calc(ras, fun = function(x) ifelse(is.na(x), 0, x))
  
  #   ras[missingIdx(ras)] <- 0
  mask(ras, raster)
}

# bufferMask <- function(feature, raster, km = 100, val = NULL, maxn = 1000,
#   parallel = FALSE, dissolve = FALSE)
#   # Returns a mask giving the non-NA areas of 'raster' which are within
#   # 'km' kilometres of 'feature'. If val is specified, non-NA values are
#   # assigned val, if null the values in (the first layer of) raster are
#   # used. gBuffer does poorly with very complex feature, so buffering is
#   # carried out in batches of size 'maxn'. If a snowfall cluster is
#   # running, setting 'parallel = TRUE', runs the buffering in parallel. 
#   # Uses rgeos::gBuffer and rgdal::spTransform
# 
#   # NOTE: bgDistance is MUCH FASTER & more memory safe for generating
#   pseudo-absences
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

runABRAID <- function (occurrence_path,
                       extent_path,
                       admin_path,
                       covariate_path,
                       discrete = rep(FALSE,
                                      length(covariate_path)),
                       verbose = TRUE,
                       max_cpus = 32,
                       load_seegSDM = function(){ library(seegSDM) }) {
  
  # Given the locations of: a csv file containing disease occurrence data
  # (`occurrence_path`, a character), an ASCII raster giving the definitive
  # extents of the disease (`extent_path`, a character), and ASCII rasters
  # giving standardised admin units (`admin_path`; these *must* be in
  # ascending order) and the ASCII rasters giving the covariates to use
  # (`covariate_path`, a character vector). Run a predictive model to produce
  # a risk map (and associated outputs) for the disease.
  # The file given by `occurrence_path` must contain the columns 'Longitude',
  # 'Latitude' (giving the coordinates of points), 'Weight' (giving the degree
  # of weighting to assign to each occurrence record), 'Admin' (giving the
  # admin level of the record - e.g. 1, 2 or 3 for polygons or -999 for points)
  # and 'GAUL' (the GAUL code corresponding to the admin unit for polygons, or
  # NA for points).
  # To treat any covariates as discrete variables, provide a logical vector
  # `discrete` with `TRUE` if the covariate is a discrete variable and `FALSE`
  # otherwise. By default, all covariates are assumed to be continuous.
  # Set the maximum number of CPUs to use with `max_cpus`. At present runABRAID
  # runs 64 bootstrap submodels, so the number of cpus used in the cluster will
  # be set at `min(64, max_cpus)`.
  
  # ~~~~~~~~
  # lambda functions  
  sub <- function(i, pars) {
    # get the $i^{th}$ row of pars 
    pars[i, ]
  }
  
  # ~~~~~~~~
  # check inputs are of the correct type and files exist
  stopifnot(class(occurrence_path) == 'character' &&
              file.exists(occurrence_path))
  
  stopifnot(class(extent_path) == 'character' &&
              file.exists(extent_path))
  
  stopifnot(class(admin_path) == 'character' &&
              all(file.exists(admin_path)))
  
  stopifnot(class(covariate_path) == 'character' &&
              all(file.exists(covariate_path)))
  
  stopifnot(class(verbose) == 'logical')

  stopifnot(class(discrete) == 'logical' &&
              length(discrete == length(covariate_path)))
  
  stopifnot(is.function(load_seegSDM))
  
  # ~~~~~~~~
  # load data
  
  # occurrence data
  occurrence <- read.csv(occurrence_path,
                         stringsAsFactors = FALSE)
  
  # check column names are as expected
  stopifnot(sort(colnames(occurrence)) == c('Admin',
                                            'GAUL',
                                            'Latitude',
                                            'Longitude',
                                            'Weight'))
  
  # convert it to a SpatialPointsDataFrame
  # NOTE: `occurrence` *must* contain columns named 'Latitude' and 'Longitude'
  occurrence <- occurrence2SPDF(occurrence)
  
  # load the definitve extent raster
  extent <- raster(extent_path)
  
  # load the admin rasters as a stack
  admin <- stack(admin_path)
  
  # load covariate rasters into a stack
  covariates <- stack(covariate_path)
  
  # set the coordinate systems for rasters as unprojected wgs84 (lat/long)
  projection(covariates) <- wgs84(FALSE)
  projection(extent) <- wgs84(FALSE)
  
  
  # ~~~~~~~~
  # Generate pseudo-absence data  
  
  # set up range of parameters for use in `extractBhatt`
  # use a similar, but reduced set to that used in Bhatt et al. (2013)
  # for dengue.
  
  # number of pseudo-absences per occurrence
  na <- c(1, 4, 8, 12)
  
  # number of pseudo-presences per occurrence
  np <- c(0, 0.025, 0.05, 0.075)

  # maximum distance from occurrence data
  mu <- c(10, 20, 30, 40)

  # get all combinations of these
  pars <- expand.grid(na = na,
                      np = np,
                      mu = mu)
  
  # convert this into a list
  par_list <- lapply(1:nrow(pars),
                     sub,
                     pars)
  
  # get the required number of cpus
  ncpu <- min(length(par_list),
              max_cpus)
  
  # start the cluster
  sfInit(parallel = TRUE,
         cpus = ncpu)
  
  # load seegSDM and dependencies on every cluster
  sfClusterCall(load_seegSDM)
  
  # generate pseudo-data in parallel
  data_list <- sfLapply(par_list,
                        extractBhatt,
                        occurrence = occurrence,
                        covariates = covariates,
                        consensus = extent,
                        admin = admin, 
                        factor = discrete)
  
  # run BRT submodels in parallel
  model_list <- sfLapply(data_list,
                         runBRT,
                         gbm.x = 4:6,
                         gbm.y = 1,
                         pred.raster = covariates,
                         gbm.coords = 2:3,
                         verbose = verbose)
  
  # get cross-validation statistics in parallel
  stat_lis <- sfLapply(model_list,
                       getStats)
  
  # stop the cluster
  sfStop()
  
  # combine and output results
  
  # cross-validation statistics (with pairwise-weighted distance sampling)
  stats <- do.call("rbind", stat_lis)
  
  write.csv(stats,
            'statistics.csv',
            row.names = FALSE)
  
  # relative influence statistics
  relinf <- getRelInf(model_list,
                      plot = FALSE)
  
  write.csv(relinf,
            'relative_influence.csv',
            row.names = TRUE)
  
  # effect plots
  
  # set up plotting layout
  rowcol <- n2mfrow(nlayers(covariates))
  
  # open plotting device
  png('effect_plots.png',
      height = rowcol[1] * 300,
      width = rowcol[2] * 300)
  
  # plot marginal effect curves
  getEffectPlots(model_list,
                 plot = TRUE)
  
  # turn off the plotting device
  dev.off()
  
  # get summarized prediction raster layers

  # lapply to extract the predictions into a list
  preds <- lapply(model_list,
                  function(x) {x$pred})

  # coerce the list into a RasterStack
  preds <- stack(preds)
  
  # summarize predictions
  preds <- combinePreds(preds)
  
  # get the width of the 95% confidence envelope as a metric of uncertainty
  uncertainty <- preds[[4]] - preds[[3]]
  
  # save the mean predicitons and uncerrtainty as rasters
  writeRaster(preds$mean,
              'mean_prediction',
              format = 'ascii')
  
  writeRaster(uncertainty,
              'prediction_uncertainty',
              format = 'ascii')

  return (NULL)
}