# function file for seegSDM package

# editing SEEG SDM
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
  
  # if any elements of df are NAs, return NAs
  if (any(is.na(df))) {
    
    results <- c(deviance = NA,
                 rmse = NA,
                 kappa = NA,
                 auc = NA,
                 sens = NA,
                 spec = NA,
                 pcc = NA,
                 kappa_sd = NA,
                 auc_sd = NA,
                 sens_sd = NA,
                 spec_sd = NA,
                 pcc_sd = NA,
                 thresh = NA)
    
  } else {
    
    # add an id column (needed for PresenceAbsence functions)
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
    
    
  }
  
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
        
        # handle the case that pwdSample returns NAs
        if (length(keep) == 0) {
          
          # if so, return NAs too
          preds[[i]] <- data.frame(PA = rep(NA, 3),
                                   pred = rep(NA, 3))
          
          # and issue a warning
          warning (paste0('failed to carry out pwd sampling in submodel ',
                          i))
          
        } else {
          
          # add an evaluation dataframe to list
          preds[[i]] <- data.frame(PA = y.data[keep],
                                   pred = pred[keep])
        }
        
      }
      
      # calculate cv statistics for all folds
      stats <- t(sapply(preds, calcStats))
      
      # return the mean of these
      return (colMeans(stats, na.rm = TRUE))
      
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
          return (colMeans(stats, na.rm = TRUE))
          
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

occurrence2SPDF <- function (occurrence, crs=wgs84(TRUE)) {
  # helper function to convert an occurrence dataframe
  # i.e. one which passes checkOccurrence into a SpatialPointsDataFrame object
  
  # get column numbers for coordinates
  coord_cols <- match(c('Longitude', 'Latitude'), colnames(occurrence))
  
  # check dates
  if ("Date" %in% names(occurrence) && class(occurrence$Date) != "Date") {
    occurrence$Date <- as.Date(occurrence$Date)
  }
  
  # convert to SPDF
  occurrence <- SpatialPointsDataFrame(occurrence[, coord_cols],
                                       occurrence,
                                       coords.nrs  = coord_cols,
                                       proj4string = crs)
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
      
      # define a function to mask out pixels which aren't to be reclassified
      keep <- function(cells, ...) {
        ifelse(cells %in% level_GAULs,
               cells,
               NA)
      }
      
      # use the function to mask out unwanted regions and speed up `zonal`
      ad <- calc(ad, keep)
      
      # extract values for each zone, aggregating by `fun`
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

extractBatch <- function(batch, covariates, factor, admin, admin_mode="average", load_stack=stack) {
  ## Extract a batch of occurrence data
  ## Takes account of synoptic or temporally resolved covariates, as well as, point and admin data
  
  ## Support functions
  classifyCovaraites <- function (covs) {
    parseDate <- function (string, partIndex) {
      # Parse a covariate sub-file name assuming "YYYY-MM-DD" format
      if (typeof(string) != "character") {
        return (NA)
      } else {
        raw_parts <- strsplit(string, "-")[[1]]
        parts <- strtoi(raw_parts)
        if (length(parts) > 3 || length(parts) < 1 || any(is.na(parts)) || any(parts <= 0)) {
          return (NA)
        } else {
          if (length(parts) >= 2 && (nchar(raw_parts[[2]]) != 2 || parts[[2]] > 12)) {
            return (NA)
          }
          if (length(parts) >= 3 && (nchar(raw_parts[[3]]) != 2 || parts[[2]] > 31)) {
            return (NA)
          }
          return (parts)
        }
      }
    }
    
    classifyDate <- function (string) {
      # Work out the temporal type of a covariate sub-file y=year, m=month, d=day
      parts <- parseDate(string)
      if (any(is.na(parts))) {
        return (NA)
      } else {
        numb_parts <- length(parts)   
        if (numb_parts == 1) {
          return ("y")
        } else if (numb_parts == 2) {
          return ("m")
        } else if (numb_parts == 3) {
          return ("d")
        } else {
          return (NA)
        }
      }
    }
    
    classifications <- lapply(covs, function (cov) {
      # Work out the temporal type of a covariate based on the names of its sub-files, s=single, y=year, m=month, d=day
      if (typeof(cov) == "list") {
        dateClasses <- lapply(names(cov), classifyDate)
        if (any(is.na(dateClasses))) {
          return (NA)
        } else if (all(dateClasses == "y")) { 
          return ("y") 
        } else if (all(dateClasses == "m")) {
          return ("m")
        } else if (all(dateClasses == "d")) {
          return ("d")
        } else {
          return (NA) 
        }
      } else if (typeof(cov) == "character" || class(cov)[[1]] == "RasterLayer"){ #RasterLayer or filepath string
        return ("s")
      } else {
        return (NA)
      }
    })  
    if (any(is.na(classifications))) {
      stop(simpleError("Not all covariates could be classified", call = NULL))
    }
    return (classifications)
  }
  
  extractSubBatch <- function (sub_batch, sub_batch_covs, sub_batch_factor) {
    # Create subBatch results matrix
    sub_batch_covs_values <- matrix(NA, nrow = nrow(sub_batch), ncol = length(sub_batch_covs))
    colnames(sub_batch_covs_values) <- names(sub_batch_covs) 
    
    points <- sub_batch$Admin == -999
    
    # if there are points
    if (any(points)) {
      # extract and add to the results
      sub_batch_covs_values[points, ] <- extract(load_stack(sub_batch_covs), sub_batch[points, ])
    }
    
    # if there are any polygons
    if (any(!points)) {
      # extract them, but treat factors and non-factors differently
      factor_covs <- sub_batch_factor == TRUE
      if (admin_mode == "average") {
        if (any(factor_covs)) {
          # if there are any factors, get mode of polygon
          sub_batch_covs_values[!points, factor_covs] <- extractAdmin(sub_batch[!points, ],
                                                                      load_stack(sub_batch_covs[which(factor_covs)]),
                                                                      admin, fun = 'modal')
        }
        if (any(!factor_covs)) {
          # if there are any continuous, get mean of polygon
          sub_batch_covs_values[!points, !factor_covs] <- extractAdmin(sub_batch[!points, ],
                                                                       load_stack(sub_batch_covs[which(!factor_covs)]),
                                                                       admin, fun = 'mean')
        }
      } else if (admin_mode == "random") {
        # Freya's "Random pixel" stuff here?
      } else if (admin_mode == "latlong") {
        sub_batch_covs_values[!points, ] <- extract(load_stack(sub_batch_covs), sub_batch[!points, ])
      } else {
        stop(simpleError("Unknown mode for admin covariate value extraction", call = NULL))
      }
    }
    return (sub_batch_covs_values)
  }
  
  extractTemporalValues <- function (batch_covs_values, cov_class_value, date_format, timestep_size) {
    # Perform time aware covariate extraction
    
    pickTemporalCovariate <- function (temporal_covariates, time=NA) {
      # Select the appropriate sub-file of each covariate for a given timestep
      suitable_covs <- temporal_covariates[which(names(temporal_covariates) <= time)]
      if (length(suitable_covs) == 0) {
        return (temporal_covariates[[min(names(temporal_covariates))]])
      } else {
        return (temporal_covariates[[max(names(suitable_covs))]])
      }
    }
    
    if (any(cov_class == cov_class_value)) {
      relevant_covs_idx <- which(cov_class == cov_class_value)
      covariate_of_interest <- covariates[relevant_covs_idx]
      covariate_of_interest_factor <- factor[relevant_covs_idx]
      time_min <- as.Date(cut(min(batch$Date), timestep_size))
      time_max <- as.Date(cut(max(batch$Date), timestep_size))
      timesteps <- seq(time_min, time_max, by=timestep_size)
      for (i in seq_along(timesteps)) {
        timestep <- format(timesteps[i], date_format)
        subset <- format(batch$Date, date_format) == timestep
        if (any(subset)) {
          # Pick the covariate subfiles for this timestep
          covariates_for_timestep <- lapply(covariate_of_interest, pickTemporalCovariate, time=timestep)
          # Extract the values
          batch_covs_values[subset, relevant_covs_idx] <- extractSubBatch(batch[subset, ], covariates_for_timestep, covariate_of_interest_factor)
        }
      }
    }
    return(batch_covs_values)
  }
  
  ## Classify the covariates
  cov_class <- classifyCovaraites(covariates)  
  
  ## Create an empty matrix for the extracted covariate records
  batch_covs_values <- matrix(NA, nrow = nrow(batch), ncol = length(covariates))
  colnames(batch_covs_values) <- names(covariates)
  
  ## Handle non-temporal covariates
  if (any(cov_class == "s")) {
    relevant_covs_idx <- which(cov_class == "s")
    batch_covs_values[, relevant_covs_idx] <- extractSubBatch(batch, covariates[relevant_covs_idx], factor[relevant_covs_idx])
  }
  
  ## Handle temporal covariates
  batch_covs_values <- extractTemporalValues(batch_covs_values, "y", "%Y", "year")
  batch_covs_values <- extractTemporalValues(batch_covs_values, "m", "%Y-%m", "month") 
  batch_covs_values <- extractTemporalValues(batch_covs_values, "d", "%Y-%m-%d", "day") 
  
  # build output data structure
  results <- cbind(PA = batch$PA,
                    Weight = batch$Weight,
                    batch@coords,
                    batch_covs_values)
  results <- as.data.frame(results)
  
  # if there are any factor covariates, convert the columns in the dataframe
  factor_vector <- unlist(factor)
  if (any(factor_vector)) {
    facts <- names(which(factor_vector))
    for (name in facts) {
      results[, name] <- factor(results[, name])
    }
  }
  
  return (results)
} 

selectLatestCovariates <- function(covariates, load_stack=stack) {
  ## For a mixed set of temporal and non-temporal raster paths, build a stack containing the most recent covariate sub-file for each covariate
  selected_covariates <- lapply(covariates, function (cov) {
    if (typeof(cov) == "list") {
      return (cov[[max(names(cov))]])
    } else {
      return (cov)
    }
  })
  
  return (load_stack(selected_covariates))
}

runBRT <- function (data,
                    gbm.x,
                    gbm.y,
                    pred.raster = NULL,
                    gbm.coords = NULL,
                    wt = NULL,
                    max_tries = 5,
                    verbose = FALSE,
                    tree.complexity = 4,
                    learning.rate = 0.005,
                    bag.fraction = 0.75,
                    n.trees = 10,
                    n.folds = 10,
                    max.trees = 10000,
                    step.size = 10,
                    method = c('step', 'perf', 'gbm'),
                    family = 'bernoulli',
                    gbm.offset = NULL,
                    ...)

# wrapper to run a BRT model with Sam's defaults
# and return covariate effects, relative influence,
# and a prediction map (on the probability scale).
# background points are weighted at 4 * presence points,
# mimicking prevalence of 0.2
{
  
  # get the required method
  method <- match.arg(method)
  
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
  
  # if using gbm.step
  if (method == 'step') {
    
    # we'll need to use a while loop to tweak the parameters ontil the
    # algorithm converges
    
    # set up for the while loop
    m <- NULL
    tries <- 0
    
    # if there's an offset, get it as a vector
    if (!is.null(gbm.offset)) {
      offset <- data[, gbm.offset]
    } else {
      # otherwise set it to null
      offset <- NULL
    }
    
    
    # try 'tries' times
    while (is.null(m) & tries < max_tries) {
      # fit the model, if it fails m will be NULL and the loop will continue
      m <- gbm.step(data = data,
                    gbm.x = gbm.x,
                    gbm.y = gbm.y,
                    offset = offset,
                    step.size = 10,
                    tree.complexity = tree.complexity,
                    verbose = verbose,
                    learning.rate = learning.rate,
                    bag.fraction = bag.fraction,
                    n.trees = n.trees,
                    n.folds = n.folds,
                    max.trees = max.trees,
                    plot.main = FALSE,
                    site.weights = wt,
                    keep.fold.models = TRUE, 
                    keep.fold.vector = TRUE,
                    keep.fold.fit = TRUE,
                    family = family,
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
    
  } else if (method == 'perf') {
    
    # set up formula
    # if the family is poisson and there's an offset, add it to the formula
    if (!is.null(gbm.offset)) {
      f <- data[, gbm.y] ~ . + offset(data[, gbm.offset])
    } else {
      f <- data[, gbm.y] ~ .    
    }
    
    # if using gbm.perf, run a single BRT with max.trees trees
    
    m <- gbm(f,
             distribution = family,
             data = data[, gbm.x],
             n.trees = max.trees,
             cv.folds = n.folds,
             interaction.depth = tree.complexity,
             verbose = verbose,
             shrinkage = learning.rate,
             bag.fraction = bag.fraction,
             weights = wt,
             ...)
    
    # and run the gbm.perf procedure to select the optimal
    # number of trees
    ntree <- gbm.perf(m,
                      plot.it = FALSE,
                      method = 'cv')
    
    # if the best number of trees is also the maximum number of trees,
    # issue a warning
    
    if (ntree == n.trees) {
      warning(paste0('The optimal number of trees by cross fold validation\
                     using gbm.perf was ',
                     ntree,
                     ', the same as the maximum number of trees.
                     You may get a better model fit if you increase n.trees.'))
    }
    
    # set the number of trees in the model to the optimal
    m$n.trees <- ntree
    
    } else {
      
      # set up formula
      # if the family is poisson and there's an offset, add it to the formula
      if (!is.null(gbm.offset)) {
        f <- data[, gbm.y] ~ . + offset(data[, gbm.offset])
      } else {
        f <- data[, gbm.y] ~ .    
      }
      
      # otherwise run a single model with n.trees trees
      m <- gbm(f,
               distribution = family,
               data = data[, gbm.x],
               n.trees = n.trees,
               cv.folds = n.folds,
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
  if (!(is.null(pred.raster))){
    pred = predict(pred.raster,
                   m,
                   type = 'response',
                   n.trees = m$n.trees)
  } else{
    pred <- NULL
  }
  
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

combinePreds <- function (preds,
                          quantiles = c(0.025, 0.975),
                          parallel = FALSE,
                          ncore = NULL)
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
  
  # parallel version
  if (parallel) {
    skipClusterManagment <- FALSE
    
    # if the user spcified the number of cores
    if (!is.null(ncore)) {
      
      # check the maximum number of cores
      max_cores <- detectCores()
      
      if (ncore > max_cores) {
        warning(paste0('ncore = ',
                       ncore,
                       'specified, but this machine only appears to have ',
                       max_cores,
                       ' cores so using ',
                       max_cores,
                       ' instead'))
        
        ncore <- max_cores
        
      }
      
    } else if (sfIsRunning()) {
      # if user didn't specify the number of cores and
      # if a snowfall cluster is running, use that number of cores
      
      skipClusterManagment <- TRUE
      
      # get the number of cpus
      ncore <- sfCpus()
      
      cat(paste0('\nrunning combinePreds on ',
                 ncore,
                 ' cores\n\n'))
    } else {
      # if no user specified valu or snowfall cluster running, run on 1 core
      warning("ncore wasn't specified and a snowfall cluster doesn't appear to be running\
              , so only running on one core")
      
      ncore <- 1
      
    }
    
    # set up function
    clusterFun <- function (x) {
      calc(x, fun = combine)
    }
    
    
    
    if(!skipClusterManagment) {
      # start the new snow cluster
      beginCluster(ncore)
      
      # run the function on the new cluster
      ans <- clusterR(preds, clusterFun)
    
      # turn off the new snow cluster
      endCluster()
    } else {
      # run the function on the existing cluster
      ans <- clusterR(preds, clusterFun, cl=sfGetCluster())
    }   
  } else {
    # otherwise run the function sequentially
    ans <- calc(preds, fun = combine)
  }
  
  # assign names to output
  names(ans)[1:2] <- c('mean', 'median')
  names(ans)[3:nlayers(ans)] <- paste0('quantile_', quantiles)
  
  # and return
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
    #npres <- sum(sub[, prescol]) 
    #nabs <- sum(1 - sub[, prescol])
    npres <- nrow(sub[sub[, prescol] >= 1,])
    nabs <- nrow(sub[sub[, prescol] == 0,])
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
  pixels <- notMissingIdx(raster)
  
  if (prob) {
    prob <- getValues(raster)[pixels]
  } else {
    prob <- NULL
  }
  
  # sample does stupid things if x is an integer, so sample an index to the
  # pixel numbers instead
  idx <- sample.int(n = length(pixels),
                    size = n,
                    replace = replace,
                    prob = prob)
  points <- pixels[idx]
  
  # get as coordinates
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

xy2AbraidSPDF <- function (pseudo, crs, pa, weight, date, admin=-999, gaul=NA, disease=NA) {
  # Coverts an XY matrix to a SpatialPointsDataFrame, with the standard ABRAID column set
  colnames(pseudo) <- c("Longitude", "Latitude")
  pseudo <- data.frame(pseudo,
                       "Weight"=weight, 
                       "Admin"=admin, 
                       "GAUL"=gaul,
                       "Disease"=disease,
                       "Date"=date,
                       "PA" = pa,
                       stringsAsFactors = FALSE)
  return (occurrence2SPDF(pseudo, crs=crs))
}

abraidBhatt <- function (pars,
                         occurrence,
                         covariates,
                         consensus,
                         admin,
                         factor,
                         load_stack=stack) {
  # A clone of 'extractBhatt' for use in abraid. 
  # It behaves the same as extractBhatt, but requires a named list or nested named list of covariates, to perform time aware extraction
  # Defensive checks are skipped.
  # WARNING: Future versions will likely return a weight value column & skip pesudo-presence generation
  
  # coerce pars into a numeric vector (lapply can pass it as a dataframe)
  pars <- as.numeric(pars)
  
  # get number of occurrence records
  no <- length(occurrence)
  
  # calculate
  na <- ceiling(pars[1] * no)
  np <- ceiling(pars[2] * no)
  mu <- pars[3]
  
  # start building the combined data set
  all <- occurrence
  all$PA <- rep(1, nrow(all))
  
  # pseudo-absences
  if (na > 0) {
    # modify the consensus layer (-100:100) to probability scale (0, 1)
    abs_consensus <- 1 - (consensus + 100) / 200
    
    # sample from it, weighted by consensus
    # (more likely in -100, impossible in +100)
    p_abs <- bgDistance(na, occurrence, abs_consensus, mu, prob = TRUE, spatial=FALSE, replace=TRUE)
    p_abs <- xy2AbraidSPDF(p_abs, crs(all), 0, NA, sample(occurrence$Date, na, replace=TRUE))
    all <- rbind(all, p_abs)
  }
  
  # pseudo-presences
  if (np > 0) {
    
    # modify consensus to the 0, 1 scale and threshold at 'threshold'
    pres_consensus <- calc(consensus,
                           fun = function(x) {
                             ifelse(x <= -25,
                                    0,
                                    (x + 100) / 200)
                           })
    
    # sample from it, weighted by consensus (more likely in 100, impossible
    # below threshold)
    p_pres <- bgDistance(np, occurrence, pres_consensus, mu, prob = TRUE, spatial=FALSE, replace=TRUE)
    p_pres <- xy2AbraidSPDF(p_pres, crs(all), 1, NA, sample(occurrence$Date, np, replace=TRUE))
    all <- rbind(all, p_pres)
  }
  
  # extract covariates 
  return (extractBatch(all, covariates, factor, admin, admin_mode="average", load_stack=load_stack))
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

runABRAID <- function (mode, 
                       disease,
                       occurrence_path,
                       extent_path,
                       supplementary_occurrence_path,
                       admin0_path,
                       admin1_path,
                       admin2_path,
                       covariate_path,
                       discrete,
                       verbose = TRUE,
                       max_cpus = 32,
                       load_seegSDM = function(){ library(seegSDM) },
                       parallel_flag = TRUE) {
  
  # Given the locations of: a csv file containing disease occurrence data
  # (`occurrence_path`, a character), a GeoTIFF raster giving the definitive
  # extents of the disease (`extent_path`, a character), a csv file containing 
  # disease occurrence data for other diseases (`supplementary_occurrence_path`,
  # a character), GeoTIFF rasters giving standardised admin units (`admin0_path`,
  # `admin1_path`, `admin2_path`) and GeoTIFF rasters giving the covariates to
  # use (`covariate_path`,  a character vector). Run a predictive model to produce
  # a risk map (and associated outputs) for the disease.
  # The file given by `occurrence_path` must contain the columns 'Longitude',
  # 'Latitude' (giving the coordinates of points), 'Weight' (giving the degree
  # of weighting to assign to each occurrence record), 'Admin' (giving the
  # admin level of the record - e.g. 1, 2 or 3 for polygons or -999 for points),
  # 'GAUL' (the GAUL code corresponding to the admin unit for polygons, or
  # NA for points) and 'Disease' a numeric identifer for the disease of the occurrence.
  # The file given by `supplementary_occurrence_path` must contain the columns 
  # 'Longitude', 'Latitude' (giving the coordinates of points), 'Admin' (giving the
  # admin level of the record - e.g. 1, 2 or 3 for polygons or -999 for points),
  # 'GAUL' (the GAUL code corresponding to the admin unit for polygons, or
  # NA for points) and 'Disease' a numeric identifer for the disease of the occurrence.
  # To treat any covariates as discrete variables, provide a logical vector
  # `discrete` with `TRUE` if the covariate is a discrete variable and `FALSE`
  # otherwise. By default, all covariates are assumed to be continuous.
  # Set the maximum number of CPUs to use with `max_cpus`. At present runABRAID
  # runs 64 bootstrap submodels, so the number of cpus used in the cluster will
  # be set at `min(64, max_cpus)`.
  
  # ~~~~~~~~
  # check inputs are of the correct type and files exist
  abraidCRS <- crs("+init=epsg:4326")
  modes <- readLines(system.file('data/abraid_modes.txt', package='seegSDM'))
  stopifnot(class(mode) == 'character' &&
              is.element(mode, modes))

  stopifnot(is.numeric(disease))

  stopifnot(class(occurrence_path) == 'character' &&
              file.exists(occurrence_path))
  
  stopifnot(class(extent_path) == 'character' &&
              file.exists(extent_path) && 
              compareCRS(raster(extent_path), abraidCRS))
    
  stopifnot(class(supplementary_occurrence_path) == 'character' &&
              file.exists(supplementary_occurrence_path))

  stopifnot(file.exists(admin0_path) && 
              compareCRS(raster(admin0_path), abraidCRS))
  
  stopifnot(file.exists(admin1_path) && 
              compareCRS(raster(admin1_path), abraidCRS))
  
  stopifnot(file.exists(admin2_path) && 
              compareCRS(raster(admin2_path), abraidCRS))
  
  stopifnot(class(verbose) == 'logical')
  
  stopifnot(class(unlist(discrete)) == 'logical' &&
              length(discrete == length(covariate_path)))
  
  stopifnot(is.function(load_seegSDM))
  
  stopifnot(is.logical(parallel_flag))
  
  stopifnot(names(discrete) == names(covariate_path))
  
  stopifnot(class(unlist(covariate_path, recursive=TRUE)) == 'character' &&
              all(file.exists(unlist(covariate_path, recursive=TRUE))) &&
              all(sapply(sapply(unlist(covariate_path, recursive=TRUE), raster), compareCRS, abraidCRS)))
  
  # ~~~~~~~~
  # load data
  
  # occurrence data
  occurrence <- read.csv(occurrence_path,
                         stringsAsFactors = FALSE)
  
  # check column names are as expected
  stopifnot(sort(colnames(occurrence)) == sort(c('Admin',
                                            'Date',     
                                            'Disease',
                                            'GAUL',
                                            'Latitude',
                                            'Longitude',
                                            'Weight')))
  
  # convert it to a SpatialPointsDataFrame
  # NOTE: `occurrence` *must* contain columns named 'Latitude' and 'Longitude'
  occurrence <- occurrence2SPDF(occurrence, crs=abraidCRS)
  
  # occurrence data
  supplementary_occurrence <- read.csv(supplementary_occurrence_path,
                         stringsAsFactors = FALSE)
  
  # check column names are as expected
  stopifnot(sort(colnames(supplementary_occurrence)) == sort(c('Admin',
                                                          'Date',    
                                                          'Disease',
                                                          'GAUL',
                                                          'Latitude',
                                                          'Longitude')))
  # Functions to assist in the loading of raster data. 
  # This works around the truncation of crs metadata in writen geotiffs.
  abraidStack <- function(paths) {
    s <- stack(paths)
    crs(s) <- abraidCRS
    extent(s) <- extent(-180, 180, -60, 85)
    return (s)
  }
  abraidRaster <- function(path) {
    r <- raster(path)
    crs(r) <- abraidCRS
    extent(r) <- extent(-180, 180, -60, 85)
    return (r)
  }
  
  # convert it to a SpatialPointsDataFrame
  # NOTE: `occurrence` *must* contain columns named 'Latitude' and 'Longitude'
  supplementary_occurrence <- occurrence2SPDF(supplementary_occurrence, crs=abraidCRS)
  
  # load the definitve extent raster
  extent <- abraidRaster(extent_path)
  
  # load the admin rasters as a stack
  # Note the horrible hack of specifying admin 3 as the provided admin 1.
  # These should be ignored since ABRAID should never contain anything other
  # than levels 0, 1 and 2
  admin <- abraidStack(list(
    "0"=admin0_path,
    "1"=admin1_path,
    "2"=admin2_path,
    "3"=admin1_path))

  # get the required number of cpus
  nboot <- 64
  ncpu <- min(nboot,
              max_cpus)
  
  # start the cluster
  sfInit(parallel = parallel_flag,
         cpus = ncpu)
  
  # load seegSDM and dependencies on every cluster
  sfClusterCall(load_seegSDM)
  
  if (verbose) {
    cat('\nseegSDM loaded on cluster\n\n')
  }
  
  # prepare absence data
  if (mode == 'bhatt') {
    sub <- function(i, pars) {
      # get the $i^{th}$ row of pars 
      pars[i, ]
    }
    
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
    
    # generate pseudo-data in parallel
    data_list <- sfLapply(par_list,
                          abraidBhatt,
                          occurrence = occurrence,
                          covariates = covariate_path,
                          consensus = extent,
                          admin = admin, 
                          factor = discrete,
                          load_stack = abraidStack)
    if (verbose) {
      cat('extractBhatt done\n\n')
    }
  } else if (mode == "all_bias") {
    presence <- occurrence
    presence <- occurrence2SPDF(cbind(PA=1, presence@data), crs=abraidCRS)
    absence <- supplementary_occurrence
    absence <- occurrence2SPDF(cbind(PA=0, absence@data[, 1:2], Weight=1, absence@data[, 3:6]), crs=abraidCRS)
    all <- rbind(presence, absence)
    
    # create batches
    batches <- replicate(nboot, subsample(all@data, nrow(all), replace=TRUE), simplify=FALSE)
    batches <- lapply(batches, occurrence2SPDF, crs=abraidCRS)
    
    if (verbose) {
      cat('batches ready for extract\n\n')
    }
    
    # Do extractions
    data_list <- sfLapply(batches,
             extractBatch,
             covariates = covariate_path,
             admin = admin, 
             factor = discrete)
  } else {
    exit(1)
  }
  
  if (verbose) {
    cat('extraction done\n\n')
  }
  
  # run BRT submodels in parallel
  model_list <- sfLapply(data_list,
                         runBRT,
                         gbm.x = 5:ncol(data_list[[1]]),
                         gbm.y = 1,
                         pred.raster = selectLatestCovariates(covariate_path, load_stack=abraidStack),
                         gbm.coords = 3:4,
                         verbose = verbose)
  
  if (verbose) {
    cat('model fitting done\n\n')
  }
  
  # get cross-validation statistics in parallel
  stat_lis <- sfLapply(model_list,
                       getStats)
  
  if (verbose) {
    cat('statistics extracted\n\n')
  }
  
  # combine and output results
  
  # make a results directory
  dir.create('results')
  
  # cross-validation statistics (with pairwise-weighted distance sampling)
  stats <- do.call("rbind", stat_lis)
  
  # keep only the relevant statistics
  stats <- stats[, c('auc', 'sens', 'spec', 'pcc', 'kappa',
		     'auc_sd', 'sens_sd', 'spec_sd', 'pcc_sd', 'kappa_sd')]

  # write stats to disk
  write.csv(stats,
            'results/statistics.csv',
            na = "",
            row.names = FALSE)
  
  # relative influence statistics
  relinf <- getRelInf(model_list)
  
  # append the names to the results
  relinf <- cbind(name = rownames(relinf), relinf)

  # output this file
  write.csv(relinf,
            'results/relative_influence.csv',
            na = "",
            row.names = FALSE)
  
  # marginal effect curves
  effects <- getEffectPlots(model_list)
  
  # convert the effect curve information into the required format
  
  # keep only the first four columns of each dataframe
  effects <- lapply(effects,
                    function (x) {
                      x <- x[, 1:4]
                      names(x) <- c('x',
                                    'mean',
                                    'lower',
                                    'upper')
                      return(x)
                    })
  
  # paste the name of the covariate in as extra columns
  for(i in 1:length(effects)) {
    
    # get number of evaluation points
    n <- nrow(effects[[i]])
    
    # append name to effect curve
    effects[[i]] <- cbind(name = rep(names(effects)[i], n),
                          effects[[i]])
  }
  
  # combine these into a single dataframe
  effects <- do.call(rbind, effects)
  
  # clean up the row names
  rownames(effects) <- NULL
  
  # save the results
  write.csv(effects,
            'results/effect_curves.csv',
            na = "",
            row.names = FALSE)
  
  # get summarized prediction raster layers
  
  # lapply to extract the predictions into a list
  preds <- lapply(model_list,
                  function(x) {x$pred})
  
  # coerce the list into a RasterStack
  preds <- stack(preds)
  
  # summarize predictions
  preds <- combinePreds(preds, parallel=parallel_flag)
  
  # stop the cluster
  sfStop()
  
  # get the width of the 95% confidence envelope as a metric of uncertainty
  uncertainty <- preds[[4]] - preds[[3]]
  
  # save the mean predicitons and uncerrtainty as rasters
  writeRaster(preds[[1]],
              'results/mean_prediction',
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  writeRaster(uncertainty,
              'results/prediction_uncertainty',
              format = 'GTiff',
              NAflag = -9999,
              options = c("COMPRESS=DEFLATE",
                          "ZLEVEL=9"),
              overwrite = TRUE)
  
  # return an exit code of 0, as in the ABRAID-MP code
  return (0)
}

#### function for generating a master mask from stack of rasters

masterMask <- function (rasters) {
  # given a stack of rasters
  # loop through rasters masking one layer by every other layer to create a master mask
  # returns master mask
  
  master <- rasters[[1]]
  
  for (i in 1:nlayers(rasters)){
    master <- mask(master, rasters[[i]])
  }
  
  return(master)
  
} 

### function which calculates the number of points falling within each pixel in a raster

rasterPointCount <- function (rasterbrick, coords, absence = NULL, extract=FALSE){
  
  # given a rasterbrick and a two-column matrix of coordinates 'coords' (in the order: x, y)
  # counts the number of points falling within each pixel in the rasterbrick
  # if absence = NULL: returns a three-column dataframe containing x and y coordinates for each pixel 
  # and the frequency of occurrence points for that pixel. The frequency for pixels containing no points is 0 
  # if absence != NULL: returns a three-column dataframe containing x and y coordinates for each pixel
  # and the frequency of occurrence points for that pixel, or 0 if the pixel contains a psuedo-absence point. 
  # Pixels containing no points are removed.
  # If any coordinates for occurrence and pseudo-absence points fall within the same pixel,
  # they are removed from the pseudo-absence dataset and a warning is issued.
  # If extract=TRUE, raster values for each pixel are extracted and returned in the dataframe
  
  # ~~~~~~~~~
  # make dataframe of cells containing at least one occurrence point
  # get raster cell ID for each set of coordinates 
  cell <- cellFromXY(rasterbrick, coords)
  # make table of counts for each cell ID
  tab <- table(cell)
  # convert table to a dataframe
  occ <- data.frame(tab)
  
  # ~~~~~~~~~
  # make dataframe of cells containing pseudo-absence points
  if (!(is.null(absence))) {
    
    # get raster cell ID for each set of coordinates
    cell_0 <- cellFromXY(rasterbrick, absence)
    # convert vector to a dataframe of unique cell IDs
    absences <- data.frame(unique(cell_0))
    # set count to 0 for all cell IDs
    absences$Freq <- 0
    # change column names to match occ
    names(absences)[names(absences)=='unique.cell_0.']<-'cell'
    
    # check if any occurrence and pseudo-absence points fall within the same pixel
    idx_overlap <- which(absences$cell %in% occ$cell)
    
    # if so, remove points from pseudo-absence dataset and issue a warning
    if (length(idx_overlap) > 0) {
      absences <- absences[-idx_overlap,]
      warning (length(idx_overlap), " pseudo-absence points were removed as they fell within the same pixels as occurrence points")
    }
    
    # combine with occurrence data
    dat_all <- rbind(absences, occ)
    
  } else {
    
    # make a dataframe of cells containing no points
    # make a vector of all raster cell IDs
    cell_NA <- c(1:ncell(rasterbrick))
    # get index of cell_IDs containing at least one point
    idx <- which(cell_NA[] %in% occ$cell)
    # remove these from cell ID vector 
    cell_NA <- cell_NA[-idx]
    # convert cell ID vector to a dataframe 
    no_points <- data.frame(cell_NA)
    # set count to 0 for all cell IDs if presence only data
    no_points$Freq <- 0
    # change column names to match occ
    names(no_points)[names(no_points)=='cell_NA']<-'cell'
    
    # ~~~~~~~~
    # combine data 
    dat_all <- rbind(no_points, occ) 
    
  }
  
  # get coords of each cell
  cell_coords <- xyFromCell(rasterbrick, as.numeric(dat_all$cell))
  
  # combine with other info
  dat_all <- cbind(dat_all, cell_coords)
  
  if (extract) {
    
    # get raster values
    raster_values <- extract(rasterbrick, cell_coords)
    
    # combined with other info
    dat_all <- cbind(dat_all, raster_values)
    
  }
  
  # remove cell_ID column
  dat_all$cell <- NULL
  
  return(dat_all)
  
}

### function for predicting to rasters using a fitted model from runBRT 

makePreds <- function (object, pred_covs) {
  
  # given the output from runBRT and a rasterbrick of covariates to predict to
  # the function makes (and returns) a prediction to the rasters based on the model from runBRT
  # note that column names for covariates used in runBRT must match prediction covariate names
  
  model_var_names <- object$model$var.names
  
  if (!(all(model_var_names %in% names(pred_covs)))){
    stop('covariate column names must match prediction covariate names')
  }
  
  # get model and n.trees from runBRT output
  model <- object$model
  n.trees <- object$model$n.trees
  
  # make prediction
  pred <- predict(pred_covs,
                  model,
                  type = 'response',
                  n.trees = n.trees)
  return(pred)        
  
}

getConditionalEffectPlots <- function(models,
                                      plot = FALSE,
                                      quantiles = c(0.025, 0.975),
                                      hold = NULL,
                                      value = NULL,
                                      ...) {
  
  # given a list of BRT model bootstraps (each an output from runBRT)
  # get the mean and quantiles (95% by default) conditional effect curves 
  # and return a list giving a matrix for each covariate and optionally plot the results.
  # option to specify a covariate value for which conditional effect curves are generated  
  # the dots argument allows some customisation of the plotting outputs
  
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
  
  getConditionalPredictions <- function(models, df_preds) {
    
    # get number of trees and model
    n.trees <- models$model$n.trees
    m <- models$model
    
    # make predictions
    p_tmp <- predict(m, df_preds, n.trees = n.trees)  
    
    return(p_tmp)
    
  }
  
  getConditionalEffect <- function (cov, models, hold = NULL, value = NULL){
    
    # get covariate values
    x <- sapply(models[[1]]$effects, '[[', 1)
    
    # get the mean of all covariate values
    cov_means <- vector()
    
    # loop through x getting the mean of each covariate
    for (i in 1:length(x)) {
      mean <- mean(as.numeric(x[[i]]))
      cov_means <- append(cov_means, mean)
    }
    
    # get levels
    levels <- getLevels(cov, models)
    
    if (is.null(levels)) {
      
      # make prediction dataframe 
      matrix_means <- do.call(rbind, replicate(100, cov_means, simplify = FALSE))
      df_means <- data.frame(matrix_means)
      
      # use range of x values for this covariate
      df_pred <- df_means
      x_values <- x[[cov]]
      df_pred[,cov] <- x_values
      
    } else {
      
      # make prediction dataframe 
      matrix_means <- do.call(rbind, replicate(length(levels), cov_means, simplify = FALSE))
      df_means <- data.frame(matrix_means)  
      
      # use range of x values for this covariate
      df_pred <- df_means
      x_values <- as.numeric(levels)
      df_pred[,cov] <- x_values
      
    }
    
    # if a covariate value is to be held constant, replace column with specified value
    # if a value has not been specified, returns an error
    if (!(is.null(hold))){
      
      if (is.null(value)){
        stop('a value for hold must be specified')
        
      } else {
        df_pred[,hold] <- value
      }
    }  
    
    colnames(df_pred) <- names
    
    # get conditional predictions for this covariate for each sub-model
    pred_list <- lapply(models, getConditionalPredictions, df_pred)
    
    # convert to dataframe
    y <- as.data.frame(pred_list)
    
    # give y names
    colnames(y) <- paste0('model', 1:ncol(y))
    
    # calculate the mean
    mn <- rowMeans(y, na.rm = TRUE)
    
    # and quantiles
    qs <- t(apply(y, 1, quantile, quantiles, na.rm = TRUE))
    
    # and return these alond with the raw data
    return (cbind(covariate = x_values,
                  mean = mn,
                  qs,
                  y))
    
  }
  
  # get the number of covariates
  ncovs <- length(models[[1]]$effects)
  
  # get covariate names
  names <- sapply(models[[1]]$effects, function (x) names(x)[1])
  
  # get conditional effects 
  effects <- lapply(1:ncovs, getConditionalEffect, models)
  
  if (plot) {
    
    for (i in 1:ncovs) {
      eff <- effects[[i]]
      
      if (is.null(getLevels(i, models))) {
        # if it's a continuous covariate do a line and CI region
        
        # blank plot
        plot(eff[, 2] ~ eff[, 1], type = 'n', ylim = range(eff[, 3:4]),
             xlab = names[i], ylab = paste('f(', names[i], ')', sep = ''))
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
    
  }
  
  return(effects)
  
}

