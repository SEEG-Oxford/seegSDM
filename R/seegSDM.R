# function file for seegSDM package

require(raster)

# for testing
# occurrences <- data.frame(ID = 1:3,
#                           lat = letters[1:3],
#                           long = rnorm(3),
#                           type = 'point',
#                           Area = runif(3))#, stringsAsFactors = FALSE)

# checking inputs
checkOccurrences <- function(occurrences, evidence, evidence_threshold = -100,
                             area_threshold = 1, verbose = TRUE) {
  
  # given a dataframe "occurrence" of occurrence records
  # and a rasterLayer "evidence" giving the evidence consensus map;
  #   - check column names of occurrence
  #   - remove polygons over the area limit
  #   - remove points or polygons the outside evidence consensus
  #   - for points outside mask, try to move them in otherwise reject
    
  # expected column names and data types
  expected_names <- c('ID',
                      'lat',
                      'long',
                      'type',
                      'Area')
  
  expected_classes <- c('integer',
                        'numeric',
                        'numeric',
                        'character',
                        'numeric')
  
  # check if any expected column names are missing
  missing <- expected_names[!expected_names %in% names(occurrences)]
  
  # if they are, throw an error 
  if (length(missing) > 0) {
    stop(paste('missing columns:', paste(missing, collapse = ', ')))
  }
  
  # check column data types
  
  # get data types of columns in occurrences
  occurrence_classes <- sapply(occurrences, class)
  
  # get the expected data types in the same order
  order <- match(expected_names, names(occurrences))
  expected_classes_ordered <- expected_classes[order]
  
  # do they match up?
  classes_match <- occurrence_classes == expected_classes_ordered
  
  # if not stop with a helpful error
  if (!all(classes_match)) {

    # list problems
    message_vector <- sapply(which(!classes_match),
                             function(i) paste(names(occurrences)[i],
                                               'is',
                                               occurrence_classes[i],
                                               'but should be',
                                               expected_classes_ordered[i]))
    
    stop(paste("data types don't match\n",
               paste(message_vector, collapse = '\n')))
  }
  
  # remove polygons over area limit
  big_polygons <- occurrences$type == 'polygon' &
    occurrences$Admin != -9999 &
    occurrences$Area > area_threshold
    
  if (verbose) {
    cat(paste(sum(big_polygons),
              "polygons had areas greater than the threshold of",
              area_threshold,
              "and will be removed.\n"))
  }

  occurrences <- occurrences[!big_polygons, , drop = FALSE]
  
  
  
  # find points and polygons below evidence consensus threshold or outside mask
  vals <- extract(evidence, occurrences[, c('lat', 'long')], drop = FALSE)
  outside_mask <- is.na(vals)
  if (verbose) {
    cat(paste(sum(outside_mask),
              "points were outside mask, attempting to find nearby land...\n"))
  }
  
  
  
  # ***
  
  # routine to move points inland here
  
  # ***
  
    
  # otherwise find NAs again
  values <- extract(evidence, occurrences[, c('lat', 'long')])
  outside_mask <- is.na(values)
  if (verbose) {
    cat(paste('removing',
              sum(outside_mask),
              "points which were outside mask and couldn't be reconciled.
              points: ",
              which(outside_mask),
              '\n'))
  }
  
  # and remove from occurrences and vals
  occurrences <- occurrences[!outside_mask, , drop = FALSE]
  values <- values[!outside_mask]
  
  # remove any points with values below the evidence consensus threshold
  low_evidence <- values <= evidence_threshold 
  
  if (verbose) {
    cat(paste('removing',
              sum(low_evidence),
              "points which were in areas with evidence below the threshold of",
              evidence_threshold,
              '\n'))
  }
  
  occurrences <- occurrences[!low_evidence, , drop = FALSE]

  # return corrected occurrences dataframe
  return (occurrences)
}