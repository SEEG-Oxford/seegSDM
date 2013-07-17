# function file for seegSDM package

# checking inuts


checkOccurrences <- function(occurrences, evidence, evidence_threshold, area_threshold = 1) {
  
  # given a dataframe "occurrence" of occurrence records
  # and a rasterLayer "evidence" giving the evidence consensus map;
  #   - check column names of occurrence
  #   - remove polygons over the area limit
  #   - remove points or polygons the outside evidence consensus
  #   - for points outside mask, try to move them in otherwise reject
    
  # check column names
  expected_names <- c('ID', 'lat', 'long')
  missing_names <- expected_names[!expected_names %in% names(occurrences)]
  if (length(missing_names) > 0) {
    stop(paste(''))
  }
  
  # remove polygons over area limit
  occurrences < occurrences[occurrences$Area < area_threshold]
  
  #   polyidx <- occurrences$type == 'polygon'
  
  
  
  
  return (NULL)
}