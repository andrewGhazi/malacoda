#' Generate a distance matrix from a matrix of annotations
#'
#' Given an nxd matrix of variant annotations, produce an nxn distance matrix
#' describing the inter-variant distances in annotation space
#'
#' @param annotation_dat an n x d data frame of annotations
#' @param log_distance a logical indicating to use the log1p of the distances (TRUE) or the raw euclidean distances (FALSE)
#'
#'
#' @importFrom magrittr %>%
generate_distance_matrix = function(annotations,
                                    log_distance = TRUE){

  annotations =
  if (log_distance) {
    annotations %>%
      as.data.frame() %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix %>%
      log1p
  } else {
    annotations %>%
      as.data.frame() %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix
  }
}

find_prior_weights = function(){

}

fit_marg_prior = function(mpra_data){

}

fit_cond_prior = function(mpra_data, annotations){

}
