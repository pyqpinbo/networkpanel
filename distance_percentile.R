#' Build a network from a distance matrix and a distance percentile threshold
#'
#' This function build a network from a distance matrix
#' @param distance_matrix The distance matrix
#' @param percentile_value The percentile value from 0 to 1
#' @param normalize Whether the distance matrix is normalized
#' @keywords distance matrix
#' @export
#' @examples

network_build <- function(distance_matrix, percentile_value = 0.5, normalize = TRUE) {
    if (normalize)
        distance_matrix = dist_normalize(distance_matrix)
    distance_threshold = dist_percentile(distance_matrix, percentile_value)
    A = matrix(0, ncol(distance_matrix), nrow(distance_matrix))
    A[distance_matrix < distance_threshold] = 1
    library("igraph")
    graph.adjacency(A, mode="undirected", diag=F)
}


#' Obtain the distance percentile threshold from a distance matrix
#'
#' This function obtain the distance percentile threshold from a distance matrix
#' @param distance_matrix The distance matrix
#' @param percentile_value The percentile value from 0 to 1
#' @keywords distance matrix
#' @export
#' @examples

dist_percentile <- function(distance_matrix, percentile_value) {
    distance_matrix[is.na(distance_matrix)] = +Inf
    quantile(distance_matrix[upper.tri(distance_matrix)], probs = c(percentile_value))
}
