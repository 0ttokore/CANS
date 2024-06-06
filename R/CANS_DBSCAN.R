#' Cluster Analysis of Nucleotide Sequences. Density-based Spatial Clustering of
#' Applications with Noise (DBSCAN)
#'
#' @description The function performs cluster analysis of nucleotide sequences
#' using the DBSACN method.It takes as input a list of nucleotide sequences `sequences`,
#' e-neighborhood size `e`, the minimum number of neighbors required for an object
#' to be considered primary `k_min` and method for distance calculation `method_stringdistmatrix`.
#'
#'
#' @param sequences List of nucleotide sequences; have to be in format `matrix`.
#' @param e e-neighborhood size.
#' @param k_min The minimum number of neighbors required for an object
#' to be considered primary. Typically selected from a range of 3 to 9. The default is 3.
#' @param method_stringdistmatrix Method for distance calculation. The default is "cosine", see [stringdist()] - methods.
#'
#' @return A list of objects and cluster numbers to which they relate,
#' the `number of clusters` and the `distance` matrix.
#' @seealso [stringdist()].
#' @export
#' @importFrom stringdist stringdistmatrix
#'
#' @examples
#' result <- CANS_DBSCAN(CANS_imitation(number_of_clusters = 3), e = 4, k_min = 5, "cosine")
#' result <- CANS_DBSCAN(CANS_imitation(number_of_clusters = 4), e = 4, k_min = 3, "qgram")
CANS_DBSCAN <- function(sequences, e, k_min = 3, method_stringdistmatrix = "cosine") {
  distances <- stringdistmatrix(sequences, sequences, method = method_stringdistmatrix)
  N <- length(sequences)

  neighbors <- lapply(1:N, function(i) which(distances[i, ] <= e & distances[i, ] != 0))
  number_of_neighbors <- sapply(neighbors, length)

  obj_classes <- ifelse(number_of_neighbors >= k_min, 3,
                        ifelse(number_of_neighbors < k_min & number_of_neighbors > 0, 2, 1))

  clusters <- numeric(N)
  number_of_clusters <- 0

  while (any(obj_classes == 3 & clusters == 0)) {
    index <- sample(which(obj_classes == 3 & clusters == 0),1)
    number_of_clusters <- number_of_clusters + 1
    clusters[index] <- number_of_clusters

    selected_indices <- intersect(which(obj_classes == 3), neighbors[[index]])
    clusters[selected_indices] <- number_of_clusters

    flag <- 1

    while (flag == 1) {
      flag <- 0
      for (i in 1:N) {
        if (obj_classes[i] == 3 & clusters[i] == 0 | (any(clusters[neighbors[[i]]] == 0 & obj_classes[neighbors[[i]]] == 3))) {
          current_cluster <- unique(union(clusters[i], clusters[neighbors[[i]]][clusters[neighbors[[i]]] != 0]))[1]
          if (current_cluster != 0) {
            clusters[intersect(which(obj_classes == 3), neighbors[[i]])] <- current_cluster
            clusters[i] <- current_cluster
            flag <- 1
          }
        }
      }
    }
  }

  indices <- which(obj_classes == 2)
  clusters[indices] <- sapply(indices, function(i) unique(clusters[neighbors[[i]]][clusters[neighbors[[i]]] != 0])[1])

  return(list(clusters = clusters, number_of_clusters = number_of_clusters,
              distances = distances, neighbors = neighbors))
}
