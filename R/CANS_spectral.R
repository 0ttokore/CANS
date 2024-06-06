#' Cluster Analysis of Nucleotide Sequences. Spectral clustering.
#'
#' @description The function performs cluster analysis of nucleotide sequences
#' using the spectral method. It takes as input a list of nucleotide sequences `sequences`,
#' the number of significant neighboring objects `number_of_neighbors`, method for
#' distance calculation `method_stringdistmatrix` and the number of estimated clusters
#' `number_of_clusters`.
#'
#' @param sequences List of nucleotide sequences; have to be in format `matrix`.
#' @param number_of_neighbors The number of significant neighboring objects.
#' Recommended to select at least 10. The default is 11.
#' @param method_stringdistmatrix Method for distance calculation. The default is "cosine", see [stringdist()] - methods.
#' @param number_of_clusters The number of clusters.
#'
#' @return A list of objects and cluster numbers to which they relate,
#' eigenvectors on which the transition to the reduced feature space is performed.
#' @seealso [stringdist()], [eigen()], [kmeans()].
#' @export
#' @importFrom stringdist stringdistmatrix
#' @importFrom stats kmeans
#'
#' @examples
#' result <- CANS_spectral(CANS_imitation(number_of_clusters = 2),
#'     method_stringdistmatrix = "qgram", 2)
#' result <- CANS_spectral(CANS_imitation(number_of_clusters = 4),
#'     method_stringdistmatrix = "cosine", 4)
CANS_spectral <- function(sequences, number_of_neighbors = 11, method_stringdistmatrix = "cosine", number_of_clusters = 2) {
  distances <- stringdistmatrix(sequences, sequences, method = method_stringdistmatrix)
  N <- length(sequences)
  neighbors <- matrix(0, nrow = N, ncol = N)

  for (i in 1:N) {
    index <- order(distances[i,])[2:number_of_neighbors]
    neighbors[i,][index] <- 1
  }

  neighbors = neighbors + t(neighbors)
  neighbors[ neighbors == 2 ] = 1

  D = colSums(neighbors)

  laplacian = diag(N) - diag(D^(-1/2)) %*% neighbors %*% diag(D^(-1/2))

  eigenvectors = eigen(laplacian, symmetric = TRUE)
  n = nrow(laplacian)
  eigenvectors = eigenvectors$vectors[,(n - 2):(n - 1)]

  clusters = kmeans(eigenvectors, number_of_clusters)
  return(list(clusters = clusters$cluster, eigenvectors = eigenvectors))
}
