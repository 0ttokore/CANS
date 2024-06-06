#' Cluster Analysis of Nucleotide Sequences. K-Means method.
#'
#' @description The function performs cluster analysis of nucleotide sequences
#' using the k-means method. It takes as input a list of nucleotide sequences `sequences`,
#' the number of estimated clusters `k`, clustering accuracy `e` and the maximum
#' allowed number of iterations `m_max` and method for distance calculation `method_stringdistmatrix`.
#'
#'
#' @param sequences List of nucleotide sequences; have to be in format `matrix`.
#' @param k The number of clusters.
#' @param e Maximum value of the difference in clustering quality between two sequential iterations.
#' @param m_max The maximum number of iterations allowed.
#' @param method_stringdistmatrix Method for distance calculation. The default is "cosine", see [stringdist()] - methods.
#'
#' @return A list of objects and cluster numbers to which they relate,
#' indexes of central objects, as well as the number of iterations performed.
#' @seealso [stringdistmatrix()], [kmeans()].
#' @export
#' @importFrom stringdist stringdistmatrix
#'
#' @examples
#' result <- CANS_k_means(CANS_imitation(number_of_clusters = 2), k = 2, "qgram")
#' result <- CANS_k_means(CANS_imitation(number_of_clusters = 4), k = 4, e = 0.00001, "jw")
CANS_k_means <- function(sequences, k, e = 0.001, m_max = 100, method_stringdistmatrix = "cosine") {
  distances <- stringdistmatrix(sequences, sequences, method = method_stringdistmatrix)
  N <- length(sequences)

  Q <- 0
  Q_m1 <- Inf
  m <- 0

  centers <- sample(1:N, k, replace = FALSE)
  clusters <- matrix(0, nrow = N, ncol = 2)

  repeat {
    distances_to_centers <- t(apply(distances, 1, function(x) x[centers]))
    closest_clusters <- apply(distances_to_centers, 1, which.min)
    clusters <- cbind(1:N, closest_clusters)

    current_distance <- 0
    for (i in 1:k) {
      cluster_elements <- clusters[clusters[, 2] == i, 1]
      min_distance <- Inf

      for (j in 1:N) {
        total_distance <- sum(distances[j, cluster_elements])

        if (total_distance < min_distance) {
          min_distance <- total_distance
          centers[i] <- j
        }
      }
      current_distance <- current_distance + min_distance
    }

    Q_m1 <- Q
    Q <- current_distance
    m <- m + 1

    if(abs(Q - Q_m1) <= e | m >= m_max) {
      return(list(clusters = clusters, centers = centers, number_of_iterations = m))
    }
  }
}
