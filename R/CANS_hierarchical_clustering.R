#' Cluster Analysis of Nucleotide Sequences. Hierarchical clustering method.
#'
#' @description The function performs cluster analysis of nucleotide sequences
#' using the hierarchical clustering method. It takes as input a list of nucleotide
#'  sequences `sequences`, the number of clusters `number_of_clusters`, method for
#'  distance calculation `method_stringdistmatrix`and cluster linking method `method_hclust`.
#' @param sequences List of nucleotide sequences; have to be in format `matrix`.
#' @param number_of_clusters The number of clusters.
#' @param method_stringdistmatrix Method for distance calculation. The default is "cosine", see [stringdist()] - methods.
#' @param method_hclust The agglomeration method to be used. The default is "complete", see [hclust()] - metrics.
#'
#' @return A list of objects and cluster numbers to which they relate, cophenetic coefficient
#' @seealso [stringdistmatrix()], [hclust()], [dendextend()], [stats()].
#' @export
#' @importFrom stringdist stringdistmatrix
#' @importFrom dendextend cutree
#' @importFrom dendextend color_branches
#' @importFrom stats cophenetic
#' @importFrom stats rect.hclust
#' @importFrom stats cor
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom stats as.dendrogram
#'
#' @examples
#' \dontrun{
#' result <- CANS_hierarchical_clustering(CANS_imitation(number_of_clusters = 4),
#'  number_of_clusters = 4, "cosine")
#' result <- CANS_hierarchical_clustering(CANS_imitation(number_of_clusters = 3),
#'  number_of_clusters = 3, "jw")
#'  }
CANS_hierarchical_clustering <- function(sequences, number_of_clusters, method_stringdistmatrix = "cosine",
                                    method_hclust = "complete"){
  data <- as.dist(stringdistmatrix(sequences, sequences, method =
                                     method_stringdistmatrix))
  hcTree <- hclust(data, method = method_hclust)

  dend <- color_branches(as.dendrogram(hcTree), k = number_of_clusters)
  dend_plot <- plot(dend, main = "Dendrogram of nucleotide sequence clusters",
                    xlab = 'Objects', ylab = 'Similarity')
  rect.hclust(hcTree, k = number_of_clusters)
  return(list(clusters = cutree()(dend, number_of_clusters),
              cophenetic = cor(data, cophenetic(hcTree))))
}
