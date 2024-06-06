#' Cluster Analysis of Nucleotide Sequences. Imitation model.
#'
#' @description Generates a list of clusters of nucleotide sequences  for performing
#' cluster analysis. Takes parameters as input: `number of sequences`, their `length`,
#' `number of clusters`, alphabet of characters `nucleotides` and `randomness` parameter.
#' The function creates random starting nucleotide sequences for each cluster,
#' as well as a vector of random nucleotide positions in the given sequences. Generates
#' new sequences based on previously generated ones by changing nucleotides at
#' certain positions. Thus, clusters of similar sequences are formed, with a share
#' of similar symbols of at least 1 - `randomness`.
#'
#'
#' @param number_of_sequences The number of sequences in one cluster. The default is 100.
#' @param sequences_length The length of each sequence. The default is 100.
#' @param number_of_clusters The number of clusters. The default is 1.
#' @param nucleotides The alphabet of characters. The default is c("`A`", "`T`", "`C`", "`G`").
#' @param randomness Determines the maximum proportion of differences between the
#' elements of each cluster. Varies from 0 to 1. The default is 0.4.
#'
#' @return A list of nucleaotide sequences.
#' @export
#'
#' @examples
#' \dontrun{
#' sequences <- CANS_imitation(number_of_sequences = 100, sequences_length = 80,
#' number_of_clusters = 3, randomness = 0.3)
#'
#' sequences <- CANS_imitation()
#' }
CANS_imitation <- function(number_of_sequences = 100, sequences_length = 80, number_of_clusters = 1,
                      nucleotides = c("A", "T", "C", "G"), randomness = 0.4){
  if (randomness < 0 | randomness > 1) {
    return(print("Randomness have to be in the range from 0 to 1!"))
  }

  sequences <- vector(mode = "list", length = number_of_sequences * number_of_clusters)
  current_cluster <- 1

  for (i in 1:number_of_clusters) {
    sequence <- replicate(sequences_length, sample(nucleotides, 1, replace = TRUE))
    start_of_new_cluster <- 1 + (current_cluster - 1) * number_of_sequences
    sequences[[start_of_new_cluster]] <- paste(sequence, collapse = "")

    positions <- replicate(round(randomness * sequences_length), sample(1:sequences_length, 1))
    for (j in seq(from = start_of_new_cluster + 1, to = start_of_new_cluster + number_of_sequences - 1)) {
      sequence[positions] <- sample(nucleotides, length(positions), replace = TRUE)
      sequences[[j]] <- paste(sequence, collapse = "")
    }
    current_cluster <- current_cluster + 1
  }
  return(sequences)
}
