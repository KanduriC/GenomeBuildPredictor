#' Predicts the genome build version of genomic track files, if any of the sequence coordinates in the input track file are unique to a genome build version.
#'
#' @description \code{predict_genome_build} Predicts the genome build version of genomic track files, if any of the sequence coordinates in the input track file are unique to a genome build version.
#'
#' @param input_track An input genomic track in the form of a GRanges object; import_genomic_track function could be used to generate this.
#' @param species Character string indicating the species of the input genomic track.
#' @return A character string indicating the genome build version of the genomic track.
#' @examples
#' build_version <- predict_genome_build(my_input_track,species='human')
#' @export

predict_genome_build <- function(input_track,species=c("human", "mouse", "drosophila")){
  library(GenomicRanges)
  if (missing(input_track)) {
    stop("No genomic track provided as input")
  }
  if (missing(species)) {
    stop("need to provide the species name of the input genomic track")
  }
  if (!all(species %in% c("human", "mouse", "drosophila"))) {
    stop("species can only be one among 'human', 'mouse' and 'drosophila'")
  }
  species <- match.arg(species)
  switch(species,human={predicted_build <- suppressWarnings(predict_human_build(input_track))},mouse={predicted_build <- suppressWarnings(predict_mouse_build(input_track))},drosophila={predicted_build <- suppressWarnings(predict_drosophila_build(input_track))})
  return(predicted_build)
}
