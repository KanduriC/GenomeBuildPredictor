#' Predicts the genome build version of genomic track files for mouse genome, if any of the sequence coordinates in the input track file are unique to a genome build version.
#'
#' @description \code{predict_mouse_build} Predicts the mouse genome build version of genomic track files, if any of the sequence coordinates in the input track file are unique to a genome build version.
#'
#' @param input_track An input genomic track in the form of a GRanges object; import_genomic_track function could be used to generate this.
#' @return A character string indicating the genome build version of the genomic track.
#' @examples
#' build_version <- predict_mouse_build(my_input_track)

predict_mouse_build <- function(input_track){
  library(GenomicRanges)
  builds <- c('mm10','mm9','mm8')
  overlap_build_1 <- sum(countOverlaps(mm10, input_track,ignore.strand=T))
  overlap_build_2 <- sum(countOverlaps(mm9, input_track,ignore.strand=T))
  overlap_build_3 <- sum(countOverlaps(mm8, input_track,ignore.strand=T))
  overlapped_with_builds <- c(overlap_build_1,overlap_build_2,overlap_build_3)
  if(all(overlapped_with_builds==0)){
    predicted_build <- 'Cannot predict the genome build because none of the sequence coordinates in the input track file are unique to any of the mouse genome builds'
  } else {
    predicted_build <- builds[which.max(overlapped_with_builds)]
  }
  return(predicted_build)
}
