#' Predicts the genome build version of genomic track files for human genome, if any of the sequence coordinates in the input track file are unique to a genome build version.
#'
#' @description \code{predict_human_build} Predicts the human genome build version of genomic track files, if any of the sequence coordinates in the input track file are unique to a genome build version.
#'
#' @param input_track An input genomic track in the form of a GRanges object; import_genomic_track function could be used to generate this.
#' @return A character string indicating the genome build version of the genomic track.
#' @examples
#' build_version <- predict_human_build(my_input_track)

predict_human_build <- function(input_track){
  library(GenomicRanges)
  builds <- c('hg38','hg19','hg18')
  overlap_build_1 <- sum(countOverlaps(hg38, input_track,ignore.strand=T))
  overlap_build_2 <- sum(countOverlaps(hg19, input_track,ignore.strand=T))
  overlap_build_3 <- sum(countOverlaps(hg18, input_track,ignore.strand=T))
  overlapped_with_builds <- c(overlap_build_1,overlap_build_2,overlap_build_3)
  if(all(overlapped_with_builds==0)){
    predicted_build <- 'Cannot predict the genome build because none of the sequence coordinates in the input track file are unique to any of the human genome builds'
    } else {
  predicted_build <- builds[which.max(overlapped_with_builds)]
    }
  return(predicted_build)
}
