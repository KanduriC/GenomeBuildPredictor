#' Imports a genomic track file using the import functionality in rtracklayer package.
#'
#' @description \code{import_genomic_track} Imports a genomic track file using the import functionality in rtracklayer package.
#'
#' @param track_file Character string indicating the filename of input genomic track file
#' @param file_format Character string indicating the file format of input file like bed, gff, wig. For all formats, see import function in rtracklayer package.
#' @return A GRanges object of input track file.
#' @examples
#' mybedfile <- import_genomic_track('mybedfile.bed',format='bed')
#' @export

import_genomic_track <- function(track_file,file_format){
  library(rtracklayer)
  if (missing(track_file)) {
    stop("No genomic track file provided as input")
  }
  if (missing(file_format)) {
    stop("need to provide the genomic track file format")
  }
  if (file_format %in% c('gff','gff1','gff2','gff3','bed','wig', 'ucsc')) {

  input_track <- import(track_file, format=file_format)
    } else {
    if (file_format %in% c('narrowPeak')){
      narrowPeak_cols <- c(signalValue = "numeric", pValue = "numeric",
           qValue = "numeric", peak = "integer")
    input_track <- import(track_file, format = "bed",
                            extraCols = narrowPeak_cols)
    }
    if (file_format %in% c('broadPeak')){
      broadPeak_cols <- c(signalValue = "numeric", pValue = "numeric",
                               qValue = "numeric")
      input_track <- import(track_file, format = "bed",
                            extraCols = broadPeak_cols)
    }
    }
  if(length(grep('chr', head(seqnames(input_track)),ignore.case=T))!=6){
    stop("the sequence names should have a chr prefix like chr1, chr2")
  }
  return(input_track)
}

