# GenomeBuildPredictor
Predicts the genome build version of genomic track files, if any of the sequence coordinates in the input track file are unique to a genome build version

### Introduction

This tool builds upon the strengths of existing Bioconductor packages that enable scalable genomic analyses to predict the genome build version of genomic track files. For more information about scalable genomics with Bioconductor, see references [1,2].

### Tool usage

See the installation instructions below. The package contains two main functions `import_genomic_track` and `predict_genome_build`. `import_genomic_track` is a wrapper around rtracklayer's import function that allows the import of bed, gff and wig files as GRanges object. In addition, the `is.Peak` argument in import function takes two possible values `narrowPeak` or `broadPeak`, which allows the input of narrowPeak and broadPeak files (see in examples below).  `predict_genome_build` attempts to predict the genome build currently for three species: human, mouse, drosophila. For the specific arguments, type `?predict_genome_build` or `?import_genomic_track` in R console.

#### Examples

For this example, we use this bed file : 

http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz

```
library(GenomeBuildPredictor)
mytrack <- import_genomic_track('wgEncodeBroadHmmGm12878HMM.bed',file_format = 'bed')
predict_genome_build(mytrack,species = 'human')
```
To import narrowPeak or broadPeak files, use the `is.Peak` argument that takes either `narrowPeak` or `broadPeak` as input.

```
mytrack <- import_genomic_track('toy_narrowPeak_file.bed',file_format = 'bed',is.Peak='narrowPeak')
predict_genome_build(mytrack,species = 'human')
```

### Installation instructions (for now)

Firstly, this tool requires the latest version of R. Please update your R version if it is not 3.3.2. Secondly, I haven't tested the portability yet. At the moment, [download the repository here](https://github.com/KanduriC/GenomeBuildPredictor/archive/master.zip), unzip it, and install from source as follows: `install.packages('GenomeBuildPredictor-master/',repos=NULL,type="source")`

#### Troubleshooting:

If the installation went fine, the packages 'GenomicRanges'and 'rtracklayer' should also be installed. This can be tested by `library(GenomicRanges)` and `library(rtracklayer)`, which should not raise any error or warning message. Otherwise, install them like this: `source("https://bioconductor.org/biocLite.R")` followed by `biocLite("rtracklayer","GenomicRanges")`.

### References

1. https://www.ncbi.nlm.nih.gov/pubmed/25633503
2. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5181792/
