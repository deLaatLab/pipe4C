# pipe4C

## 4C-seq  processing pipeline

A pipeline that processes multiplexed 4C-seq reads directly from FASTQ files. It generates files in a range of widely used formats to facilitate visualization and further data analysis using standard genome browsers and tools, including our recently developed peak caller for 4C-seq data peakC.

## Prerequisites

*	A Unix like shell (e.g. Bash v3.2+).
*	Illumina base-calling software (bcl2fastq) available from http://www.illumina.com/.
*	The 4C-seq pipeline can be downloaded from https://github.com/deLaatLab/. 
*	Bowtie2 v2.3+ available from http://bowtie-bio.sourceforge.net/bowtie2/. 
*	SAMtools v1.3+ available from http://www.htslib.org/.
*	R v3.5+ available from https://www.r-project.org/.
*	The following R packages available from CRAN: 
  * optparse
  * caTools
  * config
* The following R packages available from Bioconductor:
  * shortRead
  * genomicRanges
  * genomicAlignments
  * Bsgenome of interest
* The peakC package available from https://github.com/deWitLab/peakC/.
