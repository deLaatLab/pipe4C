# pipe4C - a 4C-seq processing pipeline

A pipeline that processes multiplexed 4C-seq reads directly from FASTQ files. It generates files in a range of widely used formats to facilitate visualization and further data analysis using standard genome browsers and tools, including our recently developed peak caller for 4C-seq data peakC.

If you have any difficulties using the pipeline, please do not hesitate to contact us (delaatbioinf@hubrecht.eu).

## Prerequisites

- A Unix like shell (e.g. Bash v3.2+)
- Bowtie2 v2.3+ available from http://bowtie-bio.sourceforge.net/bowtie2/.
- SAMtools v1.3+ available from http://www.htslib.org/.
- R v3.5+ available from https://www.r-project.org/.
- The following R packages available from CRAN:
  - optparse
  - caTools
  - config
- The following R packages available from Bioconductor:
  - shortRead
  - genomicRanges
  - genomicAlignments
  - BSgenome of interest
- The peakC package available from https://github.com/deWitLab/peakC/.

## Installation

Download the latest version of the pipeline from this git repository using:

```
    $ wget https://github.com/deLaatLab/pipe4C/archive/master.zip
    $ unzip master.zip
    $ cd ./pipe4C-master
```
**note:** the pipe4C.R and functions.R files need to be placed in the same folder. 


## Files required to run the pipeline:
* Reads in (compressed) FASTQ format.
  * Illumina Sequencing Systems generate raw data files in binary base call (BCL) format. Illumina offers bcl2fastq conversion software to demultiplex (based on the index used in the non-reading primer) and convert BCL files. If multiple flow cell lanes have been used to sequence the library, a single Read 1 FASTQ file should be created per index per sequence run by combining the FASTQ files for each flow cell lane using a standard “cat” command after BCL conversion or by using the no-lane-splitting option when running bcl2fastq.
* Configuration file (conf.yml)
  * Global and system specific parameters (such as e.g. paths and genome assemblies installed) that are likely to remain constant across different runs of the pipeline are defined in the global configuration file (conf.yml). In each run the pipeline initially loads the parameters defined in this global configuration file, and then proceeds to load run specific parameters and the experiment specific data defined in a separate file (vpFile). The global configuration file can be edited using any standard text editor. 
  The list of parameters that need to be set at least once upon installation on a system in the global configuration file are shown in table 1.
  
  <BR>
  
  
| Name            | Description                                                                                                                                                                        |
|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| fragFolder      | Path to the folder containing the fragment end libraries of the reference genomes.                                                                                                 |
| normalizeFactor | Reads mapped to the 4C fragment end library are normalized to account for sequencing depth according to the normalizeFactor.                                                       |
| enzymes         | Enzyme names used in the viewpoint file and its corresponding recognition motif.                                                                                                   |
| genomes         | Genome names used in the viewpoint file plus corresponding BSgenome package.                                                                                                       |
| bowtie2         | Path to corresponding bowtie2 index of reference genome. The reference genome used to generate the index should match the reference genome that was used to generate the BSgenome. |
| maxY            | Maximal Y value in local 4C cis plot.                                                                                                                                              |
| plotView        | Number of bp to plot around viewpoint in local 4C cis plot.                                                                                                                        |
| xaxisUnit       | X-axis unit (Mb, Kb or bp).                                                                                                                                                        |
| plotType        | Plots will be saved as PDF or PNG                                                                                                                                                  |
| binSize         | Genome bin size used in the genome plot.                                                                                                                                           |
| qualityCutoff   | Q-score. Trim 3-end of all sequences using a sliding window as soon as 2 of 5 nucleotides has quality encoding less than the Q-score. 0 = no trimming.                             |
| trimLength      | Trim reads to defined capture length from 3’-end. 0 = no trimming.                                                                                                                 |
| minAmountReads  | Minimum amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed.                                                 |
| readsQuality    | Bowtie2 minimum quality mapped reads.                                                                                                                                              |
| mapUnique       | Extract uniquely mapped reads, based on the lack of XS tag.                                                                                                                        |
| cores           | Number of cores for parallelization.                                                                                                                                               |
| wSize           | The running mean window size.                                                                                                                                                      |
| nTop            | Top fragments discarded for normalization.                                                                                                                                         |
| nonBlind        | Only keep non-blind fragments.                                                                                                                                                     |
| wig             | Create wig files for all samples.                                                                                                                                                  |
| plot            | Create viewpoint coverage plot for all samples.                                                                                                                                    |
| genomePlot      | Create genomeplot for all samples (only possible if analysis is “all” in vpFile).                                                                                                  |
| tsv             | Create tab separated value file for all samples                                                                                                                                    |
| bins            | Count reads for binned regions.                                                                                                                                                    |
| mismatchMax     | The maximum number of mismatches allowed during demultiplexing.                                                                                                                    |
  
  **Table 1.** Description of parameters that can be defined in the configuration file.
  
  <BR>
  
* Viewpoint file
  * Experiment specific parameters for each 4C-seq experiment are organized in a viewpoint file. Parameters in this file are stored in a tab-delimited format, with each row containing information for a separate experiment: 
  
| expname    | primer               | firstenzyme | secondenzyme | genome | vpchr | vppos    | analysis | fastq           |
|------------|----------------------|-------------|--------------|--------|-------|----------|----------|-----------------|
| mESC_Sox2  | TTGCACCCGTCTTCTTGATC | DpnII       | Csp6I        | mm9    | 3     | 34547661 | all      | index1.fastq.gz |
| mESC_Mccc1 | GAGGGTAATTTTAGCCGATC | DpnII       | Csp6I        | mm9    | 3     | 35873313 | cis      | index1.fastq.gz |

**Table 2.** Example of a viewpoint file in which two experiments are demultiplexed from the same FASTQ file based on their primer sequence.

<BR>
 

|   Name            | Description                                                                                                                                                                                                                                                                                  |
|-------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| expname           | Unique experiment name.                                                                                                                                                                                                                                                                      |
| primer            | Primer sequence.                                                                                                                                                                                                                                                                             |
| firstenzyme       | First restriction enzyme name (nearest to reading primer)                                                                                                                                                                                                                                    |
| secondenzyme      | Second restriction enzyme name.                                                                                                                                                                                                                                                              |
| genome            | Reference genome of interest.                                                                                                                                                                                                                                                                |
| vpchr             | The chromosome that contains the viewpoint (See note 15).                                                                                                                                                                                                                                    |
| vppos             | Coordinate of viewpoint position. Any position within the VP can be used except the RE motifs (see note 15).                                                                                                                                                                                 |
| analysis          | The final output tables will contain all reads (all) or only the reads that have been mapped to the VP chromosome (cis). For most analysis cis is sufficient and the generated output files will be smaller and therefore easier to process on local computers.                              |
| fastq             | Name of the FASTQ file.                                                                                                                                                                                                                                                                      |
| spacer (optional) | Spacer length. Number of nt included as spacer in the primer to enable out of phase sequencing. Default =0. The spacer sequence will not be used for demultiplexing. If the spacer sequence is used as a barcode include the sequence in the primer sequence and set the spacer length to 0. |

**Table 3.** Description of parameters that are required in the viewpoint file for processing a 4C-seq experiment.




## Running the pipeline:

```
Rscript <path to pipe4C.R script> [any additional arguments]. 
```

A list of both required and optional parameters that are recognized by the pipe4C.R script are shown in table 4. Default values (except vpFile, fqFolder, outFolder and confFile) are stored in the configuration file. 

For example,
```
Rscript pipe4C.R –-vpFile [path to vpFile] --fqFolder [path to folder containing the FASTQ files] –-outFolder [path to output folder] --cores 8 --wig --plot --genomePlot
```
will run the pipeline using 8 cores and generates a wig file, a viewpoint plot and a genome plot as output, next to the default outputs.


| Name             | Description                                                                                                                            |
|------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| --vpFile *       | path to the viewpoint file.                                                                                                            |
| --fqFolder *     | path to the folder containing the FASTQ files.                                                                                         |
| --outFolder *    | path to the output folder.                                                                                                             |
| --confFile       | path to configuration file – default is conf.yml in folder containing the pipeline script.                                             |
| --mismatchMax    | The maximum number of mismatches allowed during demultiplexing.                                                                        |
| --qualityCutoff  | Q-score. Trim 3-end of all sequences using a sliding window as soon as 2 of 5 nucleotides has quality encoding less than the Q-score.  |
| --trimLength     | Trim reads to defined capture length from 3-end.                                                                                       |
| --minAmountReads | Minimum amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed.     |
| --readsQuality   | Bowtie2 minimum quality mapped reads.                                                                                                  |
| --mapUnique      | Extract uniquely mapped reads, based on the lack of XS tag.                                                                            |
| --cores          | Number of cores for parallelization.                                                                                                   |
| --wSize          | The running mean window size.                                                                                                          |
| --nTop           | Top fragments discarded for normalization.                                                                                             |
| --nonBlind       | Only keep non-blind fragments.                                                                                                         |
| --wig            | Create wig files for all samples.                                                                                                      |
| --plot           | Create viewpoint coverage plot for all samples.                                                                                        |
| --genomePlot     | Create genomeplot for all samples (only possible if analysis is “all” in vpFile).                                                      |
| --tsv            | Create tab separated value file for all samples                                                                                        |
| --bins           | Count reads for binned regions.                                                                                                        |

**Table 4.** Description of parameters that are recognized by the pipe4C.R script. * are required. 


# Step by step tutorial

In this tutorial we will explain how to run the pipeline and perform peakC analysis on the Sox2 4C-seq data from Geeven et al. (SRA files GSM2824300, GSM2824301 and GSM2824302). 

## Set up the pipeline

#### Download the pipeline including the example files
Download the pipeline and example files using the following command. 

    $ wget https://github.com/deLaatLab/pipe4C/archive/master.zip
    $ unzip master.zip
    $ cd ./pipe4C-master

> **note:** the pipe4C.R and functions.R files need to be placed in the same folder.
> The FASTQ files and viewpoint file can be found in the example folder.

#### Modify the configuration file (conf.yml) using any plain text editor

*	Change the location in which the fragmented genome will be generated (fragFolder).
*	Change the location in which the bowtie2 index is stored.

## Run the pipeline

    Rscript pipe4C.R --vpFile=./example/VPinfo.txt --fqFolder=./example/ --outFolder=./outF/ --cores 8 --plot --wig

This will run the pipeline using 8 cores and generates a viewpoint plot and wig file as output, next to the default outputs.

## Let's have a look at the generated files
#### The report file
A report file is generated that contains quality metrics of all the 4C-libraries processed by the pipeline.

The report file indicates that the quality of the experiment is good, as:
* ~80% of the reads map to the viewpoint chromosome (fragMappedCisPercCorr).
* ~75% of the reads mapping to the viewpoint chromosome maps within 1MB from the viewpoint (cov1Mb).
* Furthermore ~90% of the mappable DpnII fragment ends within 100kb from the viewpoint (capt100Kb) have at least one read.

#### The viewpoint plots

The viewpoint plots can be found in the PLOT folder. A coverage plot of the viewpoint region is generated based on the normalized and smoothened data. The genomic region that is visualized, the height of the Y-axis, the unit of the X-axis (Mb, kb, bp) the file type (PDF or PNG) and window size are defined and can be adjusted in the configuration file.  In addition the quality metrics are displayed in this figure for a quick impression of the data.


#### Wig files
The normalized smoothened data will be written to a wig file for visualization in genome browsers such as the UCSC Genome Browser ([https://genome.ucsc.edu/](https://genome.ucsc.edu/)) and IGV.

#### RDS files
The pipeline produces R-objects stored as rds files, which contain all mapped 4C-seq reads as well as the mapping statistics and viewpoint information relevant to the experiment. This allows for a more thorough interactive analysis and visualization of the data in R.

## Perform peak calling using PeakC in R

PeakC is a method that was designed to identify reproducible peaks in regions close to the viewpoint in 4C-seq data. It computes non-parametric statistics based on ranks of coverage of 4C fragment ends with respect to a background model. This model estimates the background contact frequency of both the region upstream and downstream of the viewpoint independently for each 4C-seq experiment individually and uses a statistical model to identify genomic regions that are significantly contacted.

To run peakC on the rds files created by the pipe4C pipeline, you need to load the peakC R-package (make sure it's installed first!) and source the pipe4C functions in R: 

```{r source}
library(peakC)
pipe4CFunctionsFile <- "functions.R"
source(pipe4CFunctionsFile)
```
Next, you need to select rds files from 4C-Seq experiments generated with the same VP that you want to analyze. 
> Select 4C-seq experiments generated using the same VP

```
resultsDir <- "./outF/RDS/"
setwd(resultsDir)
rdsFiles <- c("set_1_viewpoint_10_ESC_replicate_1.rds","set_1_viewpoint_10_ESC_replicate_2.rds", "set_1_viewpoint_10_ESC_replicate_3.rds")
```

Now run peakC, with the default parameters 

```
resPeakC <- doPeakC(rdsFiles = rdsFiles)
```

and plot the results

```
plot_C(resPeakC,y.max=750)
```
Finally, extract the identified peak regions from the peakC analysis from the objects and export the results in a BED file which can be used for visualization in a genome browser together with the generated wig files.

```
resPeaks <- getPeakCPeaks(resPeakC=resPeakC)
peaksFile <- "./outF/set_1_viewpoint_10_ESC_peakC_peaks.bed"
exportPeakCPeaks(resPeakC=resPeakC,bedFile=peaksFile,name="set_1_viewpoint_10_ESC_peakC_peaks")
```


