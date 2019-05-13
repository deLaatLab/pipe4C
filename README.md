# pipe4C - a 4C-seq processing pipeline

A pipeline that processes multiplexed 4C-seq reads directly from FASTQ files. It generates files in a range of widely used formats to facilitate visualization and further data analysis using standard genome browsers and tools, including our recently developed peak caller for 4C-seq data peakC.

## Prerequisites

- A Unix like shell (e.g. Bash v3.2+)
- Bowtie2 v2.3+ available from http://bowtie-bio.sourceforge.net/bowtie2/.
- SAMtools v1.3+ available from http://www.htslib.org/.
- R v3.5+ available from https://www.r-project.org/.
- The following R packages available from CRAN:
  - ptparse
  - caTools
  - config
- The following R packages available from Bioconductor:
  - shortRead
  - genomicRanges
  - genomicAlignments
  - Bsgenome of interest
- The peakC package available from https://github.com/deWitLab/peakC/.

## Installation

Download the latest version of the pipeline from this git repository using:

```
    $ wget https://github.com/deLaatLab/pipe4C/archive/master.zip
    $ unzip master.zip
    $ cd ./pipe4C
```
**note:** the pipe4C.R and functions.R files need to be placed in the same folder. The global configuration file (conf.yml) can be renamed or put in a different location if required.


## Files required to run the pipeline:
* Reads in (compressed) FASTQ format.
  * Illumina Sequencing Systems generate raw data files in binary base call (BCL) format. Illumina offers bcl2fastq conversion software to demultiplex (based on the index used in the non-reading primer) and convert BCL files. If multiple flow cell lanes have been used to sequence the library, a single Read 1 FASTQ file should be created per index per sequence run by combining the FASTQ files for each flow cell lane using a standard “cat” command after BCL conversion or by using the no-lane-splitting option when running bcl2fastq.
* Configuration file (conf.yml)
  * Global and system specific parameters (such as e.g. paths and genome assemblies installed) that are likely to remain constant across different runs of the pipeline are defined in the global configuration file (conf.yml). In each run the pipeline initially loads the parameters defined in this global configuration file, and then proceeds to load run specific parameters and the experiment specific data defined in a separate file (vpFile). The global configuration file can be edited using any standard text editor. 
  The list of parameters that need to be set at least once upon installation on a system in the global configuration file are shown below:
  
| Name            | Description                                                                                                                                            |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------|
| fragFolder      | Path to the folder containing the fragment end libraries of the reference genomes.                                                                     |
| normalizeFactor | Reads mapped to the 4C fragment end library are normalized to account for sequencing depth according to the normalizeFactor.                           |
| enzymes         | Enzyme names used in the viewpoint file and its corresponding recognition motif.                                                                       |
| genomes         | Genome names used in the viewpoint file plus corresponding BSgenome package.                                                                           |
| bowtie2         | Path to corresponding bowtie2 index of reference genome.                                                                                               |
| maxY            | Maximal Y value in local 4C cis plot.                                                                                                                  |
| plotView        | Number of bp to plot around viewpoint in local 4C cis plot.                                                                                            |
| xaxisUnit       | X-axis unit (Mb, Kb or bp).                                                                                                                            |
| plotType        | Plots will be saved as PDF or PNG                                                                                                                      |
| binSize         | Genome bin size used in the genome plot.                                                                                                               |
| qualityCutoff  | Q-score. Trim 3-end of all sequences using a sliding window as soon as 2 of 5 nucleotides has quality encoding less than the Q-score. 0 = no trimming. |
| trimLength     | Trim reads to defined capture length from 3-end. 0 = no trimming.                                                                                      |
| minAmountReads | Minimum amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed.                     |
| readsQuality   | Bowtie2 minimum quality mapped reads.                                                                                                                  |
| mapUnique     | Extract uniquely mapped reads, based on the lack of XS tag.                                                                                            |
| cores          | Number of cores for parallelization.                                                                                                                   |
| wSize          | The running mean window size.                                                                                                                          |
| nTop           | Top fragments discarded for normalization.                                                                                                             |
| nonBlind      | Only keep non-blind fragments.                                                                                                                         |
| wig           | Create wig files for all samples.                                                                                                                      |
| plot          | Create viewpoint coverage plot for all samples.                                                                                                        |
| genomePlot    | Create genomeplot for all samples (only possible if analysis is “all” in vpFile).                                                                      |
| tsv           | Create tab separated value file for all samples                                                                                                        |
| bins          | Count reads for binned regions.                                                                                                                        |
  
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
 

| **Name**         | **Description**                                                                                                                                                                                                                                                     |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| expname      | Unique experiment name.                                                                                                                                                                                                                                         |
| primer       | Primer sequence.                                                                                                                                                                                                                                                |
| firstenzyme  | First restriction enzyme name (nearest to reading primer)                                                                                                                                                                                                       |
| secondenzyme | Second restriction enzyme name.                                                                                                                                                                                                                                 |
| genome       | Reference genome of interest.                                                                                                                                                                                                                                   |
| vpchr        | The chromosome that contains the viewpoint (See note 13).                                                                                                                                                                                                       |
| vppos        | Coordinate of viewpoint position. Any position within the VP can be used except the RE motifs (see note 13).                                                                                                                                                    |
| analysis     | The final output tables will contain all reads (all) or only the reads that have been mapped to the VP chromosome (cis). For most analysis cis is sufficient and the generated output files will be smaller and therefore easier to process on local computers. |
| fastq        | Name of the FASTQ file.                                                                                                                                                                                                                                         |

**Table 3.** Description of parameters that are required in the viewpoint file for processing a 4C-seq experiment.




## Running the pipeline:

```
Rscript <path to pipe4C.R script> [any additional arguments]. 
```

A list of both required and optional parameters that are recognized by the pipe4C.R script are shown below. Default values (except vpFile, fqFolder, outFolder and confFile) are stored in the configuration file. 

For example,
```
Rscript pipe4C.R –-vpFile [path to vpFile] --fqFolder [path to folder containing the FASTQ files] –-outFolder [path to output folder] --cores 8 --wig --plot --genomePlot
```
will run the pipeline using 8 cores and generates a wig file, a viewpoint plot and a genome plot as output, next to the default outputs.


| Name            | Description                                                                                                                                        |
|-----------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| vpFile *       | path to the viewpoint file.                                                                                                                        |
| fqFolder *     | path to the folder containing the FASTQ files.                                                                                                     |
| outFolder *    | path to the output folder.                                                                                                                         |
| confFile       | path to configuration file – default is conf.yml in folder containing the pipeline script.                                                         |
| qualityCutoff  | Q-score. Trim 3-end of all sequences using a sliding window as soon as 2 of 5 nucleotides has quality encoding less than the Q-score. Default = 0. |
| trimLength     | Trim reads to defined capture length from 3-end. Default = 0 (no trimming).                                                                        |
| minAmountReads | Minimum amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed.                 |
| readsQuality   | Bowtie2 minimum quality mapped reads.                                                                                                              |
| mapUnique     | Extract uniquely mapped reads, based on the lack of XS tag.                                                                                        |
| cores          | Number of cores for parallelization.                                                                                                               |
| wSize          | The running mean window size.                                                                                                                      |
| nTop           | Top fragments discarded for normalization.                                                                                                         |
| nonBlind      | Only keep non-blind fragments.                                                                                                                     |
| wig           | Create wig files for all samples.                                                                                                                  |
| plot          | Create viewpoint coverage plot for all samples.                                                                                                    |
| genomeplot    | Create genomeplot for all samples (only possible if analysis is “all” in vpFile).                                                                  |
| tsv           | Create tab separated value file for all samples                                                                                                    |
| bins          | Count reads for binned regions.                                                                                                                    |

**Table 4.** Description of parameters that are recognized by the pipe4C.R script. * are required. 
