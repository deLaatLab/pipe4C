# Step by step tutorial

In this tutorial we will explain how to run the pipeline and perform peakC analysis on the Sox2 4C-seq data from Geeven et al. (SRA files GSM2824300, GSM2824301 and GSM2824302). 

## Set up the pipeline

#### Download the pipeline including the example files
Download the pipeline and example file using the following command. 

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
