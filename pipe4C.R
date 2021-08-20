VERSION <- '1.1.4'
#1.1.4. introduced prefix
#1.1.3. introduced bigwig option


#.libPaths( c( .libPaths(), "~/mnt/G_drive/Group_deLaat_bioinf/shared/4C-seq_pipeline/R/installed/") )




get_script_path <- function(path=NULL) {
  if( is.null( path ) ){
    cmdArgs = commandArgs(trailingOnly = FALSE)
    needle = "--file="
    match = grep(needle, cmdArgs)
    if (length(match) > 0) {
      # Rscript
      return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
      ls_vars = ls(sys.frames()[[1]])
      if ("fileName" %in% ls_vars) {
        # Source'd via RStudio
        return(normalizePath(sys.frames()[[1]]$fileName))
      } else {
        # Source'd via R console
        return(normalizePath(sys.frames()[[1]]$ofile))
      }
    }
  } else {
    return(path)
  }
}

#################################################################################################################
### PARSING THE INPUT ###########################################################################################
#################################################################################################################
if( !suppressMessages(require( "optparse", character.only = TRUE ) ) ) stop( "Package not found: optparse" )

option_list = list(
  make_option(c("-v", "--vpFile"), type="character", default=NULL, 
              help="path to viewpoint file [required]", metavar="/path/to/vp_file"),
  make_option(c("-f", "--fqFolder"), type="character", default=NULL, 
              help="path to folder containing the FASTQ files [required]", metavar="/path/to/FASTQ_folder/"),
  make_option(c("-o", "--outFolder"), type="character", default=NULL, 
              help="path to the output folder [required]", metavar="/path/to/output_folder/"),
  make_option(c("-c", "--confFile"), type="character", default="conf.yml", 
              help="path to configuration file [default %default]", metavar="/path/to/conf.yml/"),
  
  
  
  make_option(c("-d", "--mismatchMax"), type="integer", default=NULL,
              help="Maximum number of mismatches allowed with primer sequence during demultiplexing.",metavar="number"),
  
  make_option(c("-q", "--qualityCutoff"), type="integer", default=NULL,
              help="Q-score. Trim 3-end of all sequences using a sliding window as soon as 2 of 5 nucleotides has quality encoding less than the Q-score.",metavar="number"),
  make_option(c("-l", "--trimLength"), type="integer", default=NULL,
              help="Trim reads to defined capture length from 3-end.",metavar="number"),
  make_option(c("-m", "--minAmountReads"), type="integer", default=NULL,
              help="Minimum amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed.",metavar="number"),
  make_option(c("-r", "--readsQuality"), type="integer", default=NULL,
              help="Bowtie2 minimum quality mapped reads.",metavar="number"),
  make_option(c("-z", "--cores"), type="integer", default=NULL,
              help="Number of cores for parallelization.",metavar="number"),
  
  make_option(c("-a", "--alnStart"), action="store_true", default=FALSE,
              help="Alignment should start at first position RE1."),
  
  
  make_option(c("-s", "--wSize"), type="integer", default=NULL,
              help="The running mean window size.",metavar="number"),
  make_option(c("-n", "--nTop"), type="integer", default=NULL,
              help="Top fragments discarded for normalization.",metavar="number"),
  
  
  make_option(c("-u", "--mapUnique"), action="store_true", default=FALSE,
              help="Extract uniquely mapped reads, based on the lack of XS tag."),
  make_option(c("-i","--nonBlind"), action="store_true", default=FALSE,
              help="Only keep non-blind fragments"),
  make_option(c("-w", "--wig"), action="store_true", default=FALSE,
              help="create wig files for all samples"),
  make_option(c("-x", "--bigwig"), action="store_true", default=FALSE,
              help="create big wig files for all samples"),
  make_option(c("-p", "--plot"), action="store_true", default=FALSE,
              help="Create viewpoint coverage plot for all samples."),
  make_option(c("-g", "--genomePlot"), action="store_true", default=FALSE,
              help="Create genomeplot for all samples (only possible if analysis is all in vpFile)."),
  make_option(c("-t", "--tsv"), action="store_true", default=FALSE,
              help="Create tab seperated values file for all samples"),
  make_option(c("-b", "--bins"), action="store_true", default=FALSE,
              help="Corunt reads for binned regions.")
  
)

helptext<-"Values stored in the configuration file (conf.yml) are used by default."

opt_parser = OptionParser(option_list=option_list, description=helptext)
argsL <- parse_args(opt_parser)

if (is.null(argsL$vpFile)){
  print_help(opt_parser)
  stop("vpFile is required.\n", call.=FALSE)
}
if (is.null(argsL$fqFolder)){
  print_help(opt_parser)
  stop("fqFolder is required.\n", call.=FALSE)
}
if (is.null(argsL$outFolder)){
  print_help(opt_parser)
  stop("outFolder is required.\n", call.=FALSE)
}


#################################################################################################################
### Load packages ###############################################################################################
#################################################################################################################
message( '\n------ Loading functions and configuration file' )



if( !suppressMessages(require( "ShortRead", character.only=TRUE ) ) ) stop( "Package not found: ShortRead" )
if( !suppressMessages(require( "GenomicRanges", character.only=TRUE ) ) ) stop( "Package not found: GenomicRanges" )
if( !suppressMessages(require( "GenomicAlignments", character.only=TRUE ) ) ) stop( "Package not found: GenomicAlignments" )
if( !suppressMessages(require( "caTools", character.only=TRUE ) ) ) stop( "Package not found: caTools" )
if( !suppressMessages(require( "config", character.only=TRUE ) ) ) stop( "Package not found: config" )


source( sub( pattern='pipe4C\\.R', replacement='functions\\.R', x=get_script_path()) )

if (argsL$confFile=='conf.yml'){
  argsL$confFile<-sub( pattern='pipe4C\\.R', replacement='conf.yml', x=get_script_path()) 
}

configOpt <- createConfig( confFile=argsL$confFile )
configOpt$pipeline.version <- VERSION


if (!is.null(argsL$qualityCutoff)){
  configOpt$qualityCutoff<-argsL$qualityCutoff
}
if (!is.null(argsL$trimLength)){
  configOpt$trimLength<-argsL$trimLength
}
if (!is.null(argsL$minAmountReads)){
  configOpt$minAmountReads<-argsL$minAmountReads
}
if (!is.null(argsL$readsQuality)){
  configOpt$readsQuality<-argsL$readsQuality
}
if (!is.null(argsL$cores)){
  configOpt$cores<-argsL$cores
}
if (!is.null(argsL$wSize)){
  configOpt$wSize<-argsL$wSize
}
if (!is.null(argsL$nTop)){
  configOpt$nTop<-argsL$nTop
}
if (!is.null(argsL$mismatchMax)){
  configOpt$mmMax<-argsL$mismatchMax
}

if (argsL$mapUnique){
  configOpt$mapUnique<-argsL$mapUnique
}
if (argsL$nonBlind){
  configOpt$nonBlind<-argsL$nonBlind
}
if (argsL$wig){
  configOpt$wig<-argsL$wig
}
if (argsL$bigwig){
  message("Making bigwigs!!")
  configOpt$bigwig<-argsL$bigwig
}
if (argsL$plot){
  configOpt$cisplot<-argsL$plot
}
if (argsL$genomePlot){
  configOpt$genomePlot<-argsL$genomePlot
}
if (argsL$tsv){
  configOpt$tsv<-argsL$tsv
}
if (argsL$bins){
  configOpt$bins<-argsL$bins
}
if (argsL$alnStart){
  configOpt$alnStart<-argsL$alnStart
}


#Check whether Bowtie2 and samtools are installed..







#################################################################################################################
### PIPELINE ####################################################################################################
#################################################################################################################


embryo()

Run.4Cpipeline( 
  VPinfo.file = argsL$vpFile
  ,FASTQ.F = argsL$fqFolder
  ,OUTPUT.F = argsL$outFolder
  ,configuration=configOpt
)


#################################################################################################################
### END #########################################################################################################
#################################################################################################################
