# FUNCTIONS 4C PIPELINE


#2020.02.10: included option to ignore the chr_fix chromosomes that are currently present in the newest hg38 bsgenome package
#2020.02.05: Changed span from 100 to 1 in the exportWig function.
#option to generate BW files.
#2021.03.16: If no non-blind fragments are identified in a chr the digestion fucntion will give a warning and continue.
#2021.05.06 : If only 1 RE motif is found the digest function will not crash. The start en end of the chromosome will be use to generate 2 fragments.
#2021.05.06 :alignToFragends update:
  #Option that aligned read should start at first position RE1
  #Error when no algined fragments are identified
#2021.08.19 : add prefix for chrom names. by default chr
#2022.02.04: allow primer sequences to be only RE1 motif. This allow the use of fastq files for which the primer has been trimmed of, such as GEO uploads from other groups.
#2022.04.26:Allow IUPAC ambiquity code in RE motif
#2023 updated bigwig function to keep non std chroms

createConfig <- function( confFile=argsL$confFile ){
  message('reading config file')
  
  configF <- config::get(file=confFile )
  baseFolder <- configF$fragFolder
  normFactor <- configF$normalizeFactor
  plotView <- configF$plot$plotView
  maxY <- configF$plot$maxY
  xaxisUnit <- configF$plot$xaxisUnit
  plotType <- configF$plot$plotType
  binSize <- configF$plot$binSize
  mmMax <- configF$mismatchMax
  maxAmountReads <-configF$maxAmountReads
  qualityCutoff <- configF$qualityCutoff
  trimLength <- configF$trimLength
  minAmountReads <- configF$minAmountReads
  readsQuality <- configF$readsQuality
  cores <- configF$cores
  wSize <- configF$wSize
  nTop <- configF$nTop
  mapUnique <- configF$mapUnique
  nonBlind <- configF$nonBlind
  wig <- configF$wig
  bigwig <- configF$bigwig
  cisplot <- configF$cisplot
  genomePlot <- configF$genomePlot
  tsv <- configF$tsv
  bins <- configF$bins
  chr_random <- configF$chr_random
  chrUn <- configF$chrUn
  chrM <- configF$chrM
  chr_fix <- configF$chr_fix
  alnStart <- configF$alnStart
  prefix <- configF$prefix
  
  #Add forward slash to baseF if missing
  if(!substr(baseFolder, nchar(baseFolder), nchar(baseFolder))=='/'){
    baseFolder<-paste0(baseFolder,'/')
  }
    
  
  

  enzymes <- data.frame( 
    name=as.character( sapply( configF$enzymes, function(x) strsplit( x, split=' ' )[[1]][1] ) )
    , RE.seq=as.character( sapply( configF$enzymes, function(x) strsplit( x, split=' ' )[[1]][2] ) )
    , stringsAsFactors=FALSE 
    #, row.names=1 #This does not work if there is only 1 row.
  )
  
  row.names(enzymes) <- enzymes$name
  enzymes[1] <- NULL
  
  
  genomes <- data.frame( 
    genome=as.character( sapply( configF$genomes, function(x) strsplit( x, split=' ' )[[1]][1] ) )
    , BSgenome=as.character( sapply( configF$genomes, function(x) strsplit( x, split=' ' )[[1]][2] ) )
    , stringsAsFactors=FALSE 
    #, row.names=1
  )
  
  row.names(genomes) <- genomes$genome
  genomes[1] <- NULL
  
  
  bt2Genomes <- data.frame( 
    genome=as.character( sapply( configF$bowtie2, function(x) strsplit( x, split=' ' )[[1]][1] ) )
    , path=as.character( sapply( configF$bowtie2, function(x) strsplit( x, split=' ' )[[1]][2] ) )
    , stringsAsFactors=FALSE 
    #, row.names=1
  )
  
  row.names(bt2Genomes) <- bt2Genomes$genome
  bt2Genomes[1] <- NULL
  
  return( list( baseFolder=baseFolder
                ,normFactor=normFactor
                ,enzymes=enzymes
                ,genomes=genomes
                ,bt2Genomes=bt2Genomes
                ,plotView=plotView
                ,maxY=maxY
                ,xaxisUnit=xaxisUnit
                ,plotType=plotType
                ,binSize=binSize
                ,qualityCutoff=qualityCutoff
                ,trimLength=trimLength
                ,minAmountReads=minAmountReads
                ,maxAmountReads=maxAmountReads
                ,readsQuality=readsQuality
                ,cores=cores
                ,wSize=wSize
                ,nTop=nTop
                ,mapUnique=mapUnique
                ,nonBlind=nonBlind
                ,wig=wig
                ,bigwig=bigwig
                ,cisplot=cisplot
                ,genomePlot=genomePlot
                ,tsv=tsv
                ,bins=bins
                ,mmMax=mmMax
                ,chr_random=chr_random
                ,chr_fix=chr_fix
                ,chrUn=chrUn
                ,chrM=chrM
                ,alnStart=alnStart
                ,prefix=prefix
                
                
                
  ) )
}

Read.VPinfo<-function(VPinfo.file){
  #Indentify Fastq files, experiment names and primer sequence from VPinfo file
  VPinfo <- read.table(VPinfo.file, sep="\t", stringsAsFactors=FALSE, header=TRUE )
  
  if (ncol(VPinfo)==9){
    defaultNames<-c( 'expname', 'primer', 'firstenzyme', 'secondenzyme', 'genome', 'vpchr', 'vppos', 'analysis', 'fastq' )
    headerCount <- sum( colnames( VPinfo ) %in% defaultNames )
  }
  
  if (ncol(VPinfo)==10){
    defaultNames<-c( 'expname', 'spacer','primer', 'firstenzyme', 'secondenzyme', 'genome', 'vpchr', 'vppos', 'analysis', 'fastq' )
    headerCount <- sum( colnames( VPinfo ) %in% defaultNames ) 
  }
  
  if ( headerCount < ncol(VPinfo) ) {
    message( "Header not correct in VPinfo file." )
    return()
  }
  
  expname <- as.character(gsub("[^A-Za-z0-9_]", "", VPinfo$expname ))
  primer <-   toupper(as.character(gsub("[^A-Za-z]", "", VPinfo$primer )))
  firstenzyme <-  as.character(gsub("[^A-Za-z0-9]", "", VPinfo$firstenzyme ))
  secondenzyme <- as.character(gsub("[^A-Za-z0-9]", "", VPinfo$secondenzyme ))
  genome <- as.character(gsub("[^A-Za-z0-9_]", "", VPinfo$genome ))
  
  #vpchr <- as.character(gsub("[^0-9XYMxym]", "",VPinfo$vpchr ))
  
  #To be sure that any kind of chromosome name for different species can be used I only remove spaces.
  vpchr <- as.character(gsub(" ", "",VPinfo$vpchr ))
  #vpchr <- as.character(gsub("chr", "",vpchr ))
  #vpchr <- as.character(gsub("Chr", "",vpchr ))
  
  vppos <- as.numeric(gsub("\\D+", "", VPinfo$vppos )) #\D, matches any character that is not a decimal digit.
  analysis <- as.character(gsub("[^a-z]", "", VPinfo$analysis ))
  fastq <- as.character( VPinfo$fastq )
  
  if (ncol(VPinfo)==9){
    VPinfo$spacer<-0
  }
  
  spacer <- as.numeric(gsub("\\D+", "", VPinfo$spacer ))
  
  
  return( data.frame( expname=expname, spacer=spacer, primer=primer, firstenzyme=firstenzyme, secondenzyme=secondenzyme, genome=genome, vpchr=vpchr, vppos=vppos, analysis=analysis, fastq=fastq, stringsAsFactors=FALSE ) )         
}

demux.FASTQ <- function(VPinfo, FASTQ.F, FASTQ.demux.F, demux.log.path, overwrite = FALSE, mmMax = 0) {
  # reading functions
  if (!require("ShortRead", character.only = TRUE)) stop("Package not found: ShortRead")
  
  # Demultiplex FastQ files: 
  # 1.  This will remove weird charactes from VPinfo name and primer info.
  # 2.  Exclude experiments that have non-unique names
  # 3.  Demultiplex each Fastq file and save it as a new Fastq file which can be uploaded to GEO.
  #     If a spacer sequence is used with random nucleotids before the primer, the spacer option should be set.
  
  
  #Demultiplex per FASTQ file
  
  
  file.fastq <- sort(unique(VPinfo$fastq))
  for (i in 1:length(file.fastq)) {
    VPinfo.fastq <- VPinfo[VPinfo$fastq == file.fastq[i], ]
    
    # Remove weird characters from VPinfo file
    exp.name.all <- as.character(VPinfo$expname)
    exp.name <- as.character(VPinfo.fastq$expname)
    primer <- toupper(as.character(VPinfo.fastq$primer))
    spacer <- as.numeric(VPinfo.fastq$spacer)
    
    # Remove experiments with duplicated names
    exp.name.unique <- NULL
    primer.unique <- NULL
    spacer.unique <- NULL
    for (a in 1:length(exp.name)) {
      if (exp.name[a] %in% exp.name.all[duplicated(exp.name.all)]) {
        error.msg <- paste("      ### ERROR: Experiment name not unique for ", exp.name[a])
        message(error.msg)
        write(error.msg, demux.log.path, append = TRUE)
      } else {
        exp.name.unique <- c(exp.name.unique, exp.name[a])
        primer.unique <- c(primer.unique, primer[a])
        spacer.unique <- c(spacer.unique, spacer[a])
      }
    }
    fq.df <- data.frame()
    for (b in 1:length(exp.name.unique)) {
      outfile <- paste0(FASTQ.demux.F, exp.name.unique[b], ".fastq.gz")
      if (any(file.exists(outfile))) {
        if (overwrite == TRUE) {
          unlink(outfile)
          error.msg <- paste("      ### WARNING: File", outfile, "exists and will be overwritten.")
          write(error.msg, demux.log.path, append = TRUE)
          message(error.msg)
          newrow <- data.frame(exp.name = exp.name.unique[b], primer = primer.unique[b], spacer = spacer.unique[b])
          fq.df <- rbind(fq.df, newrow)
        } else {
          message(paste("      ### WARNING: File", outfile, "exists. continuing with exisiting file."))
        }
      } else {
        newrow <- data.frame(exp.name = exp.name.unique[b], primer = primer.unique[b], spacer = spacer.unique[b])
        fq.df <- rbind(fq.df, newrow)
      }
    }
    
    #Check whether primers hae enough distance with maximum allowed mismatches.
    #
    #
    #https://www.rdocumentation.org/packages/DescTools/versions/0.99.19/topics/StrDist
    #The function computes the Hamming and the Levenshtein (edit) distance of two given strings (sequences).
    #The Hamming distance between two vectors is the number mismatches between corresponding entries.
    #In case of the Hamming distance the two strings must have the same length.
    #
    
    if (nrow(fq.df)>1){
      message("Check whether primer can be seperated" )
      
      primers.mm<-DNAStringSet(fq.df$primer)
      names(primers.mm)<-fq.df$exp.name
      
      for (c in seq_along(primers.mm)){
        primers.mm2<-primers.mm[-c]
        #dist <- srdistance(primers.mm2,primers.mm[c])[[1]]
        dist <- srdistance(primers.mm2,primers.mm[c],method = "Hamming")[[1]]
        
        #srdistanceFilter()
        #http://manuals.bioinformatics.ucr.edu/home/ht-seq
        
        
        if (min(dist)<= mmMax){
          error.msg <- paste("      ### WARNING: primer sequence not unique for", fq.df$exp.name[c], "with", mmMax, "mismatches allowed.")
          message(error.msg)
          write(error.msg, demux.log.path, append = TRUE)
          error.msg2 <- paste("      ### WARNING: primer overlaps with primer ", names(primers.mm2[dist<=mmMax]),"\n")
          message(error.msg2)
          write(error.msg2, demux.log.path, append = TRUE)
        }
      }
    }
    
    # Read FastQ
    if (nrow(fq.df) > 0) {
      message(paste("      >>> Reading Fastq: ", file.fastq[i], " <<<"))
      
      
      if (file.exists(paste0(FASTQ.F, "/", file.fastq[i]))) {
      
      # then we stream the fastq files at 1,000,000 reads each time ( default )
      stream <- FastqStreamer((paste0(FASTQ.F, "/", file.fastq[i])))
      message(paste0("         ### ", fq.df$exp.name, " check for primer sequence ", fq.df$primer, 
                     " spacer:", fq.df$spacer, " max mismatch:", mmMax, "\n"))
      while (length(fq <- yield(stream))) {
        for (i in 1:nrow(fq.df)) {
          primer.seq <- as.character(fq.df$primer[i])
          spacer <- as.numeric(fq.df$spacer[i])
          if (mmMax == 0) {
            demultiplex.primer = srFilter(function(x) {
              substr(sread(x), spacer + 1, spacer + nchar(primer.seq)) == primer.seq
            }, name = "demultiplex.primer")
            demux.fq <- fq[demultiplex.primer(fq)]
          } else {
            shortreads <- narrow(sread(fq), start = spacer + 1, end = spacer + nchar(primer.seq))
            dist <- srdistance(shortreads, primer.seq, method = "Hamming")[[1]]
            demux.fq <- fq[dist <= mmMax]
          }
          writeFastq(demux.fq, paste0(FASTQ.demux.F, fq.df$exp.name[i], ".fastq.gz"), mode = "a")
        }
      }
      close(stream)
    
    }else{
      error.msg <- paste("      ### WARNING: File", file.fastq[i], " does not exist. Experiment skipped.")
      write(error.msg, demux.log.path, append = TRUE)
      message(error.msg)
    }
    }
  }
}


trim.FASTQ <- function( exp.name, primer, firstcutter, secondcutter, file.fastq, trim.F, cutoff, trim.length=0, log.path, min.amount.reads=1000, maxAmountReads=1e8){
  txt.tmp <- paste0( trim.F, exp.name, ".txt" )
  info.file <- paste0( trim.F, exp.name, ".info.rds" )
  if ( file.exists( txt.tmp ) & file.exists( info.file ) ) {
    error.msg <- paste( "         ### WARNING: trimmed file", exp.name, "already exists, continuing with exisiting file." )
    write( error.msg, log.path, append=TRUE )
    message( error.msg )
    return( readRDS(info.file) ) 
  } else {
    message( paste( "         ### Reading Fastq: ", file.fastq ))
    message( paste( "         ### Reading max Reads: ", maxAmountReads ))
    demux.fq <- readFastq( file.fastq, n=maxAmountReads)
    
    # Qualtity trimming
    if ( cutoff > 0 ) {
      # Trim ends with qualtiy score < cutoff
      #trimTailw remove low-quality reads from the right end using a sliding window
      #trim as soon as 2 of 5 nucleotides has quality encoding less than phred score.
      Phred.cutoff <- rawToChar(as.raw(cutoff + 33))
      demux.fq <- trimTailw( demux.fq, 2, Phred.cutoff, 2)
      # Alternatively use trimTails
      # remove a tally of (successive) nucleotides falling at or below a quality threshold (trimTails).
      # fq.trimTails.10 <- trimTails(fq, k=2, a=Phred.cutoff, successive=FALSE)
      # successive: indicating whether failures can occur anywhere in the sequence, or must be successive
      #Mean average base quality score>cutoff
      #Does this make sense after end trimming? Should the cutoff be higher?
      demux.fq <- demux.fq[(alphabetScore(demux.fq) / width(demux.fq)) > cutoff]
      #Drop reads that are less than 30nt after quality filtering. Discuss lenght with Geert/Valerio
      demux.fq <- demux.fq[width(demux.fq) >= 30]
      #Maybe do this filtering afterwards? read len=Cap length?? (or a few nt less?)
    } #close cutoff
    

    
    sequences <- sread( demux.fq )
    nReads <- length( sequences )
    
    
    
    if ( nReads < min.amount.reads) {
      error.msg <- paste0( "         ### ERROR: ",exp.name," - Less reads in FASTQ than set as minimum. Reads: ", nReads )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      message("         ### To continue alter minAmountRead argument")
      return()
    } else {
      message( paste0( "         ### Total Reads: ", nReads ) )
    }
    
    if ( nReads > 2e6) {
    message("More than 2 million reads sequences. This is for normal 4C experiments is not required.")  
    message("Depending on the memory available in your machine the pipeline may crash. Will fix this in the future.")
    message("For now a quick fix is to alter the maxAmountReads parameter in the config file. Reduce this number.")
    }
    

    
    
    #Find the most occuring position of the firstcutter
    motifPos<-sort( table( regexpr( firstcutter, sequences ) ), decreasing=TRUE)[1]
    motif.1st.pos<-as.numeric( names(motifPos))
    motifPos.perc<-round(as.numeric(motifPos/nReads*100),2)
    
    
    if ( trim.length > 0 ){
      sequences <- substr( sequences, 1, ( trim.length-1+motif.1st.pos ) )
      
      read.length.table <- sort(table(width(sequences)), decreasing=TRUE)[1]
      read.length<-as.numeric(names(read.length.table))
      read.length.perc<-round(as.numeric(read.length.table/nReads*100),2)
      
    }else{
      read.length.table <- sort(table(width(demux.fq)), decreasing=TRUE)[1]
      read.length<-as.numeric(names(read.length.table))
      read.length.perc<-round(as.numeric(read.length.table/nReads*100),2)
    }
    
    
    if ( read.length.perc < 60) {
      error.msg <- paste0( "         ### WARNING:", exp.name," - Max read length in only ", read.length.perc, "% of the reads." )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
    }
    

   
    
   
    
    
    motif.1st.pos.2nd <- FALSE
    if ( motif.1st.pos > 0 ){
      #Check whether firstcutter motif is within the first 4 nts (as part of a barcode). If so take 2nd firstcutter pos.
      #This is a problem for fastq files for which the primer has been trimmed (GEO fastq files from other groups). Check whether primer sequence is equal to RE motif
      
      
      if ( !primer==firstcutter & motif.1st.pos < 5 ){
        message("RE1 motif identified in first 4 nts, most likely used as barcode. Use 2nd RE1 motif identified in reads")
        
        motif.1st.pos.2nd <- TRUE
        
        #motif.1st.pos <- as.numeric( names( sort( table( gregexpr( firstcutter, sequences )[[1]][2] ), decreasing=TRUE ) )[1] )
        
        motifPos<-sort( table( gregexpr( firstcutter, sequences )[[1]][2] ), decreasing=TRUE )[1]
        motif.1st.pos<-as.numeric( names(motifPos))
        motifPos.perc<-round(as.numeric(motifPos/nReads*100),2)
        
        
      }
      
      if ( motifPos.perc < 90) {
        error.msg <- paste0( "         ### WARNING: ", exp.name, " - Most occuring position RE1 found in only ", motifPos.perc, "% of the reads." )
        write( error.msg, log.path, append=TRUE )
        message( error.msg )
      }
      
      
      
      
      
     
      
      
      # Trim non-blind fragend sequences based on 2nd cutter 
      motif.2ndRE.pos <- as.numeric( names( sort( table( regexpr( secondcutter, sequences ) ), decreasing=TRUE ) )[1] ) #Determine whether part of barcode
      if ( motif.2ndRE.pos > 5 | motif.2ndRE.pos == -1 ){
        sequences.no.2nd.enzyme <- sub( paste0( secondcutter, ".*$" ), secondcutter, sequences ) # Changes only the 1st pattern match per string
      } else { #part of bardcode. use the 2nd motif
        sequences.no.2nd.enzyme <- sub( paste0( "(", secondcutter, ".*?)", secondcutter, ".*$" ), paste0( "\\1", secondcutter ), sequences )
      }
      # The . matches any character. .* matches zero or more instances of any character. 
      # But the matching is greedy; it would match as much as possible. I want matching that is not greedy (only up until the second .) 
      # so I added ? to suppress the greedy match and .*? matches any group of characters up until we hit the next thing in the regex.
      # Because the first part was enclosed in parentheses (\\..*?) it is stored as \1,
      # so the substitution pattern \\1 restores everything before the second . and the second . is replaced with the secondcutter.
      # Trim 2nd first cutter motif ( Trim blind fragments )
      if ( motif.1st.pos.2nd == FALSE ){
        sequences.no.2nd.enzyme <-sub( paste0( "(", firstcutter, ".*?)", firstcutter, ".*$" ), paste0( "\\1", firstcutter ), sequences.no.2nd.enzyme )      
      }else{
        # Trim 3rd first cutter motif (Trim blind fragments) if present in barcode
        sequences.no.2nd.enzyme <-sub( paste0( "(", firstcutter, ".*?", firstcutter, ".*?)", firstcutter, ".*$"), paste0( "\\1", firstcutter ), sequences.no.2nd.enzyme )
      }
      # Trim sequences based on 1st cutter
      keep <- substr( sequences.no.2nd.enzyme, motif.1st.pos, nchar( sequences.no.2nd.enzyme ) )
    } # close if (motif.1st.pos>0)
    
    if( motif.1st.pos == -1 ){
      rm( sequences, demux.fq )
      error.msg <- paste( "         ### ERROR:", exp.name, "1st RE motif not found in reads" )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      return()
    } 
    write( keep, txt.tmp )
    rm( sequences, demux.fq )
    captureLen <- read.length - motif.1st.pos + 1
    
    #If there are many primer dimers than this goes wrong....Maybe just used maximum read length?
    
    saveRDS( list( captureLen=captureLen, nReads=nReads, motifPosperc=motifPos.perc, readlenperc=read.length.perc  ), file=info.file )
    return( list( captureLen=captureLen, nReads=nReads, motifPosperc=motifPos.perc, readlenperc=read.length.perc ) )
  }
}


makeBAM <- function( exp.name, BAM.F, NCORES, Bowtie2Folder, genome, txt.tmp, log.path, bowtie.log.path, bamFile, map.unique, readsQual=1 ) {
  
  #Only unique, use grep and samtools or use R for this?
  
  TEMPfile <- tempfile( pattern = "aln.", tmpdir = tempdir(), fileext = "" )
  check.trunc <- 1
  while ( check.trunc < 10 ) {
    message( paste0( "         ### attempt ", check.trunc ) )
    if ( map.unique == TRUE ) {
      bamFile.tmp <- paste0( BAM.F, exp.name, "_tmp.bam" )
      CMD <- paste0(
        "(bowtie2 -p ",
        NCORES,
        " -r -x ",
        Bowtie2Folder,
        genome,
        " -U ",
        txt.tmp,
        " | samtools view -hbSu - | samtools sort -T ", TEMPfile, " -o ",
        bamFile.tmp, ") 2>&1"
      )
      bowtie.output <- system(
        command=CMD
        , intern=TRUE
      )
      for( i in bowtie.output ){ message( paste0( '          ', i) ) }
      write(paste("Bowtie2 output:",exp.name), bowtie.log.path, append=TRUE) 
      write(bowtie.output, bowtie.log.path, append=TRUE)    
      #Check whether BAM in truncated (Due to pipe problems?). If truncated try 5 times.
      if (length(grep("truncated file", bowtie.output)) == 1) {
        error.msg <- paste("         ### ERROR: ",exp.name," - Truncated BAM file. Repeating Bowtie mapping...")
        message(error.msg)
        write(error.msg, log.path, append=TRUE)
        check.trunc <- check.trunc + 1
        unlink(bamFile.tmp)
        
      } else {
        check.trunc <- 10
        system(paste0("samtools view -H ", bamFile.tmp, " > header.sam")) #Output the header only
        system(
          paste0(
            "samtools view -F 4 ",bamFile.tmp,
            " | grep -v \"XS:\" | cat header.sam - | samtools view -b - > ",
            bamFile
          )
        )
        system(
          paste0(
            "samtools index ",
            bamFile
          )
        )
        unlink(bamFile.tmp)
      }
    } else{
      # map with quality, default 1
      CMD <- paste0(
        "(bowtie2 -p ",
        NCORES,
        " -r -x ",
        Bowtie2Folder,
        genome,
        " -U ",
        txt.tmp,
        " | samtools view -q ", readsQual, " -bSu - | samtools sort -T ", TEMPfile, " -o ",
        bamFile, ") 2>&1"
      )
      bowtie.output <- system(
        command=CMD
        , intern=TRUE
      )
      system(
        paste0(
          "samtools index ",
          bamFile
        )
      )
      for( i in bowtie.output ){ message( paste0( '          ', i) ) }
      write(paste("Bowtie2 output:",exp.name), bowtie.log.path, append=TRUE)
      write(bowtie.output, bowtie.log.path, append=TRUE)
      #Check whether BAM in truncated. If truncated try 10 times.
      if (length(grep("truncated file", bowtie.output)) == 1) {
        error.msg <-paste0("         ### ERROR: ",exp.name," - Truncated BAM file. Repeating Bowtie mapping...")
        message(error.msg)
        write(error.msg, log.path, append=TRUE)
        check.trunc <- check.trunc + 1
        unlink(bamFile)
        
      } else{
        check.trunc <- 10
      }
    }
  }
  unlink( TEMPfile )
}

norm4C <- function( readsGR, nReads=1e6, nTop=2, wSize=21 ) {
  readsGR$normReads <- 0
  sumTop <- sum( -sort( -readsGR$reads )[ 1:nTop ] )
  wNorm <- nReads/( sum( readsGR$reads )-sumTop )
  readsGR$normReads <- wNorm*readsGR$reads
  readsGR$norm4C <- runmean( x=readsGR$normReads, k=wSize, endrule="mean" )
  return( readsGR )
}

exportWig <- function( gR, expName, filename, vpPos, vpChr, plotView) {
  gzname<- paste0( filename, ".gz" )
  gz1 <- gzfile( gzname, "w" )
  browserPos <- paste0( vpChr,":", vpPos-plotView, "-", vpPos+plotView )
  positionLine <- paste0( "browser position ", browserPos, "\n" )
  cat( positionLine, file=gz1 )
  trackLine <- paste0( "track type=wiggle_0 name=", expName, " visibility=full autoScale=off  viewLimits=0.0:2500.0 maxHeightPixels=50:50:8 color=0,0,200  priority=10\n" )
  cat(trackLine, file=gz1, append=TRUE )
  chrs <- as.vector( unique( seqnames( gR ) ) )
  for( chr in chrs ) {
    #chrLine <- paste0( "variableStep chrom=", chr, " span=100\n" )
    #wigToBigWig requires fixed bins.
    chrLine <- paste0( "variableStep chrom=", chr, " span=1\n" )
    cat( chrLine, file=gz1, append=TRUE )
    chrGr <- gR[seqnames(gR)==chr]
    datChr <- data.frame( start=chrGr$pos, score=round(chrGr$norm4C, digits=5 ) )
    write.table( datChr, file=gz1, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE )
  }
  close(gz1)
}

olRanges <- function(query, subject) {
 
  ## Find overlapping ranges
  olindex <- as.matrix(findOverlaps(query, subject))
  query <- query[olindex[,1]]
  subject <- subject[olindex[,2]]
  olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))
  
  ## Pre-queries for overlaps
  startup <- olma[,"Sstart"] < olma[,"Qstart"]
  enddown <- olma[,"Send"] > olma[,"Qend"]
  startin <- olma[,"Sstart"] >= olma[,"Qstart"] & olma[,"Sstart"] <= olma[,"Qend"]
  endin <- olma[,"Send"] >= olma[,"Qstart"] & olma[,"Send"] <=  olma[,"Qend"]
  
  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin
  contained <- startup & enddown
  
  ## Overlap types in one vector
  OLtype <- rep("", length(olma[,"Qstart"]))
  OLtype[olup] <- "olup"
  OLtype[oldown] <- "oldown"
  OLtype[inside] <- "inside"
  OLtype[contained] <- "contained"
  
  ## Overlap positions
  OLstart <- rep(0, length(olma[,"Qstart"]))
  OLend <- rep(0, length(olma[,"Qstart"]))
  OLstart[olup] <- olma[,"Qstart"][olup]
  OLend[olup] <- olma[,"Send"][olup]
  OLstart[oldown] <- olma[,"Sstart"][oldown]
  OLend[oldown] <- olma[,"Qend"][oldown]
  OLstart[inside] <- olma[,"Sstart"][inside]
  OLend[inside] <- olma[,"Send"][inside]
  OLstart[contained] <- olma[,"Qstart"][contained]
  OLend[contained] <- olma[,"Qend"][contained]
  
  ## Absolute and relative length of overlaps
  OLlength <- (OLend - OLstart) + 1
  OLpercQ <- OLlength/width(query)*100
  OLpercS <- OLlength/width(subject)*100
  
  ## Output type
  oldf <- data.frame(Qindex=olindex[,1], Sindex=olindex[,2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
  elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf[,-c(3,4)])
  return(query)
  
}

alignToFragends <- function( gAlign, fragments, firstcut, alnStart=TRUE ) {
  
  
  
  strand( fragments ) <- ifelse( fragments$fe_strand==3, "-", "+" )
  
  
  #To do: Remove olRanges function.
  #https://rdrr.io/bioc/IRanges/man/findOverlaps-methods.html
  #If type is start or end, the intervals are required to have matching starts or ends, 
  #minoverlap=nchar(firstcut)+1
  
  
  
  ovl <- olRanges( query=fragments, subject=as( gAlign, "GRanges" ) )
  #remove sequences that overlap only with RE motif
  ovl <- ovl[ ovl$OLlength >= nchar( as.character( firstcut ) )+1 ]
  
  
  if(alnStart){
    #Aligned read should start at the RE1 motif
    ovl <- ovl[ strand( ovl )=="+" & ovl$Sstart == start( ovl ) | strand( ovl )=="-" & ovl$Send == end( ovl ) ]
    
    #Or make it less strict? SNPs or read erros may result in a read that maps a few nts later in the fragment?
    #ovl <- ovl[ strand( ovl )=="+" & ovl$Sstart >= start( ovl ) & ovl$Sstart <= start( ovl )+nchar(as.character(firstcut))| strand( ovl )=="-" & ovl$Send <= end( ovl ) & & ovl$Send >= start( ovl )-nchar(as.character(firstcut)) ]
  }else{
    #Sequence need to start within unique fragment
    ovl <- ovl[ strand( ovl )=="+" & ovl$Sstart >= start( ovl ) | strand( ovl )=="-" & ovl$Send <= end( ovl ) ]
  }
  
  
  if(length(ovl)==0){
    error.msg <- paste( "         ### ERROR: no unique fragment ends aligned" )
    message( error.msg )
    return()
  }
  
  
  nReads <- tapply( ovl$Sindex, ovl$Qindex, length )
  
  
  
  fragments$reads <- 0
  fragments$reads[ as.numeric( names( nReads ) ) ] <- nReads
  return( fragments )
}



Digest <- function( assemblyName, firstcutter_Digest, secondcutter_Digest, baseFolder_Digest, config_genomes
                    ,chr_random,chr_fix,chrUn,chrM) {

  firstcutter <- as.character( firstcutter_Digest )
  secondcutter <- as.character( secondcutter_Digest )
  
  rdsFile <- paste0( baseFolder_Digest, assemblyName, "_", firstcutter, "_", secondcutter, ".rds" )
  
  if ( !file.exists( rdsFile ) ){
    
    if(!dir.exists( baseFolder_Digest )){
      dir.create( baseFolder_Digest )
    } 
    
    
    do.call( require, args=list( config_genomes[ assemblyName, ] ) )
    assign( 'frag.genome', base::get( config_genomes[ assemblyName, ] ) )
    
    chr <- unique( seqnames( frag.genome ) )
    
    #The "hap" ones are alternate assemblies for certain regions.
    #!! Do NOT use the hap and alt files !!!!
    chr <- chr[ grep( pattern="_hap", x=chr, invert=TRUE ) ]
    chr <- chr[ grep( pattern="_alt", x=chr, invert=TRUE ) ]
    
    #Make sure Bowtie2 index does not contain these chr
    if (chr_random==FALSE){
      chr <- chr[ grep( pattern="_random", x=chr, invert=TRUE ) ]
    }
    if (chr_fix==FALSE){
      chr <- chr[ grep( pattern="_fix", x=chr, invert=TRUE ) ]
    }
    
    if (chrUn==FALSE){
      chr <- chr[ grep( pattern="chrUn_", x=chr, invert=TRUE ) ]
    }
    if (chrM==FALSE){
      chr <- chr[ grep( pattern="chrM", x=chr, invert=TRUE ) ]  
    }
    
    
    
    outFrags <- GRanges()
    for ( i in seq_along( chr ) ){
      chrom <- chr[ i ]
      print( paste( "Processing ", chrom, sep="" ) )
      
      
      RE1_pos <- matchPattern( pattern=firstcutter, subject=frag.genome[[chrom]], fixed="subject" )
      
      if (length(RE1_pos)>0){
        
        RE1 <- GRanges( seqnames=chrom, ranges( RE1_pos ) )
        
        if(length(RE1_pos)>1){
          frag.RE1 <- RE1[ 1:( length(RE1)-1 ) ]
          end( frag.RE1 )[ 1:( length(RE1)-1 ) ] <- end( RE1 )[ 2:length( RE1 ) ]
        }else{
          message("Only 1 RE1 motif found")
          chromLen<-as.numeric(as.data.frame(seqlengths(frag.genome))[chrom,1])
          RE1_A<-RE1
          start(RE1_A)<-1
          RE1_B<-RE1
          end(RE1_B)<-chromLen
          frag.RE1 <- c(RE1_A,RE1_B)
        }
        
        

        #Allow the usage of N
        #https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/matchPattern
        #IUPAC ambiguity code in the pattern can match any letter in the subject that is associated with the code, and vice versa. 
        #RE2 <- matchPattern( pattern="CTNAG", subject=BSgenome.Hsapiens.UCSC.hg38[["chr1"]], fixed="subject" )
        
        RE2 <- matchPattern( pattern=secondcutter, subject=frag.genome[[chrom]], fixed="subject" )
        

        
        
        if(length(RE2)>0){
          RE2 <- GRanges( seqnames=chrom, ranges( RE2 ) )
          
          #This is wrong
          #RE2.ol <- olRanges(frag.RE1, RE2)
          #RE2.ol <- RE2[RE2.ol[RE2.ol$OLpercS==100]$Sindex]
          #blinds <- frag.RE1[!(frag.RE1 %over% RE2.ol)]
          #nonBlinds <- frag.RE1[unique(findOverlaps(frag.RE1,RE2.ol)@from)]
        
          
          #For blind fragments only keep RE2 sites that completely overlap with the RE1 fragment.
            
          blinds <- frag.RE1[
            -(unique(
              subjectHits(
                findOverlaps(
                  query=RE2,
                  subject=frag.RE1,
                  type="within",
                  ignore.strand=FALSE)
              )
            ))
            ]
          
          nonBlinds <- frag.RE1[unique(subjectHits(findOverlaps(query=RE2,
                                                                subject=frag.RE1,
                                                                type="within",
                                                                ignore.strand=FALSE)))]
          
          
          
          
          
        }else{
          message("No RE2 motifs found")
          blinds <- frag.RE1
          nonBlinds <- GRanges()
        }
        
        
        
        
        if(length(blinds)>0){
          blinds$type <- "blind"
          
          if(length(RE1_pos)>1){
            blinds <- rep(blinds[blinds$type=="blind"], each = 2)
          }
          
          blinds$fe_strand <- rep(c(5,3), length.out = length(blinds))
        }
        
        
        
        if(length(nonBlinds)>0){
          nonBlinds$type <- "non_blind"
          nonBlinds.fe5 <- nonBlinds
          
          #end( nonBlinds.fe5 ) <- end( RE2.ol[ findOverlaps( nonBlinds, RE2.ol, select="first" ) ] )
          
          end(nonBlinds.fe5)<-end(RE2[findOverlaps(query=nonBlinds,
                                                   subject=RE2,
                                                   minoverlap=nchar(secondcutter),
                                                   select = "first",
                                                   ignore.strand=FALSE)])
          
          
          if (length(RE1_pos)==1){
            
            #Only keep the fragend downstream of the RE1

            nonBlinds.fe5<-nonBlinds.fe5[2]
          }
          
          nonBlinds.fe5$fe_strand <- 5
          
          
          nonBlinds.fe3 <- nonBlinds
          #start( nonBlinds.fe3 ) <- start( RE2.ol[ findOverlaps( nonBlinds, RE2.ol, select="last" ) ] )
          start(nonBlinds.fe3) <-start(RE2[findOverlaps(query=nonBlinds,
                                                        subject=RE2,
                                                        minoverlap=nchar(secondcutter),
                                                        select = "last",
                                                        ignore.strand=FALSE)])
          
          
          if (length(RE1_pos)==1){
            nonBlinds.fe3<-nonBlinds.fe3[1]
          }
          
          nonBlinds.fe3$fe_strand <- 3
          
          nonBlinds.fe5.start <- GRanges()
          nonBlinds.fe3.end <- GRanges()
          
          if (length(RE1_pos)>1){
            #Add first and last fragends
            
            # if ( start( head( RE1, 1 ) ) > start( head( RE2, 1 ) ) ){
            #   frag <- GRanges( seqnames=chrom, IRanges( 1, end(RE1[1]) ) )
            #   RE2.ol.start <- olRanges(frag, RE2)
            #   RE2.ol.start <- RE2[RE2.ol.start[RE2.ol.start$OLpercS==100]$Sindex]
            #   nonBlinds.fe5.start <- RE1[1]
            #   start(nonBlinds.fe5.start) <- start(RE2.ol.start[findOverlaps(frag, RE2.ol.start,select="last")])
            #   nonBlinds.fe5.start$type <- "non_blind"
            #   nonBlinds.fe5.start$fe_strand <- 5
            # }
            
            if (start(RE1)[1] > start(RE2[1])) {
              
              #-----------RE2------RE1----------------------end
              #first frament end will be added
              
              frag <- GRanges(seqnames = chrom, IRanges(1, end(RE1[1])))
              nonBlinds.fe5.start <- RE1[1]
              
              start(nonBlinds.fe5.start) <-start(RE2[findOverlaps(query=frag,
                                                                  subject=RE2,
                                                                  minoverlap=nchar(secondcutter),
                                                                  select = "last",
                                                                  ignore.strand=FALSE)])

              nonBlinds.fe5.start$type <- "non_blind"
              nonBlinds.fe5.start$fe_strand <- 5
             
            }
      
            
            
            
            
            if ( end( tail( RE1, 1 ) ) <  end( tail( RE2, 1 ) ) ){
               frag <- GRanges( seqnames=chrom, IRanges( start(RE1[length(RE1)]), end(RE2[length(RE2)]) ) )
            
            #   RE2.ol.start <- olRanges(frag, RE2)
            #   RE2.ol.start <- RE2[RE2.ol.start[RE2.ol.start$OLpercS==100]$Sindex]
            
                   nonBlinds.fe3.end <- RE1[length(RE1)]
                   end(nonBlinds.fe3.end) <-end(RE2[findOverlaps(query=frag,
                                                                 subject=RE2,
                                                                 minoverlap=nchar(secondcutter),
                                                                 select = "first",
                                                                 ignore.strand=FALSE)])
                   
              
            #   end(nonBlinds.fe3.end) <- end(RE2.ol.start[findOverlaps(frag, RE2.ol.start,select="first")])
               nonBlinds.fe3.end$type <- "non_blind"
               nonBlinds.fe3.end$fe_strand <- 3
             }
            
            
              
            
              
            
            
            
          }
          
          
          outFrags <- suppressWarnings(sort(c(outFrags,blinds, nonBlinds.fe5,nonBlinds.fe3,nonBlinds.fe5.start,nonBlinds.fe3.end)))
        }else{
          message(paste("No non-blind fragments in",chrom))
          outFrags <- suppressWarnings(sort(c(outFrags,blinds)))
        }
      }else{
        message("No first cutter motifs found")
      }
      
      
      
    }
    
    outFrags$fe_id <- 1:length( outFrags ) 
    outFrags$pos <- ifelse( outFrags$fe_strand == 5, start( outFrags )+nchar(firstcutter), end( outFrags )-nchar(firstcutter) ) 
    
    #Save new RDS
    message( 'Saving Frag genome file as: ',rdsFile )
    saveRDS( outFrags, file=rdsFile )
    return( outFrags )
  } else {
    message( 'Frag genome file already exists.' )
    return( 0 )
  }
}


getUniqueFragends <- function( fragsGR_Unique, firstcutter_Unique="GATC", secondcutter_Unique="GTAC", captureLen_Unique=50, nThreads_Unique=4, baseFolder_Unique, Bowtie2Folder, assemblyName="mm9", config_genomes ) {
  ##If captureLen is not in colnames, generate this column
  if ( paste0( "len", captureLen_Unique ) %in% colnames( elementMetadata( fragsGR_Unique ) ) == FALSE ){
    
    message(paste('          ### Create FASTA file for capture length', captureLen_Unique))
    
    fragends <- fragsGR_Unique
    refGenomeFile <- paste0( Bowtie2Folder, assemblyName )
    fragFile <- paste0( baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, ".rds" )
    fastaFile <- paste0( baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, "_", captureLen_Unique, ".fa" )
    bamFile <- paste0( baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, "_", captureLen_Unique, ".bam" )
    uniquesFile <- paste0( baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, "_", captureLen_Unique, "_unique.txt" )
    
    ## Create FASTA files for mapping back to REFGENOMEs
    
    #Extract sequence 3'ends
    start(fragends)[fragends$fe_strand == 3&end(fragends)>(start(fragends)+captureLen_Unique-1)]<-
      end(fragends)[fragends$fe_strand == 3&end(fragends)>(start(fragends)+captureLen_Unique-1)]-captureLen_Unique+1
    
    #Extract sequence 5'ends
    end(fragends)[fragends$fe_strand == 5&start(fragends)<(end(fragends)-captureLen_Unique+1)]<-
      start(fragends)[fragends$fe_strand == 5&start(fragends)<(end(fragends)-captureLen_Unique+1)]+captureLen_Unique-1
    
    do.call( require, args=list( config_genomes[ assemblyName, ] ) )
    seqs <- getSeq( x=base::get( config_genomes[ assemblyName, ] ), names=fragends )
    
    names(seqs) <- fragends$fe_id
    writeXStringSet(seqs,fastaFile)
    
    message(paste('          ### Mapping FASTA file back to',assemblyName))
    
    #cmd <- paste0("bowtie2 -p ",nThreads," -x ",refGenomeFile," -f ",fastaFile," | samtools view -q 1 -hbSu - | samtools sort -T /tmp - > ",bamFile)
    #-q 1: Skip alignments with MAPQ smaller than 1
    #-h Include the header in the output.
    #-b Output in BAM format --> why do we need to output as a BAM file? This is just a temp file.
    #-S Input is in SAM format
    
    #-u Output uncompressed BAM. 
    #This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command. 
    #No sorting is quicker
    
    cmd <-
      paste0(
        "(bowtie2 -p ",
        nThreads_Unique,
        " -x ",
        refGenomeFile,
        " -f ",
        fastaFile,
        " | samtools view -q 1 -hbSu - > ", 
        bamFile, ") 2>&1"
      )
    
    check.trunc <- 0
    while ( check.trunc < 1 ) {
      message( paste( "          ### Bowtie2:", captureLen_Unique ) )
      bowtie.output <- system( cmd, intern=TRUE )
      if ( length( grep( "[main_samview] truncated file", bowtie.output ) ) == 1 ) {
        message( paste( "          ### ERROR: Truncated BAM file: repeating Bowtie mapping" ) )
      } else { 
        check.trunc <- 1
      }
    }
    
    system( paste0( "samtools view ", bamFile, " | grep -v \"XS:\" | cut -f 1 > ", uniquesFile ) )
    ids <- as.numeric( readLines( uniquesFile ) )
    elementMetadata( fragsGR_Unique )[[paste0("len", captureLen_Unique)]] <- fragsGR_Unique$fe_id%in%ids #Use Fragend coordinates. not coordinates for unique mapping. 
    fragsGR_Unique <- fragsGR_Unique[ , order(colnames(elementMetadata(fragsGR_Unique))) ]
    saveRDS( object=fragsGR_Unique, file=fragFile )
    unlink( c( fastaFile, bamFile, uniquesFile ) )
  }
  # if the unique capture len is already in fragGR, we just use that column and rename it to 'unique'
  fragsGR2 <- fragsGR_Unique[, c("pos", "type", "fe_strand", "fe_id", paste0( "len", captureLen_Unique ) ) ]
  colnames( elementMetadata( fragsGR2 ) )[5] <- "unique"
  return( fragsGR2 )
}

getFragMap <- function( vpChr_FragMap=NULL, onlyUnique="TRUE", firstcutter_FragMap="GATC", secondcutter_FragMap="GTAC", genome_FragMap="hg19"
                        , captureLen_FragMap=60, nThreads_FragMap=10, baseFolder_FragMap, Bowtie2Folder, config_genomes
                        ,chr_random=chr_random, chr_fix=chr_fix,chrUn=chrUn, chrM = chrM) {
  
  # here we want to check if we have the genome available for this analysis, in case yes, loading it, otherwise, advice on the 
  # missing genome and switch to the next guy
  # message(paste("      >>> Create frag map for genome :",genome[i],"with re1 :",firstcutter," and re2 :",secondcutter,"and captureLen", captureLen, "<<<"))
  fragFile <- paste0( baseFolder_FragMap, genome_FragMap,"_", firstcutter_FragMap, "_", secondcutter_FragMap, ".rds" )
  
  if ( file.exists( fragFile ) ) {
    # File with unique fragends exists, use it to check for unique fragends for given capture length 
    message( '         ### Loading: ', fragFile )
    fragsGR <- readRDS( fragFile ) 
  } else {
    # Make the digest first
    message(
      paste0(
        "         ### Creating new digest for genome : ",
        genome_FragMap,
        " with re1 : ",
        firstcutter_FragMap,
        " and re2 : ",
        secondcutter_FragMap
      )
    )
    fragsGR <- Digest( assemblyName=genome_FragMap, firstcutter_Digest=firstcutter_FragMap, secondcutter_Digest=secondcutter_FragMap, baseFolder_Digest=baseFolder_FragMap, 
                       config_genomes=config_genomes, chr_random=chr_random, chr_fix=chr_fix,chrUn=chrUn, chrM = chrM) 
  }
  
  message('         ### Retrieve and store the unique fragends')
  fragsGR <- getUniqueFragends(
    fragsGR_Unique=fragsGR
    , firstcutter_Unique=firstcutter_FragMap
    , secondcutter_Unique=secondcutter_FragMap
    , captureLen_Unique=captureLen_FragMap
    , nThreads_Unique=nThreads_FragMap
    , baseFolder_Unique=baseFolder_FragMap
    , Bowtie2Folder=Bowtie2Folder
    , assemblyName=genome_FragMap
    , config_genomes=config_genomes
  )
  if ( !( is.null( vpChr_FragMap ) ) ) {
    if(onlyUnique=="TRUE"){
      message('         ### Selecting unique fragments from ', vpChr_FragMap )
      fragsGR <- subset( fragsGR, as.character( seqnames(fragsGR) ) == vpChr_FragMap & unique == TRUE )  
    }else{
      message('         ### Selecting all fragments from ', vpChr_FragMap )
      fragsGR <- subset( fragsGR, as.character( seqnames(fragsGR) ) == vpChr_FragMap)  
    }
    
  } else {
    if(onlyUnique=="TRUE"){
      message( '         ### Selecting unique fragments from whole genome', vpChr_FragMap )
      fragsGR <- fragsGR[ fragsGR$unique ]
    }else{
      message( '         ### Selecting all fragments from whole genome', vpChr_FragMap )
    }
    
  }
  
  if(length(fragsGR)==0){
    message('No unique fragment ends found for capture length:',captureLen_FragMap)
    
  }
  
  return(fragsGR)
}

plot.chroms <- function(exp,reads, cutoff=0.999, yMax=2500){
  chroms <- as.character(unique(seqnames(reads)))
  total.reads <- sum(reads$reads)
  abs.cut <- quantile(reads$reads,cutoff, na.rm=T)
  reads$reads[reads$reads > abs.cut] <- abs.cut
  reads$reads <- reads$reads/abs.cut
  reads$chr <- sub("chr","",as.character(seqnames(reads)))
  
   for (i in seq_along(chroms)){
     reads[seqnames(reads)==chroms[i]]$chr<-i
   }
   
  #reads[seqnames(reads)=="chrX"]$chr<-length(chroms)-1
  #reads[seqnames(reads)=="chrY"]$chr<-length(chroms)

  reads$chr <-as.numeric(reads$chr)
  
  par(mar=c(5, 4, 4, 2))
  plot(c(0,max(reads$pos)), c(0,length(chroms)), type='n', axes=F, xlab="Chromosome position (Mb)", ylab="", main = exp)
  segments(reads$pos,reads$chr,reads$pos,reads$chr+reads$reads)
  lab <- seq(0,ceiling(max(reads$pos)/10e6)*10, by=10)
  axis(1, at = lab*1e6, label=lab, las=2)
  axis(2, at= 1:length(chroms)+0.5, lab=chroms, lwd=0, las=2, cex.axis=0.9 )
}

make.reads.and.bins <- function( reads, assemblyName, res=25e3, config_genomes ){
  #Extract fragends with reads
  Captures <- reads[reads$reads>0]
  Captures.rds<-Captures[,1]
  Captures.rds$reads<-round(Captures$normReads,2)
  #make bins
  #Get genome info
  do.call( require, args=list( config_genomes[ assemblyName, ] ) )
  assign( 'genome', base::get( config_genomes[ assemblyName, ] ) )
  
  #Chr size
  chr <- as.character(unique(seqnames(reads)))
  x<-seqinfo(genome)
  seqlevels(x) <- chr
  
  message("making bins")
  bin.GR   <- tileGenome(x, tilewidth=res, cut.last.tile.in.chrom=TRUE)
  bin.GR$pos <- start(resize(bin.GR,width=1,fix="center"))
  
  #overlap
  message("Overlapping reads with bins")
  reads<-Captures.rds
  hits <- findOverlaps(reads,bin.GR)
  reads$bin<-0
  reads[queryHits(hits),]$bin<-subjectHits(hits)
  sum.reads<-aggregate(reads$reads, by = list(reads$bin),FUN=sum)
  bin.GR$reads<-0
  bin.GR$reads[sum.reads$Group.1]<-sum.reads$x
  
  return(list(reads=Captures.rds,bins=bin.GR))
}

export.report <- function( RDS.F, OUTPUT.F ){ 
  files <- list.files( path=RDS.F, pattern=".rds", recursive=FALSE )
  report.df <- data.frame()
  for( i in 1:length( files ) ) {    
    vpname <- gsub( "[.].*$", "", files[i] )
    rds <- readRDS( paste0( RDS.F, files[i] ) )
    report <- rds$report
    newrow <- data.frame( vpname=vpname,report )
    report.df <- rbind( report.df, newrow )
  }
  write.table( report.df, file=paste0( OUTPUT.F,"/report.txt"), row.names=FALSE, quote=FALSE )
}

createReport <- function( allReads, mapReads, demuxReads, chromosome, vpPos, normFactor=1e6, wSize, nTop
                          , motifPosperc, readlenperc){
  

  nMapped <- length( mapReads )
  uniqueReads <- sum( allReads$reads )
  uniqueCaptures <- sum( allReads$reads>0 )
  

  #cis reads
  nMappedCis <- length( mapReads[ seqnames( mapReads ) == chromosome ] )
  nMappedCisperc <- round( 100*nMappedCis/nMapped, 2 )
  
  #reads <- allReads[ as.vector( seqnames( allReads ) ) == chromosome ] #If no reads map to cis chr than this crashed the pipeline --> fix
  reads<-allReads[seqnames(allReads) == chromosome]
  
  
  if(length(reads)>0){
  cisReads <- sum( reads$reads )
  percCis<-round( 100*cisReads/uniqueReads, 2 ) #this does not exclude the Top reads..
  cisCaptures <- sum( reads$reads>0 )
  
  
  topIdx <- 1:length( reads ) %in% order( -reads$reads )[ 1:nTop ] #should we do topReads for all reads instead of cis?
  topReads <- sum( reads$reads[ topIdx ] )
  topPct <- round( 100*topReads/cisReads, digits=2 )

  readsGR <- reads[ !topIdx ]
  

  allReads.nTop <- subsetByOverlaps( allReads, reads[ topIdx ], invert=T )
  uniqueReads.nTop <- sum( allReads.nTop$reads )
  cisReads.nTop <- sum( readsGR$reads )
  percCis.nTop <- round( 100*cisReads.nTop/uniqueReads.nTop, 2 )
  
  vpWithin100Kb <- readsGR[ unique( queryHits( findOverlaps( ranges( readsGR ), resize( IRanges( vpPos, vpPos ), width=2e5, fix="center" ) ) ) ) ]
  capt100Kb <- round( 100*mean( vpWithin100Kb$reads>0 ), digits=2 )
  cov100Kb <- round( 100*sum( vpWithin100Kb$reads )/sum( readsGR$reads ), digits=2 )
  
  vpWithin1Mb <- readsGR[ unique( queryHits( findOverlaps( ranges( readsGR ), resize( IRanges( vpPos, vpPos ), width=2e6, fix="center" ) ) ) ) ]
  capt1Mb <- round( 100*mean( vpWithin1Mb$reads>0 ), digits=2 )
  cov1Mb <- round( 100*sum( vpWithin1Mb$reads )/sum( readsGR$reads ), digits=2 )
  

  report <- data.frame(
    nReads=demuxReads # total demux reads
    , motifPosperc=motifPosperc
    , readlenperc=readlenperc
    , nMapped=nMapped # Bowtie mapped read
    , nMappedCis=nMappedCis # Bowtie mapped reads in Cis
    , nMappedCisperc=nMappedCisperc # Bowtie mapped % reads in Cis
    , fragMapped=uniqueReads # Toal unique reads
    , fragMappedCis=cisReads
    , fragMappedCisPerc=percCis
    , fragMappedCisCorr=cisReads.nTop
    , fragMappedCisPercCorr=percCis.nTop
    , nCaptures=uniqueCaptures
    , nCisCaptures=cisCaptures
    , topPct=topPct
    , capt100Kb=capt100Kb
    , cov100Kb=cov100Kb
    , capt1Mb=capt1Mb
    , cov1Mb=cov1Mb
    , stringsAsFactors=FALSE
  )
  }else{
    report <- data.frame(
      nReads=demuxReads # total demux reads
      , motifPosperc=motifPosperc
      , readlenperc=readlenperc
      , nMapped=nMapped # Bowtie mapped read
      , nMappedCis=nMappedCis # Bowtie mapped reads in Cis
      , nMappedCisperc=nMappedCisperc # Bowtie mapped % reads in Cis
      , fragMapped=uniqueReads # Toal unique reads
      , fragMappedCis=0
      , fragMappedCisPerc=NA
      , fragMappedCisCorr=NA
      , fragMappedCisPercCorr=NA
      , nCaptures=uniqueCaptures
      , nCisCaptures=NA
      , topPct=NA
      , capt100Kb=NA
      , cov100Kb=NA
      , capt1Mb=NA
      , cov1Mb=NA
      , stringsAsFactors=FALSE
    )
  }
  
  return( report )
  
}

createPlot <- function( plotTitle, vpPos, chromosome, fragGR, plotLegend=NULL, plotView=1e6, maxY=2500, minY=0, xaxisUnit=c( 'Mb', 'Kb', 'bp' ), plotRegion='cis', foldOut='./', plotType=c( 'PDF', 'PNG' ) ) {
  xaxisUnit <- match.arg( xaxisUnit )
  if( plotRegion == 'cis' ){
    if( xaxisUnit == 'Mb' ){
      scaleValue <- 1e6
    } else if (  xaxisUnit == 'Kb' ) {
      scaleValue <- 1e3
    } else {
      scaleValue <- 1
    }
    zoom <- GRanges( seqnames=chromosome, IRanges( start=vpPos-plotView, end=vpPos+plotView) )
    if ( start( zoom ) < 1 ) {
      start( zoom ) <- 1
    }
    reads.zoom <- fragGR[ unique( findOverlaps( fragGR, zoom )@from ) ]
    if( plotType == 'PDF' ){
      pdf( file=paste0( foldOut, plotTitle, ".pdf" ) )
    } else {
      png( file=paste0( foldOut, plotTitle, ".png" ), width=3200, height=3200, res=300 )
    }
    plot(
      reads.zoom$pos/scaleValue
      , reads.zoom$norm4C
      , type = "h"
      , xlim = c( start(zoom)/scaleValue, end(zoom)/scaleValue )
      , ylim = c( minY, maxY )
      , xlab = paste0( chromosome," - Chromosome position (", xaxisUnit, ")" )
      , ylab = paste0( "4C Coverage / ", scaleValue, " mapped reads" )
      , frame.plot = FALSE
      , main = plotTitle
    )
    if ( !is.null(plotLegend) ){
      legend( x="topleft", bty="n", legend=paste0( names( plotLegend ), " : ", plotLegend ) )
    }
    dev.off()
  } else if( plotRegion == "all" ) {
    if( plotType == 'PDF' ){
      pdf( file=paste0( foldOut, plotTitle, ".pdf" ) )
    } else {
      png( file=paste0( foldOut, plotTitle, ".png" ), width = 1965, height = 1481, res = 300)
    }      
    plot.chroms( exp=plotTitle, reads=fragGR$bins, cutoff=0.999, yMax=maxY )
    dev.off()
  }
}

Run.4Cpipeline <- function( VPinfo.file, FASTQ.F, OUTPUT.F, configuration){
  
  

  
  cutoff = configuration$qualityCutoff
  trim.length = configuration$trimLength
  min.amount.reads=configuration$minAmountReads
  reads.quality=configuration$readsQuality
  map.unique = configuration$mapUnique
  nThreads = configuration$cores
  wSize = configuration$wSize
  nTop = configuration$nTop
  nonBlind = configuration$nonBlind
  make.wig = configuration$wig
  make.BigWig = configuration$bigwig
  make.cisplot = configuration$cisplot
  make.gwplot = configuration$genomePlot
  tsv = configuration$tsv
  bins=configuration$bins
  mmMax=configuration$mmMax
  normFactor=configuration$normFactor
  alnStart=configuration$alnStart
  prefix=configuration$prefix
  maxAmountReads=configuration$maxAmountReads
  
  

  
  #Create output folders
  
  # define folder names. 
  LOG.F <- gsub( x=paste0( OUTPUT.F, "/LOG/" ), pattern='//', replacement='/' )
  FASTQ.demux.F <- gsub( x=paste0( OUTPUT.F, "/FASTQ/" ), pattern='//', replacement='/' )
  TRIM.F <-gsub( x=paste0( OUTPUT.F, "/TRIMMED/" ), pattern='//', replacement='/' )
  BAM.F <-gsub( x=paste0( OUTPUT.F, "/BAM/" ), pattern='//', replacement='/' )
  RDS.F <-gsub( x=paste0( OUTPUT.F, "/RDS/" ), pattern='//', replacement='/' )
  RDS.BIN.F <- gsub( x=paste0( OUTPUT.F, "/RDS-BIN/" ), pattern='//', replacement='/' )
  PLOT.F <- gsub( x=paste0( OUTPUT.F, "/PLOT/" ), pattern='//', replacement='/' )
  WIG.F <- gsub( x=paste0( OUTPUT.F, "/WIG/" ), pattern='//', replacement='/' )
  BWIG.F <- gsub( x=paste0( OUTPUT.F, "/BWIG/" ), pattern='//', replacement='/' )
  GENOMEPLOT.F <- gsub( x=paste0( OUTPUT.F, "/GENOMEPLOT/" ), pattern='//', replacement='/' )
  TSV.F <- gsub( x=paste0( OUTPUT.F, "/TSV/" ), pattern='//', replacement='/' )
  

  logDirs <- list()
  logDirs$outFolder    <- if( !dir.exists( OUTPUT.F )){dir.create( OUTPUT.F )}
  logDirs$logFolder    <- if( !dir.exists( LOG.F )){dir.create( LOG.F )}
  logDirs$fastqFolder  <- if( !dir.exists( FASTQ.demux.F )){dir.create( FASTQ.demux.F )}
  logDirs$trimFolder   <- if( !dir.exists( TRIM.F )){dir.create( TRIM.F )}
  logDirs$bamFolder    <- if( !dir.exists( BAM.F )){dir.create( BAM.F )}
  logDirs$rdsFolder    <- if( !dir.exists( RDS.F )){dir.create( RDS.F )}
  logDirs$rdsbinFolder <- if( !dir.exists( RDS.BIN.F )){dir.create( RDS.BIN.F )}

    
  # Extract experiments from VPinfo, remove weird characters and save in output folder. 
  message( paste0( "\n------ Reading the VP info file: ", VPinfo.file ) )
  VPinfo <- Read.VPinfo( VPinfo.file )
  if (is.null(VPinfo)){
    stop( "viewpoint info file (vpFile) not correct." )
  }
  write.table( VPinfo, paste0( OUTPUT.F, "/VPinfo.txt" ), sep="\t", row.names=FALSE, quote=F )
  
  
  if ( make.cisplot == TRUE ){
    logDirs$plotFolder <- ifelse( !dir.exists( PLOT.F ), dir.create( PLOT.F ), FALSE )
  }
  if ( make.wig == TRUE ){
    logDirs$wigFolder <- ifelse( !dir.exists( WIG.F ), dir.create( WIG.F ), FALSE )
  }
  if ( make.BigWig == TRUE ){
    logDirs$BwigFolder <- ifelse( !dir.exists( BWIG.F ), dir.create( BWIG.F ), FALSE )
  }
  if( any( VPinfo$analysis == 'all' ) & make.gwplot ){
    logDirs$genomeplotFolder <- ifelse( !dir.exists( GENOMEPLOT.F ), dir.create( GENOMEPLOT.F ), FALSE )
  }
  
  if ( tsv == TRUE ){
    logDirs$tsvFolder <- ifelse( !dir.exists( TSV.F ), dir.create( TSV.F ), FALSE )
  }
  
  
  # Make Log files
  LOGNAME <- format( Sys.time(), '%Y_%m_%d__%H_%M_%S' )
  log.path <- paste0( LOG.F, "Log_", LOGNAME, ".txt" )
  demux.log.path <- paste0( LOG.F, "Log_demux_", LOGNAME, ".txt" )
  trim.log.path <- paste0( LOG.F, "Log_4Ctrim_", LOGNAME, ".txt" )
  bowtie.log.path <- paste0( LOG.F, "Bowtie2Log_", LOGNAME, ".txt" )
  
  # Save Run parameters to logfile
  
  run.par <- data.frame(
    param=c( "pipeline.version", "baseFolder", "VPinfo", "FASTQ.F", "OUTPUT.F", "nThreads", "Quality.cutoff", "trim.length", "mmMax", "minReads", "maxReads"
             , "reads.quality", "map.unique", "alnStart","wSize", "nTop", "normFactor", "nonBlind","make.wig","make.bigwig", "make.cisplot", "make.gwplot", 
              "tsv","bins")
    ,value=c( configuration$pipeline.version, configuration$baseFolder, VPinfo.file, FASTQ.F, OUTPUT.F, nThreads, cutoff, trim.length, mmMax, min.amount.reads,maxAmountReads
              ,reads.quality, map.unique, alnStart, wSize, nTop, normFactor
              ,nonBlind, make.wig, make.BigWig,make.cisplot, make.gwplot, tsv,bins)
  )
  
  write.table( run.par, log.path, quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE )                    
  

  
  
  
  # Demultiplex all fastq files and write in FASTQ folder in OUTPUT.F
  message("\n------ Demultiplexing Fastq files based on VPinfo file")
  demux.FASTQ( VPinfo=VPinfo, FASTQ.F=FASTQ.F, FASTQ.demux.F=FASTQ.demux.F, demux.log.path=demux.log.path, mmMax = mmMax)
  
  #4C-seq analysis
  
  exp.name <- as.character( VPinfo$expname )
  primer <-   as.character( VPinfo$primer )
  genome <- as.character( VPinfo$genome )
  firstenzyme <-  as.character( VPinfo$firstenzyme )
  secondenzyme <- as.character( VPinfo$secondenzyme )
  vpChr <- as.character( VPinfo$vpchr )
  vppos <- as.numeric( VPinfo$vppos )
  analysis <- as.character( VPinfo$analysis )
  
  message( "\n------ 4Cseq analysis" )
  
  for ( i in 1:length( exp.name ) ) {
    
    message( paste0( "   +++ experiment: ", exp.name[i], " +++" ) )
    
    
    #skip already generated files
    rdsFile<-paste0( RDS.F, exp.name[i], ".rds" )
    if (file.exists( rdsFile ) ) {
      error.msg <- paste("         ### WARNING: rds file", exp.name[i], "already exist, skipping experiment.")
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    

    primer.sequence <- primer[i]
    
    file.fastq <- paste0( FASTQ.demux.F, exp.name[i], ".fastq.gz" )
    CHR <- paste0(prefix, vpChr[i])

    if ( exp.name[i] %in% exp.name[ duplicated( exp.name ) ] ) {  
      error.msg <- paste0( "      ### ERROR: Experiment name not unique for ", exp.name[i] )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    if ( !( firstenzyme[i] %in% rownames( configuration$enzymes ) ) ) {  
      error.msg <- paste0( "      ### ERROR: firstcutter not found for ", exp.name[i] )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    firstcutter <- configuration$enzyme[ firstenzyme[i], ]
    
    if ( !( secondenzyme[i] %in% rownames( configuration$enzymes ) ) ) {  
      error.msg <- paste0( "      ### ERROR: secondcutter not found for ", exp.name[i] )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    secondcutter <- configuration$enzyme[ secondenzyme[i], ]
    
    if ( !(genome[i] %in% rownames( configuration$genomes )) ) {
      error.msg <- paste0( "      ### ERROR: genome not found for ", exp.name[i] )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    if ( !(analysis[i] %in% c("cis", "all"))) {
      error.msg <- paste0( "      ### ERROR: analysis of cis/all not found for ", exp.name[i] )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    if ( !file.exists( file.fastq ) ) {
      error.msg <- paste0( "      ### ERROR:", exp.name[i], " FASTQ file does not exist" )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    message( "      >>> Trim of the fastq <<<" )
    trim.FASTQ <- trim.FASTQ( exp.name=exp.name[i], primer=primer.sequence, firstcutter=firstcutter, secondcutter=secondcutter, 
                              file.fastq=file.fastq, trim.F=TRIM.F, cutoff=cutoff, trim.length=trim.length, log.path=trim.log.path, 
                              min.amount.reads=min.amount.reads, maxAmountReads=maxAmountReads)
    
    txt.tmp <- paste0( TRIM.F, exp.name[i], ".txt" )
    if ( !file.exists( txt.tmp ) ){
      error.msg <- paste0( "         ### ERROR:", exp.name[i], " trimmed file not found" )
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    captureLen <- trim.FASTQ$captureLen
    
    if (captureLen<=nchar(firstcutter)){
      error.msg <- paste0( "         ### ERROR:", exp.name[i], " capture length <= length firstcutter motif.")
      message("         ### Skipping experiment.")
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    nReads <- trim.FASTQ$nReads
    motifPosperc <-trim.FASTQ$motifPosperc
    readlenperc <- trim.FASTQ$readlenperc
    
    
    
    
    # 3. make BAM files
    
    message( "      >>> Alignment of reads to reference genome <<<" )
    bamFile <- paste0( BAM.F, exp.name[i], ".bam" )
    if ( !file.exists( bamFile ) ) {
      makeBAM( exp.name=exp.name[i], BAM.F=BAM.F, NCORES=nThreads, Bowtie2Folder=configuration$bt2Genomes[ genome[i], ], genome=genome[i], txt.tmp=txt.tmp, log.path=log.path, bowtie.log.path=bowtie.log.path, bamFile=bamFile, map.unique=map.unique, readsQual=reads.quality )
    } else {
      error.msg <- paste("         ### WARNING:", exp.name[i], "BAM file already exist, continuing with exisiting file.")
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
    }
    
    mappedReads <- readGAlignments( file=bamFile, index=bamFile )
    
    
    # To do: lock frag map
    #exec 3>hg19_Dpn2_Csp6I # open a file handle; this part will always succeed
    #flock -x 3      # lock the file handle; this part will block
    #To release the lock:
    #exec 3>&-       # close the file handle
    
    
    
    message( paste0("      >>> Create frag map for genome: ", genome[i], " with RE1: ", firstcutter, " and RE2: ", secondcutter, " and capture length:", captureLen, " <<<" ) )
    frags <- getFragMap(
      vpChr_FragMap=NULL
      ,onlyUnique="TRUE"
      ,firstcutter_FragMap=firstcutter
      ,secondcutter_FragMap=secondcutter
      ,genome_FragMap=genome[i]
      ,captureLen_FragMap=captureLen
      ,nThreads_FragMap=nThreads
      ,baseFolder_FragMap=configuration$baseFolder
      ,Bowtie2Folder=configuration$bt2Genomes[ genome[i], ]
      ,config_genomes=configuration$genomes
      ,chr_random=configuration$chr_random
      ,chr_fix=configuration$chr_fix
      ,chrUn=configuration$chrUn
      ,chrM=configuration$chrM
    )
    
    if(length(frags)==0){
      error.msg <- paste("         ### ERROR:", exp.name[i], "No unique fragment ends found for capture length.")
      message("         ### Skipping experiment.")
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
      
    }
    
    message("      >>> Align reads to fragments <<<")
    readsAln <- alignToFragends(
      gAlign=mappedReads
      ,fragments=frags
      ,firstcut=firstcutter
      ,alnStart=alnStart
    )
    
    if(is.null(readsAln)){
      error.msg <- paste("         ### ERROR:", exp.name[i], "No unique fragment ends aligned.")
      message("         ### Skipping experiment.")
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
      next
    }
    
    message("      >>> Compute statistics and create report <<<")
    

    reportAnalysis <- createReport( allReads=readsAln
                                    , mapReads=mappedReads
                                    , demuxReads=nReads
                                    , chromosome=CHR
                                    , vpPos=vppos[i]
                                    , normFactor=configuration$normFactor
                                    , wSize=wSize
                                    , nTop=nTop
                                    , motifPosperc=motifPosperc
                                    , readlenperc=readlenperc
    )
  
  
    vpInfo <- data.frame(
      name=exp.name[i]
      , Fastq=file.fastq
      , primer=primer.sequence
      , genome=genome[i]
      , firstcutter=firstcutter
      , secondcutter=secondcutter
      , chr=CHR
      , pos=vppos[i]
      , captureLen=captureLen
      , date=Sys.time()
      , Pipeline=configuration$pipeline.version
    )
    
    
    
    
    #Extract Cis reads.
    message("      >>> Compute normalized 4C score per fragment <<<")
    
    rdsFile<-paste0( RDS.F, exp.name[i], ".rds" )
    if (file.exists( rdsFile ) ) {
      error.msg <- paste("         ### WARNING: rds file", exp.name[i], "already exist, overwriting exisiting file.")
      write( error.msg, log.path, append=TRUE )
      message( error.msg )
    }
    
    
    if ( analysis[i] == "cis" | make.cisplot == TRUE ) {
      reads.cis <- norm4C( readsGR=readsAln[ as.vector( seqnames( readsAln ) ) == CHR ], nReads=configuration$normFactor, wSize=wSize, nTop=nTop )
      
      if ( nonBlind ) {
        reads.cis<-reads.cis[reads.cis$type=="non_blind"]
      }
      
    }
    
    if ( analysis[i] == "cis" ) {
      saveRDS( list( reads=reads.cis, report=reportAnalysis, vpInfo=vpInfo ), file=rdsFile)
    }
    
    
    
    
    if ( analysis[i] == "all" ) {
      #normalization all reads is used
      reads.all <- norm4C( readsGR=readsAln, nReads=configuration$normFactor, wSize=wSize, nTop=nTop )
    }
    
    if ( bins == TRUE | make.gwplot == TRUE) {
      if ( analysis[i] == "all" ) {
        message("      >>> Creating bins <<<")
        bin.GR <- make.reads.and.bins( reads=reads.all, assemblyName=genome[i], res=configuration$binSize, config_genomes=configuration$genomes )
        file.out <- paste0( RDS.BIN.F, exp.name[i], "_bin_res_", (configuration$binSize/1e3), "kb.rds" )
        saveRDS( object=list( reads=bin.GR$reads, bins=bin.GR$bins, report=reportAnalysis, vpInfo=vpInfo ), file=file.out, compress=TRUE )
      }
    }
    
    if ( analysis[i] == "all" ) {
      if ( nonBlind ) {
        reads.all<-reads.all[reads.all$type=="non_blind"]
      }
      saveRDS( list( reads=reads.all, report=reportAnalysis, vpInfo=vpInfo ), file=rdsFile)      
    }
    
    if ( make.cisplot == TRUE ) {
      message("      >>> Creating local 4C Plot <<<")
      createPlot( 
        plotTitle=exp.name[i]
        , vpPos=vppos[i]
        , chromosome=CHR
        , fragGR=reads.cis
        , plotLegend=reportAnalysis
        , plotView=configuration$plotView
        , maxY=configuration$maxY
        , minY=0
        , xaxisUnit=configuration$xaxisUnit
        , plotRegion='cis'
        , foldOut=PLOT.F
        , plotType=configuration$plotType
      )
    }
    
    if ( analysis[i]=="all" & make.gwplot == TRUE ) {
      message("      >>> Creating genome-wide 4C coverage Plot <<<")
      createPlot( 
        plotTitle=exp.name[i]
        , vpPos=vppos[i]
        , chromosome=CHR
        , fragGR=bin.GR
        , plotLegend=reportAnalysis
        , plotView=configuration$plotView
        , maxY=configuration$maxY
        , minY=0
        , xaxisUnit=configuration$xaxisUnit
        , plotRegion='all'
        , foldOut=GENOMEPLOT.F
        , plotType=configuration$plotType
      )
      
    }
    
    #7. make WIG
    if ( make.wig == TRUE ) {
      message("      >>> Create WIG file <<<")
      
      if ( nonBlind ) {
        wigFile <- paste0( WIG.F, exp.name[i], "_nonblind_WIN",wSize,".wig" )  
      }else{
        wigFile <- paste0( WIG.F, exp.name[i],"_WIN",wSize ,".wig" )  
      }
      
      if ( file.exists( wigFile ) ){
        error.msg <- paste( "         ### WARNING:", exp.name[i], "WIG file already exist." )
        message( error.msg )
        write( error.msg, log.path, append=TRUE )
      } else {
        if ( analysis[i] == "cis" ) {
          reads.wig <- reads.cis[ order( reads.cis$pos ) ]
        }else{
          reads.wig <- reads.all[ order( seqnames(reads.all),reads.all$pos ) ]
        }
        
        exportWig( gR=reads.wig, expName=exp.name[i], filename=wigFile, vpPos=vppos[i], vpChr=CHR, plotView=configuration$plotView)
      }
    }
    
    #make bigwig
    
    
      
    if ( make.BigWig == TRUE ) {
      message("      >>> Create Big WIG file <<<")
      
      if ( nonBlind ) {
        BwigFile <- paste0( BWIG.F, exp.name[i], "_nonblind_WIN",wSize,".bw" )  
      }else{
        BwigFile <- paste0( BWIG.F, exp.name[i],"_WIN",wSize ,".bw" )  
      }
      
      if ( file.exists( BwigFile ) ){
        error.msg <- paste( "         ### WARNING:", exp.name[i], "Big WIG file already exist." )
        message( error.msg )
        write( error.msg, log.path, append=TRUE )
      } else {
        if ( analysis[i] == "cis" ) {
          reads.bwig <- reads.cis[ order( reads.cis$pos ) ]
        }else{
          reads.bwig <- reads.all[ order( seqnames(reads.all),reads.all$pos ) ]
        }
        
       
        exportBigWig(GR=reads.bwig, OutFile=BwigFile, assemblyName=genome[i],config_genomes=configuration$genomes)
      }
    }
    
    
    
    
    
    #Make tsv
    if ( tsv == TRUE ) {
      message("      >>> Create TSV file <<<")
      
      if ( nonBlind ) {
        tsvFile <- paste0( TSV.F, exp.name[i], "_nonblind_WIN",wSize,".tsv" )  
      }else{
        tsvFile <- paste0( TSV.F, exp.name[i],"_WIN",wSize ,".tsv" )  
      }
      
      
      
      if ( file.exists( tsvFile ) ){
        error.msg <- paste( "         ### WARNING:", exp.name[i], "TSV file already exist." )
        message( error.msg )
        write( error.msg, log.path, append=TRUE )
      } else {
        if ( analysis[i] == "cis" ) {
          reads.tsv <- reads.cis[ order( reads.cis$pos ) ]
        }else{
          reads.tsv <- reads.all[ order( seqnames(reads.all),reads.all$pos ) ]
        }
        
        
        gzname<- paste0( tsvFile, ".gz" )
        gz1 <- gzfile( gzname, "w" )
        tsvdat <- data.frame(chr=seqnames(reads.tsv)
                             ,start=start(reads.tsv)
                             ,end=end(reads.tsv)
                             ,pos=reads.tsv$pos
                             ,type=reads.tsv$type
                             ,fe_strand=reads.tsv$fe_strand
                             ,fe_id=reads.tsv$fe_id
                             ,reads=reads.tsv$reads
                             ,normReads=reads.tsv$normReads
                             ,norm4C=reads.tsv$norm4C)
        
        
        write.table( tsvdat, file=gz1, append=FALSE, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE )
      }
      close(gz1)
    }
    
    
    message('\n')
  }
  #make report from all RDS files
  message( "------ Writing final report\n" )
  export.report( RDS.F, OUTPUT.F )
}

getVPReads <- function(rds,vpRegion=2e6) {
  
  reads <- rds$reads
  vppos <- rds$vpInfo$pos
  vpChr <- rds$vpInfo$chr
  zoom <- GRanges(seqnames=vpChr, resize(IRanges(vppos,vppos),width=vpRegion,fix="center") )
  vpGR <- reads[unique(queryHits(findOverlaps(reads,zoom)))]
  
  peakCDat <- data.frame(pos=vpGR$pos,reads=vpGR$reads)
  
  return(peakCDat)
  
}

getPeakCPeaks <- function(resPeakC,min.gapwidth=4e3) {
  
  if(length(resPeakC$peak)>0){
    
    vpChr <- resPeakC$vpChr
  
    peakRanges <- reduce(IRanges(resPeakC$peak,resPeakC$peak),min.gapwidth=min.gapwidth)
    peakGR <- GRanges(seqnames=vpChr,ranges=peakRanges)
  
  } else {
  
    peakGR <- NULL
    
  }
  
  return(peakGR)
  
}

exportPeakCPeaks <- function(resPeakC,bedFile,name=NULL,desc=NULL,includeVP=TRUE,min.gapwidth=4e3) {
  
  vpChr <- resPeakC$vpChr
  vpPos <- resPeakC$vpPos
  
  if(is.null(name)) {
    
    name <- "peakC_track"
    
  }
  
  if(is.null(desc)) {
    
    desc <- paste0("peakC peaks on ",vpChr)
    
  }
  
  
  
  vpGR <- resize(GRanges(seqnames=vpChr,IRanges(vpPos,vpPos)),width=1e3,fix="center")
  
  browserPos <- paste0(vpChr,":",vpPos-1e6,"-",vpPos+1e6)
  browserPosLine <- paste("browser position",browserPos,"\n")
  cat(browserPosLine,file=bedFile)
  
  trackLine <- paste0('track name=\"',name,'\" description=\"',desc,'\" visibility=2 itemRgb=\"On\"',"\n")
  cat(trackLine,file=bedFile,append=TRUE)
  
  if(min.gapwidth>0) {
    
    resPeakC$exportPeakGR <- getPeakCPeaks(resPeakC,chr=resPeakC$vpChr,min.gapwidth=min.gapwidth)
    
  }
  
  if(includeVP) {
    
    if(!is.null(resPeakC$exportPeakGR)){
    
      bedDF <- as.data.frame(resPeakC$exportPeakGR)[,1:3]
      colnames(bedDF)[1] <- "chr"
      write.table(bedDF,file=bedFile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    
    }
    
  } else {
    
    if(!is.null(resPeakC$exportPeakGR)){
      
      bedDF <- as.data.frame(resPeakC$exportPeakGR)[,1:3]
      colnames(bedDF)[1] <- "chr"
      write.table(bedDF,file=bedFile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    
    }
  
  }
  
  
}

doPeakC <- function(rdsFiles, vpRegion=2e6, wSize=21,alphaFDR=0.05,qWd=1.5,qWr=1,minDist=15e3,min.gapwidth=4e3) {
  
  if( !suppressMessages(require( "peakC", character.only=TRUE ) ) ) stop( "Package not found: peakC" )
  if( !suppressMessages(require( "GenomicRanges", character.only=TRUE ) ) ) stop( "Package not found: GenomicRanges" )
  
  # Analyze num.exp replicated 4C-Seq experiments with peakC 
  
  num.exp <- length(rdsFiles)
  
  message("Performing peakC on ",num.exp, " experiments.")
  
  # IF num.exp>1, FIRST CHECK IF vppos is equal for all exps !! 
  
  if(num.exp>1) {
    
    peakCDat <- list()
    vppos <- vector()
    
    for(i in 1:num.exp) {
      
      if(file.exists(rdsFiles[i])) {
        message("Loading data for experiment: ", rdsFiles[i])
        rds <- readRDS(rdsFiles[i])
        vppos[i] <- rds$vpInfo$pos
        vpChr <- as.vector(rds$vpInfo$chr)
        peakCDat[[i]] <- getVPReads(rds=rds,vpRegion=vpRegion)
      }else {
        stop( paste0("File not found: ",rdsFiles[i]) )
      }
      
    }
    
    
    vppos <- unique(vppos)
    vpChr <- unique(vpChr)
    
    message("viewpoint position: ",vppos)
    
    
    if(length(vppos)==1){
      resPeakC <- suppressWarnings(combined.analysis(data=peakCDat,num.exp=num.exp,vp.pos=vppos,wSize=wSize,alphaFDR=alphaFDR,qWr=qWr,minDist=minDist))
    }else{
      stop(paste0("Viewpoint positions not unique"))
    }
    
    
  }
  
  if(num.exp==1){
    
    rds <- readRDS(rdsFiles[1])
    vppos <- rds$vpInfo$pos
    vpChr <- as.vector(rds$vpInfo$chr)
    peakCDat <- getVPReads(rds=rds,vpRegion=vpRegion)
    resPeakC <- suppressWarnings(single.analysis(data=peakCDat,vp.pos=vppos,wSize=wSize,qWd=qWd,qWr=qWr,minDist=minDist))
    
    if(length(resPeakC$peak)==0){
    
        message("No significant peaks are found, returning empty peak list.")
      
    }
  }
  
  resPeakC$vpPos <- vppos
  resPeakC$vpChr <- vpChr
  
  if(length(resPeakC$peak)>0){
  
  	resPeakC$exportPeakGR <- getPeakCPeaks(resPeakC,min.gapwidth=min.gapwidth)
  
  } else{
    
    message("No significant peaks are found, returning empty peak list.")   
  	resPeakC$exportPeakGR <- NULL
  	
  }
  
  return(resPeakC)
  
}

embryo <- function(){
  
  
  logo<-c(" "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","i",",",".","r"," ","r","r","s","i",";",":",";","s","i","r"," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",";",":",";",";","X","X","5","X","X","2","2","A","2","2","s",";","2",",","i",";","."," "," "," "," "," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","i","s",";","A","A","3","X","A","3","h","h","5","A","3","M","2","2","2","A","h","2","A","r","5","i",":"," "," "," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","i",".","r","X","5","r","M","2","M","G","s","H","H","S","h","M","2","G","H","2","M","r","S","h","3","5","i","r","s"," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",";","r",".","2","i","5","G","h","G","#","3","@","B","2","G","h","G","&","2","h","A","G","G","h","3","M","A","2","2",".","X"," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",",",":","X","2","X","h","X","9","H","5","#","G","#","9","H","H","M","G","S","M","h","M","M","M","2","H","2","5","X","3","A","r"," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","X","s","i","h","2","S","5","3","M","#","#","5","3","#","#","9","M","5","9","G","#","G","A","&","B","3","3","5","H","2","h","i","i",";"," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","r",":","2","h","2","h","#","#","#","S","5","#","G","5","B","h","M","&","A","B","#","9"," ","2","3","&","h","A","G","s","H","3",",","s",","," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","s",".","2","r","5","M","s","S","S","2","&","h","S","B","h","S","#","3","&","h","H","B","M","&","S","2","#","M","H","H","A","3","5","i",";","."," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",".","s",";","r","s","5","G","2","3","&","H","#","5","M","B","G","3","h","h","9","h","2","G","M","9","G","3","S","G","A","h","s","A","A",",",","," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",":"," ","r","A","r","5","2","H","M","M","3","G","B","2","h","h","G","#","5","5","S","S","H","M","h","9","G","M","M","3","9","3","s","5","i","i"," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","r",".",";","2","i","3","X","2","3","3","3","M","M","H","M","A","B","3","2","H","r","M","M","s","3","2","3","H","A","9","h","A","2",";",","," ",
          " "," "," "," "," "," "," ",";",".",":","r"," "," "," "," "," "," "," "," ",";",":","i",";","2","X","2","3","5","5","M","A","3","A","5","M","A","M","3","2","2","A","s","i","5","5","h","A","H","h","3","5","i",";","."," ",
          " "," "," "," "," "," "," ",";","s",";",":",";","."," "," "," "," "," "," "," "," ",":",":",",","A","s","i","X","s","A","s","s","A","A","A","X","r","r","X",";",".",",","i",";","r","X","2","M","h","r","3",",","s"," "," ",
          " "," "," "," "," "," ",".",";","s","s","r",".",";"," "," "," "," "," "," "," "," "," ",","," ",",",":",",","A"," ",";","i",";","A","i",",",";",";",":",","," "," ","."," ",";","i","5","h","A","G","5",":","X",":"," "," ",
          " "," "," "," "," "," "," ","s","i","2","h"," ","i",","," "," "," "," "," "," "," "," "," "," ",":"," ",",",":",":","r"," ","s",",",":","i",".",";"," ",":"," ",".",".",".",",","A","h","2","M","X","3","5",";","i"," "," ",
          " "," "," "," "," "," "," ",";","X",";","M","3",":",":",";",","," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",",",".",":","2","X","h","5","S","h","s","5",";","r"," "," ",
          " "," "," "," "," "," "," ",".","5","r","s","X","s","r",".",";"," ","."," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",":",";","X","2","s","H","h","2","A","2",":",":"," "," "," ",
          " "," "," "," "," "," "," "," "," ","r","s",";","5","r",".","i"," ",";"," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",","," ",";","r","3","A","h","h","2","h","A",":","s"," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," ","A","i","i","2","X","r",";",":"," ",".","."," "," "," "," "," "," "," "," "," "," "," "," ",".","."," ",":",".","3",":","2","X","5","H","M","2","A","s","s","r",","," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," ",":",";","3","5","X","2","i","r","r",".","s",":",",",".",","," ",".","."," ",","," ",","," ",":","r",",","5",";","r","5","5","M","M","i","A","5","s","A","."," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," ",".","X","X","3","2","s","s",":","r",".","i",".",".",":"," ",":"," "," ","i",".","i",":",".","X","r","3","M","A","3","A","5","A","A","i",",","s"," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," ","i","r",":","X","2","A","A","s",":",";","s","r",";",",",";","r","A",",","X","X","5","5","A","5","A","3","X","X","X",":","i","r",",",":","."," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","i",".","3","A",",","h","h","2","r","r","2","5","X","r","A","2","M","A","X","5","2","H","2","s",";","5",",","2","."," "," "," "," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",":",".","X","i",".","5",";","h","5",",","3","A","r","r","i",";","A","i","A",";"," ","5"," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",";"," ",";",".",",",";",":",":","r"," ",";","i",".",";"," ",".",":"," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",
          " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",","," ",".",","," ",","," ",",",","," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "
  )
  
  logo.mat<-matrix(logo, nrow=28, ncol=55, byrow = T)
  
  for (i in 1:28) {
    for (j in 1:55) {
      cat(logo.mat[i,j])
    }
    cat('\n')
  }
  
}


exportBigWig <- function(GR, OutFile, assemblyName,config_genomes){
  
  #load the package
  require(rtracklayer)
 
  do.call( require, args=list( config_genomes[ assemblyName, ] ) )
  assign( 'genome', base::get( config_genomes[ assemblyName, ] ) )
  
  
  message("Converting RDS to BigWig file: ",OutFile)
    
  #Maybe better to keep original frag end coordinates?
  
    start(GR)<-GR$pos
    end(GR)  <-GR$pos
    GR$score <-GR$norm4C
    
    #seqinfo(GR)<-keepStandardChromosomes(seqinfo(genome))  
    seqlevels(GR)<-seqlevels(genome)
    seqinfo(GR)<-seqinfo(genome)
    
    rtracklayer::export(GR, OutFile)
    
    
    
    
    

    
    
    
    
 
}