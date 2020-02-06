# Runs ExomeDepth over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runExomeDepth.r [ExomeDepth_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ExomeDepth))
library(methods)
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

# returns gene for position matching ROI in bed data
auxGetGene <- function(bedData, chr, pos){
  gene <- bedData[bedData$chr == chr & (pos >= bedData$start & pos <= bedData$end),"gene"]
  if (length(gene) == 0)
    return(NULL)
  else
    return(gene)
}
  

# Process ExomeDepth alg. body
processExomedepthBody <- function(testCountsDF, controlCountsDF, countsDef, params){
  all <- data.frame()
  
  for (i in 1:ncol(testCountsDF)) {
    
    # build necessary matrix
    sampleName = colnames(testCountsDF)[i]
    if (sampleName %in% colnames(controlCountsDF)){
      controlCounts <- as.matrix(controlCountsDF)[, -i]
    } else 
      controlCounts <- as.matrix(controlCountsDF)
    testCounts <- as.matrix(testCountsDF)[,i]
    
    # define final controls set
    references = select.reference.set(test.counts = testCounts,
                                      reference.count = controlCounts,
                                      bin.length = (countsDF$end - countsDF$start) / 1000, 
                                      n.bins.reduced = 10000,
                                      phi.bins = params$phi.bins)
    controls = apply(X = as.matrix(controlCounts[, references$reference.choice]), MAR=1, FUN=sum)

    
    # call cnvs
    all_exons = new('ExomeDepth',
                    test = testCounts,
                    reference = controls,
                    formula = 'cbind(test, reference) ~ 1')
    all_exons = CallCNVs(x = all_exons,
                         transition.probability = params$transition.probability,
                         expected.CNV.length = params$expected.CNV.length,
                         chromosome = countsDef$space,
                         start = countsDef$start,
                         end = countsDef$end,
                         name = countsDef$names)
    
    if (nrow(all_exons@CNV.calls) > 0){
      # add sample column: remove first char "X"
      all_exons@CNV.calls$sample <- sampleName
      
      # add gene column
      for(i in 1:nrow(all_exons@CNV.calls)) {
        row <- all_exons@CNV.calls[i,]
        row$gene <- auxGetGene(bedData, row$chromosome, row$end)
        if (is.null(row$gene))
          stop(paste("Error: gene not found for chromosome", row$chromosome, ", start", row$start, ", end", row$end))
        
        all_exons@CNV.calls[i, "Gene"] <- row$gene
      }
    }
    
    
    all <- rbind(all, all_exons@CNV.calls)
  }
  
  return(all)
}


# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  exomedepthParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  exomedepthParamsFile <- "algorithms/exomedepth/exomedepthParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
}

#Load the parameters file
params <- yaml.load_file(exomedepthParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# print params
print(paste("Params for this execution:", list(params)))


# go over datasets and run ExomeDepth for those which are active
for (name in names(datasets)) {

  dataset <- datasets[[name]]
  
  if (dataset$include){
    print(paste("Starting ExomeDepth for", name, "dataset", sep=" "))
    
    
    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    bedData <- read.table(bedFile, sep="\t", stringsAsFactors=FALSE, col.names = (c("chr", "start", "end", "gene")))
    fastaFile <- file.path(dataset$fasta_file)
    
    # set readlength from algorithm params if defined
    if (!is.null(params$readLength)){
      readLength <- params$readLength
    } else
      readLength <- dataset$read_length
    
    # build output folder and file
    if (!is.null(params$outputFolder)) {
      outputFolder <- params$outputFolder
    } else
      outputFolder <- file.path(getwd(), "output", paste0("exomedepth-", name))  
    if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
      unlink(outputFolder, recursive = TRUE);
      dir.create(outputFolder)
    }
    outputFile <- file.path(outputFolder, "all_cnv_calls.txt")
    failedRegionsFile <- file.path(outputFolder, "failedROIs.csv")
    
    # Do pre-calc part of the algorithm
    if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
      
      # read bam counts   
      bamFiles <- list.files(bamsDir, "*.bam$", full.names=T)
      counts <- getBamCounts(bed.file = bedFile,
                             bam.files = bamFiles,
                             read.width = readLength,
                             referenceFasta = fastaFile)
      countsDF  <- as.data.frame(counts)
      names(countsDF) <- gsub(".bam", "", names(countsDF)) # remove .bam from sample name
      
      # Fix rare bug: sometimes X is added at the beginning of input samples
      parts <- strsplit(bamFiles[1], "/")
      original <- gsub(".bam", "", parts[[1]][length(parts[[1]])])
      if (substr(original, 1, 1) != "X" && substr(names(countsDF)[7], 1, 1) == "X")
        names(countsDF)[7:ncol(countsDF)] <- substring(names(countsDF)[7:ncol(countsDF)] , 2)
      
      if (!is.null(params$execution) && params$execution == "onlyPrecalcPhase") {
        # save results to be used by other executions
        dir.create(params$precalcFolder)
        saveRDS(countsDF, file.path(params$precalcFolder, "countsDF.rds"))
        
        print(paste("ExomeDepth (Only pre-calc phase) for", name, "dataset finished", sep=" "))
        cat("\n\n\n")
        quit()
      }
    } else {  # skipPrecalcPhase mode: read previous results
      print(paste("ExomeDepth Skipping pre-calc phase for", name, "dataset finished", sep=" "))
      countsDF <- readRDS(file.path(params$precalcFolder, "countsDF.rds"))
    }

    
    # Process each sample
    all <- data.frame()
    if (!is.null(dataset$clinical_indication) && dataset$clinical_indication != ""){
      sample_indications <- read.table(dataset$clinical_indication, header = T, sep = "\t", stringsAsFactors=F)[,c(1,2)]
      sample_indications <- sample_indications[sample_indications$Genes != "" & sample_indications$Genes != " ", ] # remove empty items
      indicationsComp <- readCompatibleIndications(dataset$clinical_indication) # get compatible indications for each indication

      for (indication in names(indicationsComp)) {
        print(paste("\nProcessing indication:", indication))


        # samples to test: samples matching this indication
        testSamplesNames <- sample_indications[sample_indications$Genes == indication,][["SampleID"]]
        if (length(testSamplesNames > 0)){

          # Calculate available control samples: those compatible with test indication
          compatibles <- indicationsComp[[indication]]
          samplesToExclude <- sample_indications[!sample_indications$Genes %in% compatibles,][["SampleID"]]
          controlSamplesNames <- sample_indications[!sample_indications$SampleID %in% samplesToExclude, ][["SampleID"]]

          # Call other algorithm steps
          results <- processExomedepthBody(testCountsDF = countsDF[testSamplesNames],
                                           controlCountsDF = countsDF[controlSamplesNames],
                                           countsDef = countsDF[1:6],
                                           params = params)
          all <- rbind(all, results)
        }
      }

    } else {
      message("Clinical indication not found, using all samples potentially as controls")

      all <- processExomedepthBody(testCountsDF = countsDF[7:ncol(countsDF)],
                                   controlCountsDF = countsDF[7:ncol(countsDF)],
                                   countsDef = countsDF[1:6],
                                   params = params)
    }

    # write results
    write.table(all, file = outputFile, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    df <- read.table(outputFile,sep="\t", stringsAsFactors=FALSE, header = TRUE) # REMOVE
    
    # Save results in GRanges format
    message("Saving GenomicRanges results")
    saveResultsFileToGR(outputFolder, "all_cnv_calls.txt", chrColumn = "chromosome", sampleColumn = "sample",
                        startColumn = "start", endColumn = "end", cnvTypeColumn = "type")
    
    # Save empty failed regions (not available in this algorithm)
    fr <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(fr) <- c("SampleID", "Chr", "Start", "End", "Gene")
    write.table(fr, failedRegionsFile, sep="\t", row.names=FALSE, quote = FALSE)
    
    print(paste("ExomeDepth for", name, "dataset finished", sep=" "))
    cat("\n\n\n")   
  }
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)
