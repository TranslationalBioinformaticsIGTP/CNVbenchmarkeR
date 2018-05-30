# Runs CoNVaDING over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runConvading.R [convading_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(stringr))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

# translates DEL/DUP into a common format
auxConvert <- function(x) {
  if (x == "DEL") return("deletion") 
  else if (x == "DUP") return("duplication")
  else return("")
}

# Saves csv file with all failed exons (in a common format)
saveFailedROIs <- function(outputFolder){
  outputFile <- file.path(outputFolder, "failedROIs.csv")
  output <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(output) <- c("SampleID", "Chr", "Start", "End", "Gene")
  
  for(file in dir(outputFolder, pattern="best.score.totallist.txt")) {
    # read data and filter for ROIs
    sampleData <- read.table(file.path(outputFolder, file), sep="\t", header=TRUE, stringsAsFactors = FALSE)
    sampleData <- sampleData[sampleData$QUALITY %in% c("LOW_QUALITY", "LOW_QUALITY,FAILED"),] 
    
    # get sample name
    sampleName <- strsplit(file, "\\.best")[[1]][1]
    
    # add to final output
    output <- rbind(output, data.frame(SampleID = sampleName, Chr = sampleData[, 1], Start = sampleData[, 2], End = sampleData[, 3], Gene = sampleData[, 4]))
  }
  
  # save output file
  write.table(output, outputFile, sep="\t", row.names=FALSE, quote = FALSE)
}


# Process convading alg. body
processConvadingBody <- function(launchFile, inputFolder, controlsFolder, targetQcList, params){

  # from .txt tp select best control samples
  system(paste("perl", launchFile, "-mode StartWithMatchScore", "-inputDir" , inputFolder, "-controlsDir" , controlsFolder, "-outputDir", inputFolder))

  # cnv calling
  system(paste("perl", launchFile, "-mode StartWithBestScore", "-inputDir" , inputFolder,
               "-regionThreshold" , params$regionThreshold, "-ratioCutOffLow" , params$ratioCutOffLow, "-ratioCutOffHigh" , params$ratioCutOffHigh,
               "-zScoreCutOffLow" , params$zScoreCutOffLow, "-zScoreCutOffHigh" , params$zScoreCutOffHigh, "-outputDir", inputFolder))

  # generate the list of targets and corresponding quality thresholds
  system(paste("perl", launchFile, "-mode GenerateTargetQcList", "-inputDir" , inputFolder, "-controlsDir" , controlsFolder,
               "-regionThreshold", params$regionThreshold, "-ratioCutOffLow" , params$ratioCutOffLow, "-ratioCutOffHigh",
               params$ratioCutOffHigh, "-zScoreCutOffLow" , params$zScoreCutOffLow, "-zScoreCutOffHigh" , params$zScoreCutOffHigh,
               "-sampleRatioScore", params$sampleRatioScore, "-outputDir", inputFolder))

  
  # create final list
  system(paste("perl", launchFile, "-mode CreateFinalList", "-inputDir" , inputFolder, "-targetQcList", targetQcList,
               "-percentageLessReliableTargets", params$percentageLessReliableTargets, "-outputDir", inputFolder))
  
}


# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
    convadingParamsFile <- args[1]
    datasetsParamsFile <- args[2]
} else {
    convadingParamsFile <- "algorithms/convading/convadingParams.yaml"
    datasetsParamsFile <- "datasets.yaml"
}

#Load the parameters file  
params <- yaml.load_file(convadingParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# extract convading params
convadingFolder <- file.path(params$convadingFolder)
print(paste("Params for this execution:", list(params)))


# go over datasets and run convading for those which are active
for (name in names(datasets)) {
    dataset <- datasets[[name]]
    if (dataset$include){
        print(paste("Starting convading for", name, "dataset", sep=" "))
      
        # extract fields
        bamsDir <- file.path(dataset$bams_dir)
        bedFile <- file.path(dataset$bed_file)
        fastaFile <- file.path(dataset$fasta_file)

        # Create output folder
        if (!is.null(params$outputFolder)) {
          outputFolder <- params$outputFolder
        } else
          outputFolder <- file.path(getwd(), "output", paste0("convading-", name))  
        if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
          unlink(outputFolder, recursive = TRUE);
          dir.create(outputFolder)
        }

        # build input/output file paths
        launchFile <- file.path(convadingFolder, "CoNVaDING.pl")
        controlsFolder <- file.path(outputFolder, "controls") # Where controls will be created
                
        if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
          
          # Call convading...
          # from .BAM to normalized coverage .txt
          system(paste("perl", launchFile, "-mode StartWithBam", "-inputDir" , bamsDir, "-useSampleAsControl -controlsDir" , controlsFolder, "-bed", bedFile, "-outputDir", outputFolder))
          
          
          if (!is.null(params$execution) && params$execution == "onlyPrecalcPhase") {
            print(paste("CoNVaDING (Only pre-calc phase) for", name, "dataset finished", sep=" "))
            cat("\n\n\n")
            quit()
          }
          
        } else { # skipPrecalcPhase mode: copy previous results from precalc folder to output folder
          # copy files to ouput folder
          flist <- list.files(params$precalcFolder,full.names = TRUE)
          dir.create(outputFolder, recursive = TRUE)
          file.copy(flist, outputFolder)
          flist <- list.files(file.path(params$precalcFolder, "controls"),full.names = TRUE)
          dir.create(controlsFolder)
          file.copy(flist, controlsFolder)
        }
      


        if (dataset$validated_results_file_format == "panelcn" && dataset$validated_results_file != "" && params$defineControlsByIndication == TRUE){
          sample_indications <- readIndicationsForSamples(dataset$validated_results_file)
          indicationsComp <- readCompatibleIndications(dataset$validated_results_file) # get compatible indications for each indication
          allFiles <- list.files(outputFolder, "*coverage.txt$", full.names=T)
          
          # obtain samples withot any indication
          samplesWithoutIndication <- c()
          for(file in allFiles){
            
            # extract sample name from file
            sName <- tail(strsplit(strsplit(file, "\\.")[[1]][1], "/")[[1]], 1)
            
            # add test sample to copy later
            if (!sName %in% sample_indications$SampleID)
              samplesWithoutIndication <- c(samplesWithoutIndication, sName)
          }
          
          # process each indication
          for (indication in names(indicationsComp)) {
            print(paste("Processing indication:", indication))

            # Create folders for each indication
            indicationFolder <- file.path(outputFolder, paste0("indication", str_replace_all(indication, pattern=" ", repl="")))
            dir.create(indicationFolder)
            indicationControlsFolder <- file.path(indicationFolder, "controls")
            dir.create(indicationControlsFolder)

            targetQcList <- file.path(indicationFolder, "targetQcList.txt")

            # samples to test: samples matching this indication
            testSamplesNames <- sample_indications[sample_indications$Genes == indication,][["SampleID"]]
            if (length(testSamplesNames > 0)){

              # Copy test sample files
              filesToCopy <- c()
              for(file in allFiles){

                # extract sample name from file
                sName <- tail(strsplit(strsplit(file, "\\.")[[1]][1], "/")[[1]], 1)

                # add test sample to copy later
                if (sName %in% testSamplesNames)
                  filesToCopy <- c(filesToCopy, file)
              }
              file.copy(filesToCopy, indicationFolder)

              # copy control files (samples with compatible indications)
              compatibles <- indicationsComp[[indication]]
              samplesToExclude <- sample_indications[!sample_indications$Genes %in% compatibles,][["SampleID"]]
              controlSamplesNames <- sample_indications[!sample_indications$SampleID %in% samplesToExclude, ][["SampleID"]]
              filesToCopy <- c()
              for(file in allFiles){

                # extract sample name from file
                sName <- tail(strsplit(strsplit(file, "\\.")[[1]][1], "/")[[1]], 1)

                # add test sample to copy later
                if (sName %in% controlSamplesNames | sName %in% samplesWithoutIndication)
                  filesToCopy <- c(filesToCopy, file)
              }
              file.copy(filesToCopy, indicationControlsFolder)
            }

            processConvadingBody(launchFile, indicationFolder, indicationControlsFolder, targetQcList, params)
            file.copy(list.files(indicationFolder, "*list.txt$", full.names=T), outputFolder)
          }

        } else {
          targetQcList <- file.path(outputFolder, "targetQcList.txt")
          processConvadingBody(launchFile, outputFolder, controlsFolder, targetQcList, params)
        }

        # Add cnv.Type column to match common format
        for (file in list.files(outputFolder, "*list.txt$", full.names=TRUE)) {
          print(paste("Adding CNV.type column to file ", file))
          data <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=T)
          data$CNV.type <- lapply(data$ABBERATION, function(x) sapply(x, auxConvert)) # Add CNV.type column before storing file
          data$CNV.type <- as.factor(unlist(data$CNV.type))
          write.table(data, file, sep="\t", row.names=FALSE, quote = FALSE)
        }
        
        # Save results in GRanges format
        message("Saving GenomicRanges results")
        saveResultsFolderToGR(outputFolder, "best.score.longlist.txt", chrColumn = "CHR", startColumn = "START", endColumn = "STOP")
        
        # Save failed ROIs in common format
        saveFailedROIs(outputFolder)
        
        print(paste("convading for", name, "dataset finished", sep=" "))
        cat("\n\n\n")
    }
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)