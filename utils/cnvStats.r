# Stats classes definitions, used to generate specificity, sensibility and other stats results
# Used by summary.r
suppressPackageStartupMessages(library(GenomicRanges))
source(if (basename(getwd()) == "utils") "utils.r" else if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

#####  SummaryStats class definition for a dataset #####
SummaryStats <- setRefClass("SummaryStats",

                            fields = list(datasets = "list",  # list of datasets storing a list of samples each one. Key: dataset name; Value: list fo samples
                                          algorithms = "character",  # array of loaded algorithms (list of characters)
                                          globalStats = "list" ),  # list of list of global stats.
                            # First Key: dataset name; Second key: alg name. Value: datata frame (cols: TP,TN,FP,FN,total,no_calls,no_call_rate,sensitivity,specificity,PPV,NPV,F1,MCC); rows: roi,gene,whole_strategy)


                            methods = list(


                              # Load validated dataset results
                              loadValidatedResults = function(path, datasetName, bedFile) {
                                cat(paste("\nLoading", datasetName, "Dataset validated results\n"))

                                # load data depending on format
                                cols <- colnames(read.csv(path, sep="\t", header=TRUE, stringsAsFactors = FALSE))
                                if (all(cols[1:3] == c("SampleID", "Indication", "Truth"))) {  # panelcnDataset format
                                  data <- read.csv(path, sep="\t", header=TRUE, stringsAsFactors = FALSE, nrows = 176)
                                  data <- subset(data, SampleID != "")[, 1:5]
                                  cascades <- auxReadCascades(path)
                                  data <- transformToICRformat(data, cascades)
                                } else   # ICR96 format
                                  data <- read.csv(path, sep="\t", header=TRUE, stringsAsFactors = FALSE)

                                # add dataset
                                datasets[[datasetName]] <<- getSamplesWithResults(data, bedFile)
                                globalStats[[datasetName]] <<- list()
                              },


                              # Returns cascades from panecln dataset format file
                              auxReadCascades = function(path){
                                cascades <- read.csv(path, sep="\t", stringsAsFactors = FALSE, header = TRUE, skip = 191)
                                cascades <- cascades[!cascades[,1] %in% c(" ", "", NA),]
                                cascades <- cascades[cascades[,2] != "" & cascades[,2] != " ",]
                                row.names(cascades) <- cascades[,1]

                                return(cascades)
                              },


                              # transforms panecln format data frame to ICR format data frame
                              transformToICRformat = function(data, cascades) {
                                # Calculate exon positions
                                allCheckedGenes <- c()
                                for (i in 1:nrow(cascades)) {
                                  cas <- cascades[i, "Genes"]
                                  cas <- strsplit(cas, ", ")[[1]]
                                  allCheckedGenes <- unique(c(allCheckedGenes, cas))
                                }
                                exonsPositionsByGene <- getExonsPositions(allCheckedGenes)

                                # Create new data frame
                                newData <- data.frame(matrix(ncol = 6, nrow = 0))
                                colnames(newData) <- c("SampleID", "Chromosome", "X5PrimeExon37", "X3PrimeExon37", "ExonCNVType", "Gene")

                                # Fill new data frame
                                sampleIds <- unique(data["SampleID"])[[1]]
                                for(id in sampleIds) {

                                  rows <- subset(data, SampleID == id)
                                  currentCascade <- cascades[rows[1, "Indication"], "Genes"]
                                  currentGenes <- strsplit(currentCascade, ", ")[[1]]

                                  for (gene in currentGenes){
                                    exonsPositions <- exonsPositionsByGene[[gene]]

                                    # Get exon positions
                                    pos <- which(rows$Truth == gene)
                                    if (!identical(pos, integer(0))){
                                      exons <- rows[pos,]$Exons
                                      if (grepl("-", exons)){ # multi-exon
                                        exonStart <- strsplit(exons, "-")[[1]][1]
                                        exonEnd <- strsplit(exons, "-")[[1]][2]
                                      } else { # single-exon
                                        exonStart <- exons
                                        exonEnd <- exons
                                      }

                                      # Get start end position from exons info
                                      exonStart <- as.numeric(substr(exonStart,2, nchar(exonStart)))
                                      exonEnd <- as.numeric(substr(exonEnd,2, nchar(exonEnd)))
                                    } else {
                                      exonStart <- min(exonsPositions$rank)
                                      exonEnd <- max(exonsPositions$rank)
                                    }

                                    # Fix exon position in case it was not coding
                                    while(!exonStart %in% exonsPositions$rank) {
                                      exonStart <- exonStart + 1
                                      if (exonStart > 100)
                                        stop(paste("Error: exon not found for gene", gene))
                                    }
                                    while(!exonEnd %in% exonsPositions$rank){
                                      exonEnd <- exonEnd - 1
                                      if (exonEnd < 0)
                                        stop(paste("Error: exon not found for gene", gene))
                                    }

                                    # get start end coords
                                    if (unique(exonsPositions$strand) == "1"){
                                      start <- subset(exonsPositions, rank==exonStart)$genomic_coding_start
                                      end <- subset(exonsPositions, rank==exonEnd)$genomic_coding_end
                                    } else if (unique(exonsPositions$strand) == "-1"){
                                      start <- subset(exonsPositions, rank==exonEnd)$genomic_coding_start
                                      end <- subset(exonsPositions, rank==exonStart)$genomic_coding_end
                                    } else
                                      stop("Strand is not - or +.")

                                    # get chr
                                    chr <- unique(exonsPositions$chromosome_name)

                                    # set cnvType
                                    if (identical(pos, integer(0))){
                                      cnvType <- ""
                                    } else if (rows[pos,]$MLPAResult %in% c("CN0", "CN1")) {
                                      cnvType <- "Deletion"
                                    } else if (rows[pos,]$MLPAResult %in% c("CN3", "CN4"))
                                      cnvType <- "Duplication"

                                    row <- data.frame(SampleID = as.character(id), Chromosome = chr, X5PrimeExon37 = start,
                                                      X3PrimeExon37 = end, ExonCNVType = as.character(cnvType), Gene = gene)
                                    newData <- rbind(newData, row)
                                  }
                                }

                                i <- sapply(newData, is.factor)
                                newData[i] <- lapply(newData[i], as.character)
                                return(newData)
                              },


                              # Returns a list of Samples with validated results
                              getSamplesWithResults = function(data, bedFile){
                                samples <- list()

                                # load bed file regions
                                bedData <- read.csv(bedFile, header = T, sep = "\t", stringsAsFactors=F)
                                bedDataGR <- GRanges(seqnames = bedData[,1], ranges = IRanges(start=bedData[,2], end=bedData[,3]))

                                # extract sample names
                                sampleIds <- unique(data["SampleID"])[[1]]

                                # build Sample object for each sample name
                                for(id in sampleIds) {
                                  validatedGenes <- subset(data, SampleID == id)[["Gene"]]

                                  allPositives <- data.frame()
                                  allNegatives <- data.frame()
                                  for (gene in (validatedGenes)){

                                    allGenePositives <- data.frame()
                                    allGeneNegatives <- data.frame()
                                    for(cnvType in c("Deletion", "Duplication", "")) {

                                      # get matching data
                                      events <- subset(data, SampleID == id & ExonCNVType == cnvType & Gene == gene)
                                      grEvents <- GRanges(seqnames = events[["Chromosome"]], ranges = IRanges(start=events[["X5PrimeExon37"]], end=events[["X3PrimeExon37"]]))
                                      overlaps <- countOverlaps(bedDataGR, grEvents, type = "any")
                                      aux <- bedData[overlaps != 0,]

                                      if (nrow(aux) > 0) {
                                        aux$cnvType <- cnvType
                                      } else
                                        aux$cnvType <- character(0)

                                      # store it
                                      if (cnvType %in% c("Deletion", "Duplication")) {
                                        allGenePositives <- rbind(allGenePositives, aux)
                                      } else
                                        allGeneNegatives <- rbind(allGeneNegatives, aux)
                                    }

                                    # if postive was found, other ROIs belonging to the gene without CNV have to be included as negative
                                    if (nrow(allGenePositives) > 0) {
                                      posFoundsGR <- GRanges(seqnames = allGenePositives[[1]], ranges = IRanges(start=allGenePositives[[2]], end=allGenePositives[[3]]))
                                      overlaps <- countOverlaps(bedDataGR, posFoundsGR, type = "any")
                                      aux <- bedData[overlaps == 0 & bedData[,4] == gene,]
                                      if (nrow(aux) > 0){
                                        aux$cnvType <- ""
                                        allGeneNegatives <- rbind(allGeneNegatives, aux)
                                      }
                                    }

                                    allPositives <- rbind(allPositives, allGenePositives)
                                    allNegatives <- rbind(allNegatives, allGeneNegatives)
                                  }

                                  colnames(allPositives) <- c("chr", "start", "end", "gene", "cnvType")
                                  colnames(allNegatives) <- colnames(allPositives)
                                  allPositives$cnvType <- gsub("Deletion", "deletion", allPositives$cnvType)
                                  allPositives$cnvType <- gsub("Duplication", "duplication", allPositives$cnvType)
                                  allNegatives$cnvType <- gsub("", "none", allNegatives$cnvType)
                                  all <- rbind(allPositives, allNegatives)

                                  # sort by chr and start
                                  all <- all[with(all, order(chr, start)), ]


                                  # Build and store Sample object
                                  s <- SampleStats()
                                  s$new(id = as.character((id)),
                                        validatedROIs = all,
                                        validatedGenes = validatedGenes)
                                  samples[[as.character((id))]] <- s
                                }

                                return(samples)
                              },



                              # Loads algorithm results in GRanges format from outputFolder
                              # datasetName: dataset name for which algorithm results has been calculated
                              loadAlgorithmResults = function(outputFolder, algorithmName, datasetName, bedFile){

                                cat(paste("  Loading", algorithmName, "results for", datasetName, "dataset\n"))

                                # Load ROIs and filter them (only samples validated by MLPA)
                                if(file.exists(file.path(outputFolder, "failedROIs.csv"))){
                                  failedROIsData <- read.csv(file.path(outputFolder, "failedROIs.csv"), header = T, sep = "\t", stringsAsFactors=F)
                                  # failedROIsData[,1] <- unlist(lapply(failedROIsData[,1], getSampleName)) # process sample names
                                  failedROIsData <- failedROIsData[failedROIsData$SampleID %in%  getSamplesNames(datasetName),]
                                } else {
                                  failedROIsData <- data.frame(matrix(ncol = 5, nrow = 0))
                                  names(failedROIsData) <- c("SampleID", "Chr", "Start", "End", "Gene")
                                }

                                # load bed file regions
                                bedData <- read.csv(bedFile, header = T, sep = "\t", stringsAsFactors=F)
                                bedDataGR <- GRanges(seqnames = bedData[,1], ranges = IRanges(start=bedData[,2], end=bedData[,3]))

                                # Load results
                                samples <- datasets[[datasetName]]
                                resultsFile <- file.path(outputFolder, "grPositives.rds")
                                grPositives <- readRDS(resultsFile)

                                # Fix failedRois: translate them into ROIs of bed file
                                if (nrow(failedROIsData) > 0 ){
                                  fixedFailedRois <- data.frame()
                                  failedROIsGR <- GRanges(seqnames = failedROIsData$Chr, ranges = IRanges(start=failedROIsData$Start, end=failedROIsData$End))
                                  mcols(failedROIsGR)$sample <- failedROIsData$SampleID
                                  for (sample in unique(failedROIsData$SampleID)){
                                    failed <- failedROIsGR[mcols(failedROIsGR)[,"sample"] == sample ]
                                    overlaps <- countOverlaps(bedDataGR, failed, type = "any")
                                    aux <- bedData[overlaps != 0,]
                                    if (nrow(aux) > 0 ){
                                      aux$sample <- sample
                                      fixedFailedRois <- rbind(fixedFailedRois, aux)
                                    }
                                  }
                                  failedROIsData2 <- data.frame(fixedFailedRois$sample, fixedFailedRois[, 1],  fixedFailedRois[, 2], fixedFailedRois[, 3], fixedFailedRois[, 4])
                                  names(failedROIsData2) <- names(failedROIsData)
                                  failedROIsData <- failedROIsData2
                                }

                                # Fix found positives: translate them into ROIs of bed file
                                fixedPositives <- data.frame()
                                for (sample in unique(mcols(grPositives)[, "sample"])){

                                  for (cnvType in c("deletion", "duplication")){
                                    positives <- grPositives[mcols(grPositives)[,"sample"] == sample & mcols(grPositives)[,"cnvType"] == cnvType]
                                    overlaps <- countOverlaps(bedDataGR, positives, type = "any")
                                    aux <- bedData[overlaps != 0,]
                                    if (nrow(aux) > 0 ){
                                      aux$sample <- sample
                                      aux$cnvType <- cnvType
                                      fixedPositives <- rbind(fixedPositives, aux)
                                    }
                                  }
                                }
                                colnames(fixedPositives) <- c("chr", "start", "end", "gene", "sample", "cnvType")

                                # extract sample names
                                sampleIds <- unique(fixedPositives$sample)


                                # build Sample object for each sample name
                                for (s in samples){
                                  foundPositives <- data.frame()

                                  # Extract positives. Look only in those genes checked by validated technique and belonging to current sample
                                  if (s$id %in% sampleIds)
                                    foundPositives <- fixedPositives[fixedPositives$gene %in% s$validatedGenes & fixedPositives$sample == s$id,]

                                  # get ROIs only for validated genes and current sample
                                  failedSampleROIs <- failedROIsData[failedROIsData$Gene %in% s$validatedGenes & failedROIsData$SampleID == s$id,]

                                  # mark positives labeled as duplication and deletion as failed
                                  if (nrow(foundPositives) > 1) {
                                    wrong <- duplicated(foundPositives[, c("chr", "start")])
                                    if (TRUE %in% wrong){
                                      wrongRows <- foundPositives[wrong,]
                                      for (i in 1:nrow(wrongRows)){
                                        r <- wrongRows[i,]
                                        foundPositives <- foundPositives[!(foundPositives$chr == r$chr & foundPositives$start == r$start), ]

                                        if (nrow(failedSampleROIs[failedSampleROIs$Chr == r$chr & failedSampleROIs$Start == r$start, ]) == 0)
                                          failedSampleROIs[nrow(failedSampleROIs) + 1,] = c(s$id, r$chr, r$start, r$end, r$gene)

                                        print(paste("ROI in chr", r$chr, "and starting at", r$start, " for sample", s$id, "was labeled as failed region because was defined by the algorithm as deletion and duplication."))
                                      }
                                    }
                                  }

                                  # check found results against validated results
                                  s$loadResults(algorithmName, foundPositives, failedSampleROIs)
                                }

                                # add algorithm
                                algorithms <<- c(algorithms, algorithmName)
                                fillGlobalStats(datasetName, algorithmName)
                              },


                              # Writes summary report to file
                              writeSummary = function(path, datasetName){
                                cat(paste("  Generating summary file at", path, "\n"))

                                # Create summary folder if not exists
                                outputFolder <- dirname(path)
                                dir.create(outputFolder, showWarnings = FALSE)

                                # Create and fill summary output file
                                summaryText <- ""
                                for (alg in algorithms)
                                  summaryText <- .self$updateSummaryText(summaryText, datasetName, alg)
                                summaryText <- paste(summaryText, "\n\nSummary file generated at", Sys.time())
                                write(summaryText, path)
                              },


                              # Writes csv results file for a desired dataset (gene level)
                              writeFalseNegativeResults = function(path, datasetName){
                                cat(paste("  Generating False Negatives summary file at", path, "\n"))

                                allFns <- data.frame(sample = character(), gene = character(), cnvType = character(), alg = character(), stringsAsFactors=FALSE)
                                samples <- datasets[[datasetName]]
                                for (algorithm in names(globalStats[[datasetName]])){
                                  for(s in samples) {
                                    if (s$algorithms[[algorithm]]$nFNgene > 0 ){
                                      fns <- s$algorithms[[algorithm]]$FNdata
                                      for (i in 1:nrow(fns)) {
                                        allFns[nrow(allFns) + 1,] <- list(s$id, fns[i, "gene"], fns[i, "cnvType"], algorithm)
                                      }
                                    }
                                  }
                                }

                                write.table(allFns, path, sep="\t", row.names=FALSE, quote = FALSE)  # write output file
                              },


                              # Writes csv results file for a desired dataset
                              writeCSVresults = function(path, datasetName){
                                cat(paste("  Generating CSV results file at", path, "\n"))

                                # Create folder if not exists
                                outputFolder <- dirname(path)
                                dir.create(outputFolder, showWarnings = FALSE)

                                # Create and fill summary output file
                                all <- NULL
                                for (alg in names(globalStats[[datasetName]])){
                                  algResults <- globalStats[[datasetName]][[alg]]
                                  algResults$algorithm <- alg
                                  algResults$level <- rownames(algResults)

                                  if (is.null(all))
                                    all <- algResults
                                  else
                                    all <- rbind(all, algResults)
                                }


                                write.table(all, path, sep="\t", row.names=FALSE, quote = FALSE)  # write output file
                              },



                              # Fill global stats for desired dataset and algortithm
                              fillGlobalStats = function(datasetName, alg){
                                # create data.frame and fill with 0
                                df <- data.frame(matrix(ncol = 15, nrow = 4))
                                colnames(df) <- c("TP","TN","FP","FN","total","no_calls","no_call_rate","sensitivity","specificity","PPV","NPV","F1","MCC","kappa","accuracy")
                                rownames(df) <- c("roi","roi_incNoCall", "gene","whole_strategy") # stats level
                                for (i in 1:nrow(df))
                                  df[i,] <- rep(0, 15)

                                # Calculate basic stats
                                samples <- datasets[[datasetName]]
                                for(s in samples) {
                                  df["roi", "TP"] <- df["roi", "TP"] + s$algorithms[[alg]]$nTP
                                  df["roi", "TN"] <- df["roi", "TN"] + s$algorithms[[alg]]$nTN
                                  df["roi", "FP"] <- df["roi", "FP"] + s$algorithms[[alg]]$nFP
                                  df["roi", "FN"] <- df["roi", "FN"] + s$algorithms[[alg]]$nFN
                                  df["roi", "total"] <- df["roi", "total"] + s$algorithms[[alg]]$nROIs
                                  df["roi", "no_calls"] <- df["roi", "no_calls"] + s$algorithms[[alg]]$nFailedROIs
                                  df["roi_incNoCall", "TP"] <- df["roi_incNoCall", "TP"] + s$algorithms[[alg]]$nTP_incNoCall
                                  df["roi_incNoCall", "TN"] <- df["roi_incNoCall", "TN"] + s$algorithms[[alg]]$nTN_incNoCall
                                  df["roi_incNoCall", "FP"] <- df["roi_incNoCall", "FP"] + s$algorithms[[alg]]$nFP_incNoCall
                                  df["roi_incNoCall", "FN"] <- df["roi_incNoCall", "FN"] + s$algorithms[[alg]]$nFN_incNoCall
                                  df["roi_incNoCall", "total"] <- df["roi", "total"] + s$algorithms[[alg]]$nROIs
                                  df["roi_incNoCall", "no_calls"] <- df["roi", "no_calls"] + s$algorithms[[alg]]$nFailedROIs
                                  df["gene", "TP"] <- df["gene", "TP"] + s$algorithms[[alg]]$nTPgene
                                  df["gene", "TN"] <- df["gene", "TN"] + s$algorithms[[alg]]$nTNgene
                                  df["gene", "FP"] <- df["gene", "FP"] + s$algorithms[[alg]]$nFPgene
                                  df["gene", "FN"] <- df["gene", "FN"] + s$algorithms[[alg]]$nFNgene
                                  df["gene", "total"] <- df["gene", "total"] + s$algorithms[[alg]]$nCheckedGenes
                                  df["gene", "no_calls"] <- df["gene", "no_calls"] + s$algorithms[[alg]]$nFailedGenes
                                  df["whole_strategy", "TP"] <- df["whole_strategy", "TP"] + s$algorithms[[alg]]$nTPgene_incNoCall
                                  df["whole_strategy", "TN"] <- df["whole_strategy", "TN"] + s$algorithms[[alg]]$nTNgene_incNoCall
                                  df["whole_strategy", "FP"] <- df["whole_strategy", "FP"] + s$algorithms[[alg]]$nFPgene_incNoCall
                                  df["whole_strategy", "FN"] <- df["whole_strategy", "FN"] + s$algorithms[[alg]]$nFNgene_incNoCall
                                  df["whole_strategy", "total"] <- df["whole_strategy", "total"] + s$algorithms[[alg]]$nCheckedGenes
                                  df["whole_strategy", "no_calls"] <- df["whole_strategy", "no_calls"] + s$algorithms[[alg]]$nFailedGenes
                                }


                                # Calculate other stats
                                for (row in rownames(df)){
                                  df[row, "no_call_rate"] <- round(df[row, "no_calls"] / df[row, "total"] * 100, 4)
                                  df[row, "sensitivity"] <- round(df[row, "TP"] / (df[row, "TP"] + df[row, "FN"]), 4)
                                  df[row, "specificity"] <- round(df[row, "TN"] / (df[row, "TN"] + df[row, "FP"]), 4)
                                  df[row, "PPV"] <- round(df[row, "TP"] / (df[row, "TP"] + df[row, "FP"]), 4)
                                  df[row, "NPV"] <- round(df[row, "TN"] / (df[row, "TN"] + df[row, "FN"]), 4)
                                  df[row, "F1"] <- round(2 * df[row, "TP"] / (2 * df[row, "TP"] + df[row, "FP"] + df[row, "FN"]), 4)
                                  df[row, "MCC"] <- round((df[row, "TP"] * df[row, "TN"] - df[row, "FP"] * df[row, "FN"]) /
                                                            sqrt((df[row, "TP"] + df[row, "FP"]) * (df[row, "TP"] + df[row, "FN"]) * (df[row, "TN"] + df[row, "FP"]) * (df[row, "TN"] + df[row, "FN"])), 4)
                                  df[row, "accuracy"] <- round((df[row, "TP"] + df[row, "TN"]) / (df[row, "TN"] + df[row, "FP"] + df[row, "FN"] + df[row, "TP"]), 4)
                                  expectedAccuracy <- ((df[row, "TP"] + df[row, "FP"]) * (df[row, "TP"] + df[row, "FN"]) + (df[row, "TN"] + df[row, "FP"]) * (df[row, "TN"] + df[row, "FN"])) / ((df[row, "TN"] + df[row, "FP"] + df[row, "TP"] + df[row, "FN"])^2)
                                  df[row, "kappa"] <-   round((df[row, "accuracy"] - expectedAccuracy)/(1 - expectedAccuracy), 4)
                                }

                                # save stats
                                globalStats[[datasetName]][[alg]] <<- df
                              },


                              getSamplesNames = function(datasetName){
                                samples <- datasets[[datasetName]]
                                result <- c()

                                for(s in samples) {
                                  result <- append(result, s$id)
                                }
                                return(result)
                              },


                              # Prints false negatives for all samples of last checked results
                              printFalseNegatives_geneLevel = function(datasetName, algorithm, printable = TRUE){
                                samples <- datasets[[datasetName]]
                                text <- "\n\n\tFalse Negatives per gene:"

                                for(s in samples) {
                                  if (s$algorithms[[algorithm]]$nFNgene > 0 ){
                                    text <- paste(text, "\n\t\tSample", s$id, s$algorithms[[algorithm]]$nFNgene, "false negative(s)")
                                  }
                                }

                                if (printable)
                                  cat(text)
                                else
                                  return(text)
                              },


                              # Prints true negatives for all samples of last checked results
                              printTrueNegatives_geneLevel = function(datasetName, algorithm, printable = TRUE){
                                samples <- datasets[[datasetName]]
                                text <- "\n\n\tTrue Negatives per gene:"

                                for(s in samples) {
                                  if (s$algorithms[[algorithm]]$nTNgene > 0 )
                                    text <- paste(text, "\n\t\tSample", s$id, s$algorithms[[algorithm]]$nTNgene, "true negative(s)")
                                }

                                if (printable)
                                  cat(text)
                                else
                                  return(text)
                              },


                              # Prints false positives for all samples of last checked results
                              printFalsePositives_geneLevel = function(datasetName, algorithm, printable = TRUE){
                                samples <- datasets[[datasetName]]
                                text <- "\n\n\tFalse Positives per gene:"

                                for(s in samples) {
                                  if (s$algorithms[[algorithm]]$nFPgene > 0 )
                                    text <- paste(text, "\n\t\tSample", s$id, s$algorithms[[algorithm]]$nFPgene, "false positives(s)")
                                }

                                if (printable)
                                  cat(text)
                                else
                                  return(text)
                              },


                              # Prints false negatives for all samples (at gene level and including no-calls)
                              printWholeStrategyFalseNegatives = function(datasetName, algorithm, printable = TRUE){
                                samples <- datasets[[datasetName]]
                                text <- "\n\n\tFalse Negatives at whole strategy level:"

                                for(s in samples) {
                                  if (s$algorithms[[algorithm]]$nFNgene_incNoCall > 0 )
                                    text <- paste(text, "\n\t\tSample", s$id, s$algorithms[[algorithm]]$nFNgene_incNoCall, "false negative(s)")
                                }

                                if (printable)
                                  cat(text)
                                else
                                  return(text)
                              },


                              # Prints false positives for all samples (at gene level and including no-calls)
                              printWholeStrategyFalsePositives = function(datasetName, algorithm, printable = TRUE){
                                samples <- datasets[[datasetName]]
                                text <- "\n\n\tFalse Positives at whole strategy level:"

                                for(s in samples) {
                                  if (s$algorithms[[algorithm]]$nFPgene_incNoCall > 0 )
                                    text <- paste(text, "\n\t\tSample", s$id, s$algorithms[[algorithm]]$nFPgene_incNoCall, "false positives(s)")
                                }

                                if (printable)
                                  cat(text)
                                else
                                  return(text)
                              },


                              updateSummaryText = function(summaryText, datasetName, alg){
                                if (!identical(summaryText, character(0)))
                                  summaryText <- paste0(summaryText, "\n\n")
                                stats <- globalStats[[datasetName]][[alg]]

                                cText <- paste(alg,"results on dataset", datasetName,
                                               "\n\tStats per ROI:",
                                               "\n\t\tSensitivity:", stats["roi", "sensitivity"],
                                               "\n\t\tSpecificity:", stats["roi", "specificity"],
                                               "\n\t\tNo-call rate:", stats["roi", "no_call_rate"], "%  from a total of", stats["roi", "total"], "ROIs",
                                               "\n\tStats per gene:",
                                               "\n\t\tSensitivity:", stats["gene", "sensitivity"],
                                               "\n\t\tSpecificity:", stats["gene", "specificity"],
                                               "\n\t\tNo-call rate:", stats["gene", "no_call_rate"], "%  from a total of", stats["gene", "total"], "genes",
                                               "\n\t\tF1 Score:", stats["gene", "F1"],
                                               "\n\t\tPositive Predictive Value:", stats["gene", "PPV"],
                                               "\n\t\tNegative Predictive Value:", stats["gene", "NPV"],
                                               "\n\t\tMatthews correlation coefficient:", stats["gene", "MCC"],
                                               "\n\t\tCohen's kappa:", stats["gene", "kappa"],
                                               "\n\tStats at whole strategy level (per gene and including no-calls):",
                                               "\n\t\tSensitivity:", stats["whole_strategy", "sensitivity"],
                                               "\n\t\tSpecificity:", stats["whole_strategy", "specificity"],
                                               " -> ", stats["whole_strategy", "FP"], "extra gene(s) have to be covered by MLPA,", stats["whole_strategy", "TN"], "genes saved")

                                cText <- paste(cText, printFalseNegatives_geneLevel(datasetName, alg, FALSE))
                                cText <- paste(cText, printFalsePositives_geneLevel(datasetName, alg, FALSE))
                                cText <- paste(cText, printWholeStrategyFalseNegatives(datasetName, alg, FALSE))
                                cText <- paste(cText, printWholeStrategyFalsePositives(datasetName, alg, FALSE))

                                summaryText <- paste0(summaryText, cText)
                                return(summaryText)
                              },


                              # Returns desired stat metric. Ex: getStat("icr96", "decon", "sensitivity", "roi")
                              getStat = function(datasetName, alg, metric, level){
                                return(globalStats[[datasetName]][[alg]][level, metric])
                              }

                            ))




#####  SampleStats class definition. Represents sample stats for a certain dataset  #####

SampleStats <- setRefClass("SampleStats",

                           fields = list(id = "character",
                                         validatedROIs = "ANY", # GRanges object containing validated events per ROI: deletion, duplication or none
                                         validatedGenes = "character", # genes that have been validated by MLPA for this sample
                                         algorithms = "list", # list of algorithms stats for this sample (SampleAlgStats object)
                                         nROIs = "integer"),  # number of ROIs validated for this sample



                           methods = list(

                             # To be called after constructed. Sets initial values
                             new = function(id, validatedROIs, validatedGenes){
                               id <<- id
                               validatedROIs <<- validatedROIs
                               validatedGenes <<- validatedGenes
                               nROIs <<- nrow(validatedROIs)
                             },


                             # Check found positives against validated positivies. Sets stats vars
                             loadResults = function(algorithm, foundPositives, failedROIs) {
                               sas <- SampleAlgStats()
                               algorithms[[algorithm]] <<- sas
                               algorithms[[algorithm]]$new(foundPositives, nROIs, failedROIs)

                               # Add alg column with alg info
                               validatedROIs$alg <<- NA
                               for(i in 1:nrow(validatedROIs)){
                                 row <- validatedROIs[i,]
                                 if (nrow(foundPositives) > 0) {
                                   algRoi <- subset(foundPositives, chr == row$chr & start == row$start)

                                   if (nrow(algRoi) > 0) {
                                     validatedROIs[i, "alg"] <<- algRoi$cnvType
                                   } else
                                     validatedROIs[i, "alg"] <<- "none"
                                 } else
                                   validatedROIs[i, "alg"] <<- "none"
                               }

                               # Calculate ROIs stats
                               validatedROIs$eval <<- NA
                               for(i in 1:nrow(validatedROIs)){
                                 row <- validatedROIs[i,]

                                 if ( (row$cnvType == "deletion" && row$alg == "deletion") || (row$cnvType == "duplication" && row$alg == "duplication") ){
                                   sas$nTP <- sas$nTP + 1
                                   validatedROIs[i, "eval"] <<- "tp"
                                 } else if (row$cnvType == "none" && row$alg == "none"){
                                   sas$nTN <- sas$nTN + 1
                                   validatedROIs[i, "eval"] <<- "tn"
                                 } else if ((row$cnvType == "none" && row$alg != "none")
                                            || (row$cnvType == "duplication" && row$alg == "deletion")
                                            || (row$cnvType == "deletion" && row$alg == "duplication")){
                                   sas$nFP <- sas$nFP + 1
                                   validatedROIs[i, "eval"] <<- "fp"
                                 } else if (row$cnvType != "none" && row$alg == "none"){
                                   sas$nFN <- sas$nFN + 1
                                   validatedROIs[i, "eval"] <<- "fn"
                                 }
                               }


                               # Calculate ROIs stats including no-calls
                               for(i in 1:nrow(validatedROIs)){
                                 row <- validatedROIs[i,]
                                 isFailedROI <- (nrow(subset(sas$failedROIs, Start == row$start)) > 0)

                                 if ( validatedROIs[i, "eval"] == "tp" | isFailedROI){
                                   sas$nTP_incNoCall <- sas$nTP_incNoCall + 1
                                 } else if (validatedROIs[i, "eval"] == "fn") {
                                   sas$nFN_incNoCall <- sas$nFN_incNoCall + 1
                                 } else if (validatedROIs[i, "eval"] == "fp" | isFailedROI){
                                   sas$nFP_incNoCall <- sas$nFP_incNoCall + 1
                                 } else {
                                   sas$nTN_incNoCall <- sas$nTN_incNoCall + 1
                                 }
                               }

                               # Calculate stats at gene level
                               for (g in validatedGenes){
                                 geneValidatedROIs <- subset(validatedROIs, gene == g)

                                 if ("tp" %in% geneValidatedROIs$eval){  # TP if at least one roi was found by alg
                                   sas$nTPgene <- sas$nTPgene + 1
                                 } else if ("fn" %in% geneValidatedROIs$eval){  # FN if there was a positive but was not found
                                   sas$nFNgene <- sas$nFNgene + 1
                                   type <- geneValidatedROIs[geneValidatedROIs$eval == "fn", "cnvType"][1]
                                   sas$FNdata[sas$nFNgene, ] <- list(g, type)
                                 } else if ("fp" %in% geneValidatedROIs$eval){  # FP if there was no positive but algorithm found it
                                   sas$nFPgene <- sas$nFPgene + 1
                                 } else  # TN
                                   sas$nTNgene <- sas$nTNgene + 1
                               }

                               # Calculate stats at gene level and including no-calls
                               for (g in validatedGenes){
                                 geneValidatedROIs <- subset(validatedROIs, gene == g)
                                 geneFailedROIs <- subset(sas$failedROIs, Gene == g)

                                 if ("tp" %in% geneValidatedROIs$eval |  # TP if at least one roi was found by alg
                                     nrow(subset(geneFailedROIs, Start %in% subset(geneValidatedROIs, eval == "fn")$start)) > 0 ){  # and include no-calls
                                   sas$nTPgene_incNoCall <- sas$nTPgene_incNoCall + 1
                                 } else if ("fn" %in% geneValidatedROIs$eval){  # FN if there was a positive but was not found
                                   sas$nFNgene_incNoCall <- sas$nFNgene_incNoCall + 1
                                 } else if ("fp" %in% geneValidatedROIs$eval |  # FP if there was no positive but algorithm found it
                                            nrow(subset(geneFailedROIs, Start %in% subset(geneValidatedROIs, eval == "tn")$start)) > 0 ){
                                   sas$nFPgene_incNoCall <- sas$nFPgene_incNoCall + 1
                                 } else  # TN
                                   sas$nTNgene_incNoCall <- sas$nTNgene_incNoCall + 1
                               }

                               sas$countCheckedGenes()
                             }

                           ))





#####  SampleAlgStats class definition. Represents sample stats for an algorithm checked against a certain dataset   #####

SampleAlgStats <- setRefClass("SampleAlgStats",

                              fields = list(# ROI stats
                                nTP = "numeric",
                                nTN = "numeric",
                                nFN = "numeric",
                                nFP = "numeric",
                                nTP_incNoCall = "numeric",
                                nTN_incNoCall = "numeric",
                                nFN_incNoCall = "numeric",
                                nFP_incNoCall = "numeric",

                                # gene stats
                                nTPgene = "numeric",
                                nFPgene = "numeric",
                                nTNgene = "numeric",
                                nFNgene = "numeric",
                                FNdata = "ANY", # data.frame with two columns: gene and cnvType. Contains false negatives for this samples at gene level

                                # gene stats including no calls
                                nTPgene_incNoCall = "numeric",
                                nFPgene_incNoCall = "numeric",
                                nTNgene_incNoCall = "numeric",
                                nFNgene_incNoCall = "numeric",

                                nCheckedGenes = "numeric", # number of gene that have been validated for this sample
                                nROIs = "numeric", # number of total ROIs analyzed for this sample (failed and not failed ROIs)
                                nFailedROIs = "numeric", # number of failed ROIs (ROIs where CNV call was not possible due to algortihm quality controls)
                                nFailedGenes = "numeric", # number of failed Genes (Genes with AT LEAST one failed ROI)

                                failedROIs = "ANY", # GRanges() containing all failed ROIs (ROIs where CNV call was not possible due to algortihm quality controls)
                                foundPositives = "ANY"), # GRanges() object containing CNVs found by algorithms


                              methods = list(

                                # To be called after constructed. Sets initial values
                                new = function(foundPositives, nROIs, fROIs){
                                  foundPositives <<- foundPositives
                                  nTP <<- 0; nFP <<- 0; nTN <<- 0; nFN <<- 0
                                  nTP_incNoCall <<- 0; nFP_incNoCall <<- 0; nTN_incNoCall <<- 0; nFN_incNoCall <<- 0
                                  nTPgene <<- 0; nFPgene <<- 0; nTNgene <<- 0; nFNgene <<- 0
                                  nTPgene_incNoCall <<- 0; nFPgene_incNoCall <<- 0; nTNgene_incNoCall <<- 0; nFNgene_incNoCall <<- 0
                                  
                                  nROIs <<- nROIs
                                  nFailedROIs <<- nrow(fROIs)
                                  nFailedGenes <<- length(unique(fROIs[, "Gene"]))
                                  failedROIs <<- fROIs
                                  FNdata <<- data.frame(gene=character(), cnvType=character(),stringsAsFactors=FALSE)
                                },
                                
                                
                                countCheckedGenes = function(){
                                  nCheckedGenes <<- nTPgene + nTNgene + nFPgene + nFNgene
                                }
                              ))
