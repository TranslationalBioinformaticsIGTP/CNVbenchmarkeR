# Collection of utils 

# Stores tab results file to GenomicRanges file (named as "grPositives.rds")
saveResultsFileToGR = function(folderPath, resultsFileName, geneColumn = "Gene", sampleColumn = "Sample", 
                               chrColumn = "Chr", startColumn = "Start", endColumn = "End", cnvTypeColumn = "CNV.type") {
  suppressPackageStartupMessages(library(GenomicRanges))

  # load results file path
  resultsFilePath <- file.path(folderPath, resultsFileName)
  data <- read.csv(resultsFilePath, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  
  # extract sample names
  sampleIds <- unique(data[sampleColumn])[[1]]

  all <- GRanges()
  for(id in sampleIds) {

    # Extract positives
    # foundPositives <- subset(data, Sample == id)
    foundPositives <- data[data[, sampleColumn] == id, ]
  
    # Convert to GR
    grPositives = obtainGRpostives(foundPositives, id, chrColumn, startColumn, endColumn, geneColumn, cnvTypeColumn)
    
    # store
    all <- c(all, grPositives)
  }
  
  # save all results
  saveRDS(all, file.path(folderPath, "grPositives.rds"))
}


# Stores results folder to GenomicRanges file (named as "grPositives.rds"). Obly contained tab files with specified pattern name will be readed.
saveResultsFolderToGR = function(folderPath, pattern, geneColumn = "GENE", sampleColumn = "Sample", chrColumn = "Chr", 
                                 startColumn = "Start", endColumn = "End", cnvTypeColumn = "CNV.type") {
  suppressPackageStartupMessages(library(GenomicRanges))
  
  all <- GRanges()
  for(file in dir(folderPath, pattern=pattern)) {
    foundPositives <- read.csv(file.path(folderPath, file), sep="\t", header=TRUE, stringsAsFactors = FALSE)
    
    if (nrow(foundPositives) > 0){  

      # Convert positives to GR
      grPositives = obtainGRpostives(foundPositives, file, chrColumn, startColumn, endColumn, geneColumn, cnvTypeColumn)
      
      # store
      all <- c(all, grPositives)
    } 
  }
  
  
  # save all results
  saveRDS(all, file.path(folderPath, "grPositives.rds"))
}

# Returns dataframe with columns (SampleID, Genes). Only works form panelcnDataset style format
# Returns NULL if format was not valid
# file: MLPA results file
readIndicationsForSamples = function(file){

  # Read cascades
  cascades <- read.csv(file, sep="\t", stringsAsFactors = FALSE, header = TRUE, skip = 191)[,c(1,2)]
  cascades <- cascades[cascades[,1] != "" & cascades[,1] != " ",]
  cascades <- cascades[cascades[,2] != "" & cascades[,2] != " ",]
  row.names(cascades) <- cascades[,1]

  # Return NULL if not valid format
  if (!all(colnames(cascades) == c("X", "Genes")))
    return(NULL)
    
  # Read samples and indications
  sample_indications <- read.csv(file, sep="\t", stringsAsFactors = FALSE, header = TRUE, nrows = 176)[,c(1,2)]
  sample_indications <- sample_indications[sample_indications[,1] != "" & sample_indications[,2] != "" & sample_indications[,2] != " ",]  
  
  # Return NULL if not valid format
  if (!all(colnames(sample_indications) == c("SampleID", "Indication")))
    return(NULL)
  
  # Set columns names and replace indication
  colnames(sample_indications) <- c("SampleID", "Genes")
  for (i in 1:nrow(sample_indications)){
    sample_indications[i, "Genes"] <-  cascades[cascades$X == sample_indications[i,]$Genes, ]$Genes
  }
  
  return(sample_indications)
}


# returns list with indication description as key and compatible indications as value
# file: MLPA results file
readCompatibleIndications = function(file){
  cascades <- read.csv(file, sep="\t", stringsAsFactors = FALSE, header = TRUE, skip = 191)[,c(1,2)]
  cascades <- cascades[cascades[,1] != "" & cascades[,1] != " ",]
  cascades <- cascades[cascades[,2] != "" & cascades[,2] != " ",]
  indications <- unique(cascades[, "Genes"])
  result <- list()
  
  for(i in 1:length(indications)) {
    genes <- strsplit(indications[i], ", ")[[1]]
    compatibles <- c()
    for(j in 1:length(indications)) {
      genesToCheck <- strsplit(indications[j], ", ")[[1]]
      
      compatible <- T
      for (k in 1:length(genes)){
        if (genes[k] %in% genesToCheck){
          compatible <- F
          break
        }
      }
      
      if (compatible){
        if (!(indications[j] %in% compatibles))
          compatibles[length(compatibles) + 1] <- indications[j]
      }
    }
    result[[indications[i]]] <- compatibles
  }

  return(result)
}


# Returns exons positions for canonical transcript (biggest CDS) of grch37 genome. Returns a list whose index is gene name
# geneSymbols: vector of gene symbols
getExonsPositions = function(genesSymbols){
  suppressPackageStartupMessages(library("biomaRt"))
  human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2012.archive.ensembl.org")
  
  result <- list()
  for (gene in genesSymbols){
    data <- getBM(attributes=c('transcript_biotype','ensembl_transcript_id', 'external_transcript_id'),
                  filters=c("hgnc_symbol", 'with_ccds'),
                  value=list(gene, TRUE),
                  mart=human)
    
    # Filter: only finishing 001 - protein coding
    data <- data[data$transcript_biotype == "protein_coding",]
    if (nrow(data[grepl('001$', data$external_transcript_id),]) > 0) {
      data <- data[grepl('001$', data$external_transcript_id),]  # ends with 001
    } else
      data <- data[grepl('01$', data$external_transcript_id),]  # ends with 01
    
    if (nrow(data) == 0){
      stop(paste("Error when obtaining exons positions for gene", gene, ": no protein coding transcript name ending with 001 or 01 was found"))
    }
    
    data <- getBM(attributes = c("chromosome_name", "genomic_coding_start", "genomic_coding_end", "ensembl_transcript_id", "ensembl_exon_id", "rank", "strand"), 
                  filters = c("ensembl_transcript_id"),
                  value = list(data$ensembl_transcript_id),
                  mart=human)
    data <- data[!is.na(data$genomic_coding_start) & !is.na(data$genomic_coding_end),]
    if (nrow(data) == 0){
      stop(paste("Error after filtering fields when obtaining exons positions for gene", gene, ": no protein coding transcript name ending with 001 was found"))
    }

    result[[gene]] <- data
  }
  
  return(result)
}



# Returns gen positions (chr - start - end - gene) as data.frame
# sorted bedData data frame containeing all ROIs (chr - start - end - gene)
getGenePositions <- function(bedData){
  
  genePositions <- data.frame()
  for (gene in unique(bedData[,4])){
    geneData <- bedData[bedData[,4] == gene,]
    row <- data.frame(geneData[1, 1], geneData[1, 2], geneData[nrow(geneData), 3], gene)
    names(row) <- c("chr", "start", "end", "gene")
    
    # Add new row
    genePositions <- rbind(genePositions, row)
  }

  return(genePositions)  
}


### AUXILIARY FUNCTIONS ###


# Processes full name to extract sample name
getSampleName = function(fullName){
  extractedId <- NULL
  if (grepl("\\.", fullName)){
    extractedId <- strsplit(fullName, "\\.")[[1]][1]
  } else {
    extractedId <- fullName
  }
  return(extractedId)
}



# Returns GenomicRanges object with gene and sample metadata columns
obtainGRpostives = function(foundPositives, fullSampleName, chrColumn, startColumn, endColumn, geneColumn, cnvTypeColumn){
  # Convert to GRanges
  start <- foundPositives[[startColumn]]
  end <- foundPositives[[endColumn]]
  
  # check if any algorithm put start bigger than end (bug)
  idx <- which(start > end)
  if (length(idx) > 0) {
    # swap values
    aux <- start[idx]
    start[idx] <- end[idx]
    end[idx] <- aux
  }
    
  grPositives = GRanges(seqnames = foundPositives[[chrColumn]],
                        ranges = IRanges(start=start, end=end))
  
  # extract sample name
  extractedId <- getSampleName(fullSampleName)


  # Add metadata columns: gene, sample info and cnv type
  mcols(grPositives)$gene <- foundPositives[[geneColumn]]
  mcols(grPositives)$sample <- extractedId
  mcols(grPositives)$cnvType <- foundPositives[[cnvTypeColumn]]
  
  return(grPositives)
}