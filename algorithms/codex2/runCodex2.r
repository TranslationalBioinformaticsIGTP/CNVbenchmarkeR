# Runs codex2 over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runCodex2.r [codex2_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(CODEX2))
source(if (basename(getwd()) == "optimizers") "../utils/segment_targeted.R" else "utils/segment_targeted.R") # Load utils functions
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

# translates del/dup to deletion/duplication
auxCNname <- function(x) {
  if (x == "del") return("deletion") 
  else if (x == "dup") return("duplication")
}

# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
    codex2ParamsFile <- args[1]
    datasetsParamsFile <- args[2]
} else {
    codex2ParamsFile <- "codex2Params.yaml"
    datasetsParamsFile <- "../../datasets.yaml"
}

#Load the parameters file  
params <- yaml.load_file(codex2ParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
print(paste("Params for this execution:", list(params)))

# extract codex2 params
codex2Folder <- file.path(params$codex2Folder)


# go over datasets and run codex2 for those which are active
for (name in names(datasets)) {
    dataset <- datasets[[name]]
    if (dataset$include){
        print(paste("Starting codex2 for", name, "dataset", sep=" "))
      
        # extract fields
        bamsDir <- file.path(dataset$bams_dir)
        bedFile <- file.path(dataset$bed_file)
        fastaFile <- file.path(dataset$fasta_file)

        # Create output folder
        if (!is.null(params$outputFolder)) {
          outputFolder <- params$outputFolder
        } else
          outputFolder <- file.path(getwd(), "output", paste0("codex2-", name))  
        if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
          unlink(outputFolder, recursive = TRUE);
          dir.create(outputFolder)
        }
        
        files <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)
        
        # Do pre-calc part of the algorithm
        if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
          
          # GET COVERAGE, GC, MAPP
          # get bam directories, read in bed file, get sample names
          chr <- 1
          sampname <- as.matrix(unlist(strsplit(files,"\\.bam")))
          bambedObj <- getbambed(bamdir = files,
                                 bedFile = bedFile,
                                 sampname = sampname,
                                 projectname = "projectname",
                                 chr)
          ref <- bambedObj$ref
          
          # get raw depth of coverage, gc content and mappability
          coverageObj <- getcoverage(bambedObj, mapqthres = 20)
          Y <- coverageObj$Y
          readlength <- coverageObj$readlength
          gc <- getgc(chr, ref)
          mapp <- getmapp(chr,ref)
          
          # loop
          ref.all <- ref; Y.all <- Y; gc.all <- gc; mapp.all <- mapp
          chr.all <- rep(chr,length=length(mapp))
          targ.chr <- unique(as.matrix(read.table(bedFile, sep = "\t")[,1]))
          for(chr in 2:23) {
            if(chr == 23)
              chr='X'
            
            if(!is.element(chr,targ.chr))
              next
            
            # get bam directories, read in bed file, get sample names
            bambedObj <- getbambed(bamdir = files,
                                   bedFile = bedFile,
                                   sampname = sampname,
                                   projectname = "projectname",
                                   chr)
            ref <- bambedObj$ref
            
            # get raw depth of coverage, gc content and mappability
            coverageObj <- getcoverage(bambedObj, mapqthres = 20)
            Y <- coverageObj$Y
            readlength <- coverageObj$readlength
            gc <- getgc(chr, ref)
            mapp <- getmapp(chr,ref)
            
            ref.all <- c(ref.all,bambedObj$ref)
            Y.all <- rbind(Y.all,coverageObj$Y)
            gc.all <- c(gc.all,gc)
            mapp.all <- c(mapp.all,mapp)
            chr.all <- c(chr.all,rep(chr,length=length(mapp)))
          }

          if (!is.null(params$execution) && params$execution == "onlyPrecalcPhase") {
            dir.create(params$precalcFolder)
            saveRDS(ref.all, file.path(params$precalcFolder, "ref.all.rds"))
            saveRDS(Y.all, file.path(params$precalcFolder, "Y.all.rds"))
            saveRDS(gc.all, file.path(params$precalcFolder, "gc.all.rds"))
            saveRDS(mapp.all, file.path(params$precalcFolder, "map.all.rds"))
            saveRDS(chr.all, file.path(params$precalcFolder, "chr.all.rds"))
            
            print(paste("codex2 (Only pre-calc phase) for", name, "dataset finished", sep=" "))
            cat("\n\n\n")
            quit()
          }
        } else {  # skipPrecalcPhase mode: read previous results
          print(paste("codex2 Skipping pre-calc phase for", name, "dataset finished", sep=" "))
          ref.all <- readRDS(file.path(params$precalcFolder, "ref.all.rds"))
          Y.all <- readRDS(file.path(params$precalcFolder, "Y.all.rds"))
          gc.all <- readRDS(file.path(params$precalcFolder, "gc.all.rds"))
          mapp.all <- readRDS(file.path(params$precalcFolder, "map.all.rds"))
          chr.all <- readRDS(file.path(params$precalcFolder, "chr.all.rds"))
        }
        
        # QUALITY CONTROL
        
        # prepare vars
        sampname <- as.matrix(unlist(strsplit(files,"\\.bam"))) # REMOVE
        gene.all <-  as.matrix(read.table(bedFile, head=F, sep='\t')[,4])
        Y <-  Y.all; ref <- ref.all; gc <- gc.all; mapp <- mapp.all; gene <- gene.all; chr <- chr.all
        qcObj <- qc(Y, sampname, chr, ref, mapp, gc, 
                    cov_thresh = c(params$cov_thresh_down, params$cov_thresh_up),
                    length_thresh = c(params$length_thresh_down, params$length_thresh_up),
                    mapp_thresh = params$mapp_thresh,
                    gc_thresh = c(params$gc_thresh_down, params$gc_thresh_up))
        Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc; gc_qc <- qcObj$gc_qc; mapp_qc <- qcObj$mapp_qc; ref_qc <- qcObj$ref_qc ; qcmat <- qcObj$qcmat
        gene_qc <- gene[which(as.logical(qcmat[,4])==TRUE)]
        chr_qc <- chr.all[which(as.logical(qcmat[,4])==TRUE)]

        # sometimes need to exclude samples with very few reads (due to capture failure)
        sampfilter <- apply(Y_qc, 2, median) >= params$sample_reads_median_limit
        sampname_qc <- sampname_qc[sampfilter]
        Y_qc <- Y_qc[,sampfilter]
        rm(qcObj)

        
        # NORMALIZATION
        
        normObj <- normalize2(Y_qc, gc_qc , K=1:15, normal_index=1:length(sampname)) 
        Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC; RSS <- normObj$RSS; K <- normObj$K
        
        
        # SEGMENTATION
        
        optK <- which.max(BIC)
        finalcall <- matrix(ncol=14)
        colnames(finalcall) <- c('sample_name','chr','gene','cnv','st_bp','ed_bp','length_kb','st_exon','ed_exon','raw_cov','norm_cov','copy_no','lratio','mBIC')

        for(genei in unique(gene_qc)){
          if (!is.na(genei)){
            cat('Segmenting gene', genei, '\n')
            geneindex <- which(gene_qc == genei)
            yi <- Y_qc[geneindex,]
            yhati <- Yhat[[optK]][geneindex,]
            refi <- ref_qc[geneindex]
            chri <- chr_qc[geneindex][1]
            finalcalli <- segment_targeted(yi, yhati, sampname_qc, refi, genei, chri, lmax=length(geneindex), mode='fraction')
            finalcall <- rbind(finalcall,finalcalli)
          }
        }

        finalcall <- finalcall[-1,]
        cn <- (as.numeric(as.matrix(finalcall[,'copy_no'])))
        cn.filter <- (cn <= params$cn_del_factor) | (cn >= params$cn_dup_factor)
        finalcall <- finalcall[cn.filter,]
        length_exon <- as.numeric(finalcall[,'ed_exon']) - as.numeric(finalcall[,'st_exon']) + 1
        finalcall <- cbind(finalcall[,1:7], length_exon, finalcall[,10:14])
        
        # set right sample name and cnv naming
        for (i in 1:nrow(finalcall)){
          parts <- strsplit(finalcall[i, "sample_name"], "/")[[1]]
          finalcall[i, "sample_name"] <- parts[length(parts)]
          finalcall[i, "cnv"] <- auxCNname(finalcall[i, "cnv"])
        }

        # save results
        write.table(finalcall, file = file.path(outputFolder, "cnvFounds.csv"), sep='\t', quote=F, row.names=F)
        
        # Save results in GRanges format
        message("Saving CNV GenomicRanges")
        saveResultsFileToGR(outputFolder, "cnvFounds.csv",
                            geneColumn = "gene",
                            sampleColumn = "sample_name", 
                            chrColumn = "chr",
                            startColumn = "st_bp",
                            endColumn = "ed_bp",
                            cnvTypeColumn = "cnv")
        
        print(paste("CODEX2 for", name, "dataset finished", sep=" "))
        cat("\n\n\n")
    }
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)