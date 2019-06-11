#USAGE: Rscript optimizer.r [optimizer_params_file]   (call from optimizers folder)
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
library(methods)
options(scipen = 999)  # to disable scientific number notation
source("../utils/cnvStats.r")  # Load utils functions
source("../utils/optimizerUtils.r")  # Load opt. utils functions


# Returns an inteval up to 9 Values (class object) containing default one
# value: object containing default, min, max.
getIntervalValues <- function(value, name, low = FALSE){
  
  if (low == TRUE) {
    print(paste("Low resolution for param", name))
    exp <- c(0.6, 0.85, 0.92, 0.97, 1, 1.03, 1.08, 1.15, 1.4)
  } else
    exp <- c(0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.02, 1.04, 1.06, 1.08, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75)

  defaultValue <- as.numeric(value$default)
  result <- list()
  for (e in exp){
    v <- abs(defaultValue) ^ e
    if (defaultValue < 0)
      v <- -v
    
    if (v >= value$min && v <= value$max){
      newVal <- Value()
      newVal$new(v)
      result[[length(result) + 1]] <- newVal
    }
  }
  
  values <- ValuesList()
  values$new(result, name, defaultValue, TRUE)
  
  return(values)
}

# Return available options for a categorical param
getFixedValues = function(param, name){
  vl <- ValuesList()
  vals <- list()
  for (v in param$options){
    newVal <- Value()
    newVal$new(v)
    vals[[length(vals) + 1]] <- newVal
  }
  vl$new(vals, name, param$default, FALSE)
  
  return(vl)
}


#Load the parameters file  
args <- commandArgs(TRUE)
print(args)

params <- yaml.load_file(args[1])
algName <- params$algorithm$name
algFolder <- file.path(getwd(), algName)

# create logs and input folder to be used by jobs
logsFolder <- file.path(algFolder, "logs")
unlink(logsFolder, recursive = TRUE);
dir.create(logsFolder, showWarnings = FALSE)
inputFolder <- file.path(algFolder, "input")
unlink(inputFolder, recursive = TRUE);
dir.create(inputFolder, showWarnings = FALSE)
outputFolder <- file.path(algFolder, "output")
unlink(outputFolder, recursive = TRUE);
dir.create(outputFolder, showWarnings = FALSE)

# Create SummaryStats object and load current dataset validated results
ss <- SummaryStats()
ss$loadValidatedResults(params$dataset[[1]]$validated_results_file, names(params$dataset) , params$dataset[[1]]$bed_file)


### ALGORITHM FIRST STEPS: PRECALC PHASE  ###

# Create precalc folder
precalcFolder <- file.path(algFolder, "precalc")
unlink(precalcFolder, recursive = TRUE);
dir.create(precalcFolder, showWarnings = FALSE)


# create dataset yaml file
datasetParamsFile <- file.path(precalcFolder, "datasetParams.yaml")  
datasetContent <- strsplit(readChar(args[1], file.info(args[1])$size), "dataset:")[[1]][2]
datasetContent <- gsub("\n  ", "\n", datasetContent)
write(datasetContent, file = datasetParamsFile)

# create algorithm yaml file
algParamsFile <- file.path(precalcFolder, "algorithmParams.yaml")  
content <- 'execution: "onlyPrecalcPhase"\n'
content <- paste0(content, 'precalcFolder: ', precalcFolder, '\n')
content <- paste0(content, 'outputFolder: ', precalcFolder, '\n')
if (!is.null(params$algorithm$convadingFolder))
  content <- paste0(content, 'convadingFolder: ', params$algorithm$convadingFolder, '\n')  
if (!is.null(params$algorithm$deconFolder))
  content <- paste0(content, 'deconFolder: ', params$algorithm$deconFolder, '\n')  
write(content, file = algParamsFile)


# call algorithm (only precalc phase)
logPrecalcFile <- file.path(precalcFolder, "precalc.log")
cmd <- paste("Rscript", params$algorithm$path, algParamsFile, datasetParamsFile, ">", logPrecalcFile, "2>&1")
print(cmd)
system(cmd)

# copy job.sh script to algorithm folder
script <- readChar("job.sh", file.info("job.sh")$size)
script <- gsub("\\*rootPath\\*", dirname(getwd()), script)
script <- gsub("\\*algorithmScript\\*", params$algorithm$path, script)
file <- file(file.path(algFolder, "job.sh"))
write(script, file)


#### THROW OPTIMIZE ALG: CLUSTER JOBS TO CALCULATE CNVs WITH DIFFERENT PARAMS ####

#Load algortihm parameters values  
algorithmParamsFile <- file.path(algFolder, paste0(algName, "Params.yaml"))
algParams <- yaml.load_file(algorithmParamsFile)
paramsValues <- list()  # list of values to be tested for each param
for (i in 1:length(algParams)){
  p <- algParams[i]
  paramName <- names(p)
  
  if (is.null(p[[paramName]]$options)){ # numerical value
    paramsValues[[i]] <- getIntervalValues(p[[paramName]], paramName)
  } else
    paramsValues[[i]] <- getFixedValues(p[[paramName]], paramName)
}


# call optimizer algorithm
optimize(params_vals = paramsValues,
         stats = ss,
         datasetName = names(params$dataset),
         datasetContent = datasetContent,
         optimizerParams = params$params,
         bedFile =  params$dataset[[1]]$bed_file,
         algFolder = algFolder,
         algOwnParams = params$algorithm)
