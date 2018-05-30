# Generates summary file
#USAGE: Rscript summary.r [algortihms_params_file] [datasets_params_file]
source(if (basename(getwd()) == "utils") "cnvStats.r" else "utils/cnvStats.r") # Load class definitions
suppressPackageStartupMessages(library(yaml))

# load algorithms and datasets from yaml files
args <- commandArgs(TRUE)
cat("\n"); print(args)
if(length(args)>0) {
  algorithmsParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  algorithmsParamsFile <- "algorithms.yaml"
  datasetsParamsFile <- "datasets.yaml"
}

#Load the parameters file
algorithms <- yaml.load_file(algorithmsParamsFile)$algorithms
datasets <- yaml.load_file(datasetsParamsFile)

# Run summary for all datasets and algorithms
for (dName in names(datasets)) {
  dataset <- datasets[[dName]]

  if (dataset$include){

    # Create SummaryStats object and load current dataset validated results
    ss <- SummaryStats()
    ss$loadValidatedResults(dataset$validated_results_file, dName, dataset$bed_file)  

    # Calculate metrics
    for (alg in names(algorithms)) {
      useAlg <- algorithms[[alg]]
      if (useAlg){
        outputFolder <- file.path(getwd(), "output", paste0(alg, "-", dName))
        ss$loadAlgorithmResults(outputFolder, alg, dName, dataset$bed_file)
      }
    }
    
    # Write results to summary file
    ss$writeSummary(file.path(getwd(), "output", "summary", paste0("summary-", dName, ".txt")), dName)
    ss$writeCSVresults(file.path(getwd(), "output", "summary", paste0("results-", dName, ".csv")), dName)
  }
}  