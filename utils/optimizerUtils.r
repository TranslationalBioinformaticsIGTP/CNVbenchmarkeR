suppressPackageStartupMessages(library(methods))


#####  Value class: represents a value of a param executed along with fixed values from the rest of params   #####
Value <- setRefClass("Value",
                     
                     fields = list(status = "character",  # done, executing, error, none
                                   value = "ANY",
                                   valuesSequence = "character",  # combination of all params values Ej: 5277867554
                                   
                                   # stats achieved after execution
                                   roi_sensitivity = "numeric", 
                                   gene_nocallsrate = "numeric",
                                   gene_sensitivity = "numeric",  
                                   gene_specificity = "numeric", 
                                   ws_sensitivity = "numeric",  
                                   ws_specificity = "numeric",  
                                   
                                   cmd = "character",  # command used to launch the sequenece
                                   jobId = "character",
                                   paramsContent = "character"), # character describing used params
                     
                     methods = list(
                       
                       # To be called after constructed. Sets initial values
                       new = function(val){
                         value <<- val
                         status <<- "none"
                       },
                       
                       # Returns a copy of current object
                       clone = function(){
                         cloned <- Value()
                         cloned$new(value)
                         cloned$status <- status
                         cloned$valuesSequence <- valuesSequence
                         cloned$roi_sensitivity <- roi_sensitivity
                         cloned$gene_nocallsrate <- gene_nocallsrate
                         cloned$gene_sensitivity <- gene_sensitivity
                         cloned$gene_specificity <- gene_specificity
                         cloned$ws_sensitivity <- ws_sensitivity                         
                         cloned$ws_specificity <- ws_specificity
                         cloned$cmd <- cmd
                         cloned$paramsContent <- paramsContent
                         cloned$jobId <- jobId
                         
                         return(cloned)
                       },
                       
                       # execute on cluster
                       # params(list): Key is param name, value is param value
                       launch = function(params, valuesSequenceIdxs, algFolder, datasetContent, optimizerParams, algOwnParams){
                         out <- system(paste0("qstat -u \"", optimizerParams$clusterUser, "\""), intern = TRUE)
                         nJobs <- length(out) - 2 
                         
                         if(nJobs < optimizerParams$clusterJobsLimit){ # maximum number of all user jobs
                           
                           # store combination of params idx to recognize this execution
                           valuesSequence <<- valuesSequenceIdxs
                           cat(paste("\nParams sequence to execute: ", valuesSequence, "\n"))
                           
                           # create inputs files / folders
                           inputFolder <- file.path(algFolder, "input", valuesSequence)
                           dir.create(inputFolder, showWarnings = FALSE)
                           outputFolder <- file.path(algFolder, "output", valuesSequence)
                           dir.create(outputFolder, showWarnings = FALSE)
                           precalcFolder <- file.path(algFolder, "precalc")
                           algParamsFile <- file.path(inputFolder, "algorithm.yaml")
                           datasetParamsFile <- file.path(inputFolder, "dataset.yaml")
                           write(datasetContent, file = datasetParamsFile)
                           
                           # create params file
                           content <- 'execution: "skipPrecalcPhase"\n'
                           content <- paste0(content, 'precalcFolder: ', precalcFolder, '\n')
                           content <- paste0(content, 'outputFolder: ', outputFolder, '\n')
                           for (name in names(algOwnParams)){
                             if (!name %in% c("name", "path")){
                               content <- paste0(content, name, ': ', algOwnParams[[name]], '\n')
                             }
                           }
                           paramsContent <<- ""
                           for (paramName in names(params)){
                             paramsContent <<- paste0(paramsContent, paramName, ": ", params[paramName], "\n")
                           }
                           content <- paste0(content, paramsContent)
                           write(content, file = algParamsFile)
                           
                           # create command to launch on cluster
                           cmd <<- "qsub"
                           cmd <<- paste0(cmd, " -v algorithmParams=", algParamsFile)
                           cmd <<- paste0(cmd," -v datasetParams=", datasetParamsFile)
                           logFile <- file.path(algFolder, "logs", paste0(valuesSequence, ".log"))
                           cmd <<- paste(cmd, "-e", logFile, "-o", logFile)
                           cmd <<- paste(cmd, file.path(algFolder,"job.sh"))
                           print(cmd)
                           
                           # launch command and store jobId
                           out <- system(cmd, intern = TRUE)
                           jobId <<- strsplit(out, " ")[[1]][3]

                           status <<- "executing"  
                           return(TRUE)
                         } else
                           return(FALSE)
                       },
                       
                       
                       # Checks if current val is better or equal (sensitiviity at three levels) than anotherVal
                       isBetter = function(anotherVal){
                         if (is.null(anotherVal))
                           return(TRUE)
                         
                         if (ws_sensitivity > anotherVal$ws_sensitivity){
                           return(TRUE)
                         } else if(ws_sensitivity == anotherVal$ws_sensitivity & gene_sensitivity > anotherVal$gene_sensitivity) {
                          return(TRUE)
                         } else if(ws_sensitivity == anotherVal$ws_sensitivity & gene_sensitivity == anotherVal$gene_sensitivity & roi_sensitivity > anotherVal$roi_sensitivity) {
                           return(TRUE)
                         } else {
                           return(FALSE)
                         }
                       },
                       
                       
                       # Sets a Value status as finished id sequence is already done in other execution
                       # Completes Value info with anotherValue
                       alreadyDone = function(anotherValue){
                         if (anotherValue$status == "done"| anotherValue$status == "error"){
                           cat(paste(anotherValue$valuesSequence,"sequence already executed\n\n"))
                           roi_sensitivity <<- anotherValue$roi_sensitivity
                           gene_nocallsrate <<- anotherValue$gene_nocallsrate
                           gene_sensitivity <<- anotherValue$gene_sensitivity
                           gene_specificity <<- anotherValue$gene_specificity
                           ws_sensitivity <<- anotherValue$ws_sensitivity                         
                           ws_specificity <<- anotherValue$ws_specificity
                           valuesSequence <<- anotherValue$valuesSequence
                           cmd <<- anotherValue$cmd
                           paramsContent <<- anotherValue$paramsContent
                           status <<- anotherValue$status
                         }
                       },
                       
                       
                       # Returns character to print the sequence
                       getStr = function(){
                         str <-  paste("Sequence", valuesSequence, "with", ws_sensitivity, "sensitivity and ", ws_specificity,
                                       "specificity and", gene_nocallsrate, "% no call-rate\nGene sensitivity: ", gene_sensitivity,
                                       "\nGene specificity: ", gene_specificity,
                                       "\nROI sensitivity: ", roi_sensitivity,
                                       "\nParams used:\n", paramsContent)
                         return(str)
                       },
                       
                       
                       # generates a yaml file with sequence value params
                       generateYaml = function(algFolder){
                         resultsFolder <- file.path(algFolder, "output", "results")
                         dir.create(resultsFolder, showWarnings = FALSE)
                         outputFile <- file.path(resultsFolder, paste0(valuesSequence, ".yaml"))
                         write(paramsContent, file = outputFile)
                       },
                       
                       
                       # Checks if execution is finished. Updates sensitivity and specificity if finished.
                       checkFinished = function(summaryStats, datasetName, algFolder, bedFile){
                         
                         out <- system(paste0("qstat -j ", jobId, " 2>&1"), intern = TRUE)
                         if (grepl("do not exist", out[1])){

                           # load alg results in SummaryStats object
                           outputFolder <- file.path(algFolder, "output", valuesSequence)
                           tryCatch({
                             summaryStats$loadAlgorithmResults(outputFolder, valuesSequence, datasetName, bedFile)
                             
                             # Calculate stats (at gene level including no calls)
                             roi_sensitivity <<- summaryStats$getStat(datasetName, valuesSequence, "sensitivity", "roi")
                             gene_nocallsrate <<- summaryStats$getStat(datasetName, valuesSequence, "no_call_rate", "gene")
                             gene_sensitivity <<- summaryStats$getStat(datasetName, valuesSequence, "sensitivity", "gene")
                             gene_specificity <<- summaryStats$getStat(datasetName, valuesSequence, "specificity", "gene")
                             ws_sensitivity <<- summaryStats$getStat(datasetName, valuesSequence, "sensitivity", "whole_strategy")
                             ws_specificity <<- summaryStats$getStat(datasetName, valuesSequence, "specificity" , "whole_strategy")
                             
                             status <<- "done"  
                           },error = function(error_condition) {
                             status <<- "error"
                             print(paste("Controlled error at sequence ", valuesSequence, "->", error_condition))
                             ws_sensitivity <<- -1
                             ws_specificity <<- -1
                             gene_sensitivity <<- -1
                             gene_specificity <<- -1
                             roi_sensitivity <<- -1
                             gene_nocallsrate <<- -1
                           })
 
                         }
                       }
                     )
)


##### Class representing a list of values to be executed for a certain param  ######
ValuesList <- setRefClass("ValuesList",
                          fields = list(listOfValues = "list",
                                        paramName = "character",
                                        fixed = "ANY",  # chosen value
                                        fixedIdx = "numeric",  # chosen value index
                                        originalDefault = "ANY",  # Original default value that is not going to be modified
                                        isDone = "logical",  # TRUE when all values have been executed
                                        isNumeric = "logical"),  # Numeric (TRUE) or character (FALSE)
                          
                          methods = list(
                            
                            # To be called after constructed. Sets initial values
                            new = function(listValues, param_name, fixedV, is_numeric){
                              listOfValues <<- listValues
                              paramName <<- param_name
                              fixed <<- fixedV
                              originalDefault <<- fixedV
                              isDone <<- FALSE
                              isNumeric <<- is_numeric

                              for (i in 1:length(listOfValues)){
                                if (listOfValues[[i]]$value == fixed){
                                  fixedIdx <<- i
                                  break
                                }
                              }
                            },
                            
                            
                            # Returns Value index in list of Values. False if not found
                            find = function(val, list){
                              for (i in 1:length(list))
                                if (list[[i]]$value == val$value)
                                  return(i)
                              
                              return(FALSE)
                            },

                            
                            # Updates fixed value taking best of existing values
                            fixBestValue = function(minSpecificity, optimizerParams){
                              
                              bests <- filterAndSortBestValues(listOfValues, minSpecificity, fixedIdx, optimizerParams)
                              
                              # select default value if all have same specificity
                              if (bests[[1]]$ws_specificity == bests[[length(bests)]]$ws_specificity){
                                defaultIdx <- find(listOfValues[[fixedIdx]], bests)
                                if (defaultIdx != FALSE) {
                                  myList <- list()
                                  myList[[1]] <- bests[[defaultIdx]]
                                  bests <- myList
                                }
                              }
                              
                              fixedIdx <<- find(bests[[1]], listOfValues)
                              fixed <<- listOfValues[[fixedIdx]]$value
                            }
                          )
)


StartPoint <- setRefClass("StartPoint",
                          fields = list(paramIdx = "numeric",
                                        paramsValues = "list",  # ValuesList object list
                                        order = "integer",  # vector of paramValues indices order
                                        orderIdx = "numeric", 
                                        finished = "logical"),  
                          
                          methods = list(
                            
                            # To be called after constructed. Sets initial values
                            new = function(pValues, execOrder){
                              paramsValues <<- list()
                              for (i in 1:length(pValues))
                                paramsValues[[i]] <<- cloneValueList(pValues[[i]]) # objects have to be cloned
                              finished <<- FALSE
                              print("Params execution order:")
                              print(execOrder)
                              order <<- execOrder
                              orderIdx <<- 1
                              paramIdx <<- order[orderIdx]
                            }
                          )
)


cloneValueList = function(vl){
  myVl <- ValuesList()
  myListOfValues <- list()
  for (i in 1:length(vl$listOfValues)){
    myVal <- Value()
    myVal$new(vl$listOfValues[[i]]$value)
    myListOfValues[[i]] <- myVal 
  }
  
  myVl$new(myListOfValues, vl$paramName, vl$fixed, vl$isNumeric)
  return(myVl)
}


# returns a sorted sublist with those values having best sensitivity. List is sorted by specificity
filterAndSortBestValues <- function(listOfValues, minSpecLimit, fixedIdx, optimizerParams){
  
  # get current default sensitivity levels
  defaultWSSensitivity <- listOfValues[[fixedIdx]]$ws_sensitivity
  defaultGeneSensitivity <- listOfValues[[fixedIdx]]$gene_sensitivity
  defaultROISensitivity <- listOfValues[[fixedIdx]]$roi_sensitivity
  
  # get best sensitivity
  best <- NULL
  bestList <- list()
  print("Filtering and sorting best values...")
  
  if (length(listOfValues) > 1){
    
    for (i in 1:length(listOfValues)){
      if (listOfValues[[i]]$isBetter(best)){
        limit <- minSpecLimit
        if (listOfValues[[i]]$ws_sensitivity > defaultWSSensitivity) {
          limit <- limit - limit * optimizerParams$allowedWSloss * 0.01 
        } else if (listOfValues[[i]]$gene_sensitivity > defaultGeneSensitivity){
          limit <- limit - limit * optimizerParams$allowedGSloss * 0.01
        } else if (listOfValues[[i]]$roi_sensitivity > defaultROISensitivity){
          limit <- limit - limit * optimizerParams$allowedRSloss * 0.01
        }
      
        if(listOfValues[[i]]$ws_specificity >= limit)
          best <- listOfValues[[i]]$clone()
      } 
    }
    
    # if no one is selected, select default one
    if (is.null(best))
      best <- listOfValues[[fixedIdx]]$clone()

    # build a list of best params sorted by specificity
    for (i in 1:length(listOfValues)) {
      if (listOfValues[[i]]$status != "error") {
        if (listOfValues[[i]]$ws_sensitivity == best$ws_sensitivity & listOfValues[[i]]$gene_sensitivity == best$gene_sensitivity & listOfValues[[i]]$roi_sensitivity == best$roi_sensitivity){
          if (length(bestList) == 0) {
            bestList[[1]] <- listOfValues[[i]]
          } else {
            idx <- 1
            while(bestList[[idx]]$ws_specificity > listOfValues[[i]]$ws_specificity) {
              idx <- idx + 1
              if (idx > length(bestList))
                break
            }
            
            bestList <- append(bestList, listOfValues[[i]], after=(idx - 1))
          }
        }
      }
    }
    
  } else{
    bestList[[1]] <- best
  }
    
  return(bestList)  
}


# returns list os params with fixed values [paramName, paramFixedValue], except specified value for current param idx
getParams = function(paramsValues, currentParamIdx, value){
  result <- list()
  for (i in 1:length(paramsValues)){
    pv <- paramsValues[[i]]
    
    if (i != currentParamIdx)
      result[[pv$paramName]] <- pv$fixed
    else
      result[[pv$paramName]] <- value
  }

  return(result)
}

# Returns execution sequence made of idxs
getIdxs = function(paramsValues, currentParamIdx, valueIdx){
  result = ""
  for (i in 1:length(paramsValues)){
    pv <- paramsValues[[i]]
    if (i != currentParamIdx)
      result <- paste0(result, "_", pv$fixedIdx)
    else
      result <- paste0(result, "_", valueIdx)
  }
  result <- substr(result, 2, nchar(result))
  
  return(result)
}

calculateMinSpec = function(params_vals, stats, datasetName, algFolder, datasetContent, optimizerParams, algOwnParams, bedFile){
  seq <- ""
  pars <- list()
  for (i in 1:length(params_vals)){
    pv <- params_vals[[i]]
    seq <- paste0(seq, "_" , pv$fixedIdx)
    pars[[pv$paramName]] <- pv$fixed
  }
  seq <- substr(seq, 2, nchar(seq))
  val <- params_vals[[1]]$listOfValues[[1]]$clone()
  val$launch(pars, seq, algFolder, datasetContent, optimizerParams, algOwnParams)

  while (val$status != "done") {
    val$checkFinished(stats, datasetName, algFolder, bedFile)
    Sys.sleep(5)
  }
  print("Initial execution:")
  print(val$getStr())

  return(val$ws_specificity)
}


# stats: SummaryStats object
# algFolder: algorithm folder where precalc folder exists and logs/output folders will be created
optimize = function(params_vals, stats, datasetName, datasetContent, optimizerParams, bedFile, algFolder, algOwnParams){
  
  # List of executing/executed params values. Used to avoid repeated. Key: sequence, value: Value
  execSequences <- list()

  # failed sequences
  failedSequences <- c()
  
  # list of executions with best sensitivity
  listOfBests <- list()
  
  # calculate default
  minSpecificity <- calculateMinSpec(params_vals, stats, datasetName, algFolder, datasetContent, optimizerParams, algOwnParams, bedFile)

  # the alg. will start from different params
  startPoints = list()

  for (i in 1:length(paramsValues)){
  # for (i in c(2, 7)){
      
    # define random order to execute params.
    execOrder <- sample(1:length(paramsValues))  
    while(execOrder[1] != i)
      execOrder <- sample(1:length(paramsValues))
    
    # Create start point
    sp <- StartPoint()
    sp$new(params_vals, execOrder)
    startPoints[[length(startPoints) + 1]] <- sp 
  }
  
  allFinished <- FALSE
  while(!allFinished){
    
    allFinished <- TRUE
    for (startPoint in startPoints){
      
      if (!startPoint$finished) {
        allFinished <- FALSE
        paramIdx <- startPoint$paramIdx
        paramsVals <- startPoint$paramsValues
        param <- paramsVals[[paramIdx]]
        print(paste("Processing param", param$paramName))
        
        paramFinished <- TRUE
        
        # try all values for this param
        for (i in 1:length(param$listOfValues)){
          val <- param$listOfValues[[i]]
          
          if (val$status == "executing") {
            val$checkFinished(stats, datasetName, algFolder, bedFile)
          } else if (val$status == "none") {
            seq <- getIdxs(paramsVals, paramIdx, i)
            if (seq %in% names(execSequences)) {
              val$alreadyDone(execSequences[[seq]]) # updates val status and stats if executed sequence is finished
            } else {
              if(val$launch(getParams(paramsVals, paramIdx, val$value), seq, algFolder, datasetContent, optimizerParams, algOwnParams)){
                execSequences[[val$valuesSequence]] <- val
              }
            }
          } 
          if (val$status == "none" | val$status == "executing") {
            paramFinished <- FALSE
          } else if (val$status == "error") {
            failedSequences <- c(failedSequences, val$valuesSequence)
          }
        }
        
        if(paramFinished){
          param$fixBestValue(minSpecificity, optimizerParams)
          cat("Param Finished\n\n\n\n")
          bestOfParam <- param$listOfValues[[param$fixedIdx]]
          cat(paste("bestOfParam", bestOfParam$getStr()))
          param$isDone <- TRUE
          
          # next param
          startPoint$orderIdx <- startPoint$orderIdx %% length(paramsVals) + 1
          startPoint$paramIdx <- startPoint$order[startPoint$orderIdx] 
          
          # if next param is done, loop finished
          if (paramsVals[[startPoint$paramIdx]]$isDone ) {
            startPoint$finished <- TRUE
            print("Start Point finished")

            if ( (length(listOfBests) == 0) || (bestOfParam$isBetter(listOfBests[[1]])) ) {  # improves current sensitivity
              listOfBests <- list(bestOfParam$clone())
            } else if (listOfBests[[1]]$ws_sensitivity == bestOfParam$ws_sensitivity & listOfBests[[1]]$gene_sensitivity == bestOfParam$gene_sensitivity & listOfBests[[1]]$roi_sensitivity == bestOfParam$roi_sensitivity){ # tie
              for (i in 1:length(listOfBests)) {
                if (listOfBests[[i]]$ws_specificity < bestOfParam$ws_specificity) {
                  listOfBests <- append(listOfBests, bestOfParam$clone(), after = (i - 1))
                } else if(i == length(listOfBests)) {  # end achieved
                  listOfBests <- append(listOfBests, bestOfParam$clone())
                }
              }
            }
            
            if (length(listOfBests) > 0)
                cat(paste("\nCurrently best params:", listOfBests[[1]]$getStr()))
          }
        }
      }
    }
    
    # avoid loop checking all time
    Sys.sleep(5)
  }
  
  print("FINISHING...")
  
  cat(paste("\n\nError sequences: "))
  cat(unique(failedSequences))
  
  if (length(listOfBests) > 0) {
    cat(paste("\n\nBest params:", listOfBests[[1]]$getStr()))
    listOfBests[[1]]$generateYaml(algFolder)
  }

  
  print(paste("Finishing at",Sys.time())) 
}