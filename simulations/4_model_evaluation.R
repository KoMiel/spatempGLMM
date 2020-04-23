
# import packages

library(rjson)
library(pROC)

# read settings file
f <- "../settings.json"
settings <- fromJSON(file = f)

# read scenario file
f <- "../scenarios.json"
scenarios <- fromJSON(file = f)

# read directories
directoryDatasets <- settings$directory$datasets
directoryModels <- settings$directory$models
directoryResults <- settings$directory$results

# read the model list
modelList <- read.table("../modelList.csv", sep = ',')

# get lists of all scenarios and landscapes
scenarioList <- unique(modelList[,1])
lscapeList <- unique(modelList[,2])

# considered options for p_c
p_cs <- c(0.64,0.8,1,1.25,1.5625)

# loop over scenarios

for (scenario in scenarioList){

  # read subpath
  subpath <- scenarios$scenarios[[scenario]]$path
  
  # create directories if not already existing
  dir.create(path = directoryResults, showWarnings = FALSE)
  dir.create(path = paste(c(directoryResults, subpath), collapse = ""), showWarnings = FALSE)
  
  # empty data frame for AUCs
  AUCs <- data.frame(matrix(nrow = length(lscapeList), ncol = 6, 0))
  names(AUCs) <- c("base", "auto", "rand", "comb", "comp", "yearly")
  
  # empty data frame for model slopes
  parameters <- data.frame(matrix(nrow = length(lscapeList), ncol = 34))
  names(parameters) <- c("base_mean(b2)", "base_sd(b2)", "base_mean(b3)", "base_sd(b3)",
                         "rand_mean(b2)", "rand_sd(b2)", "rand_mean(b3)", "rand_sd(b3)",
                         "auto_mean(b2)", "auto_sd(b2)", "auto_mean(b3)", "auto_sd(b3)", "auto_mean(bc)", "auto_sd(bc)",
                         "comb_mean(b2)", "comb_sd(b2)", "comb_mean(b3)", "comb_sd(b3)", "comb_mean(bc)", "comb_sd(bc)",
                         "comp_mean(b2)", "comp_sd(b2)", "comp_mean(b3)", "comp_sd(b3)", "comp_mean(bc)", "comp_sd(bc)", "comp_mean(b1)", "comp_sd(b1)",
                         "yearly_mean(b2)", "yearly_sd(b2)", "yearly_mean(b3)", "yearly_sd(b3)", "yearly_mean(bc)", "yearly_sd(bc)")
      
  # loop over landcapes
  for (lscape in lscapeList) {
  
    # load models with autocovariate component
    filename <- paste(c(directoryModels, subpath, lscape, "/resultsAutocov.rdata"), collapse = "")
    load(filename)
    
    # empty data frame for temporary results (used to select the best model out of the 5 candidates)
    AUC_temp <- data.frame(matrix(nrow = length(resultsAutocov), ncol = 4, 0))
    names(AUC_temp) <- c("auto", "comb", "comp", "sd")
    
    # loop over all options for the range parameter
    for (p_c in 1:length(p_cs)) {
    
      # get the predictions on unknown data
      testResults <- resultsAutocov[[p_c]][[2]]
      
      # get years/sampling events in data
      years <- unique(testResults$V4)
      
      # loop over all years
      for (year in years) {
      
        # get results of one year only
        testResultsYear <- testResults[testResults$V4 == year,]
        
        # calculate and sum up AUCs of individual years
        AUC_temp[p_c, "auto"] <- AUC_temp[p_c, "auto"] + auc(V1 ~ predAuto, data = testResultsYear) * sum(testResultsYear$V1)/sum(testResults$V1)
        AUC_temp[p_c, "comb"] <- AUC_temp[p_c, "comb"] + auc(V1 ~ predComb, data = testResultsYear) * sum(testResultsYear$V1)/sum(testResults$V1)
        AUC_temp[p_c, "comp"] <- AUC_temp[p_c, "comp"] + auc(V1 ~ predComp, data = testResultsYear) * sum(testResultsYear$V1)/sum(testResults$V1)
        
        # calculate the standard deviations (for standardization of the autocovariate beta)
        AUC_temp[p_c, "sd"] <- AUC_temp[p_c, "sd"] + sd(testResultsYear$Autocovariate)/10
      }
    }
    
    # get the best AUC and save it
    AUCs[lscape + 1, "auto"] <- max(AUC_temp[, "auto"])
    AUCs[lscape + 1, "comb"] <- max(AUC_temp[, "comb"])
    AUCs[lscape + 1, "comp"] <- max(AUC_temp[, "comp"])
    
    # get the set of estimated slopes of the best models, standardize mean and sd of autocovariate beta and store in data frame
    bestAuto <- which.max(AUC_temp[, "auto"])
    pars <- resultsAutocov[[bestAuto]][[3]]
    parsAverage <- apply(pars, 2, FUN = function(x) mean(unlist(x)))
    parsAverage[5:6] <- parsAverage[5:6] * AUC_temp[bestAuto, "sd"]/AUC_temp[3, "sd"]
    parameters[lscape + 1, 9:14] <- parsAverage[1:6]
    
    bestComb <- which.max(AUC_temp[, "comb"])
    pars <- resultsAutocov[[bestComb]][[3]]
    parsAverage <- apply(pars, 2, FUN = function(x) mean(unlist(x)))
    parsAverage[11:12] <- parsAverage[11:12] * AUC_temp[bestComb, "sd"]/AUC_temp[3, "sd"]
    parameters[lscape + 1, 15:20] <- parsAverage[7:12]
    
    bestComp <- which.max(AUC_temp[, "comp"])
    pars <- resultsAutocov[[bestComp]][[3]]
    parsAverage <- apply(pars, 2, FUN = function(x) mean(unlist(x)))
    parsAverage[17:18] <- parsAverage[17:18] * AUC_temp[bestComp, "sd"]/AUC_temp[3, "sd"]
    parameters[lscape + 1, 21:28] <- parsAverage[13:20]
    
    # get the residuals of the best models (on training data), average over sampling events
    trainResults <- resultsAutocov[[bestAuto]][[1]]
    trainResults$autoRes <- trainResults$V1 - trainResults$predAuto
    autoRes <- aggregate(autoRes ~ V2 + V3, data = trainResults, FUN = mean)
    
    trainResults <- resultsAutocov[[bestComb]][[1]]
    trainResults$combRes <- trainResults$V1 - trainResults$predComb
    combRes <- aggregate(combRes ~ V2 + V3, data = trainResults, FUN = mean)
    
    trainResults <- resultsAutocov[[bestComp]][[1]]
    trainResults$compRes <- trainResults$V1 - trainResults$predComp
    compRes <- aggregate(compRes ~ V2 + V3, data = trainResults, FUN = mean)
    
    # load models with no autocovariate component
    filename <- paste(c(directoryModels, subpath, lscape, "/resultsNoautocov.rdata"), collapse = "")
    load(filename)
    
    # get the predictions on unknown data
    testResults <- resultsNoautocov[[2]]
    
    # get years/sampling events in data
    years <- unique(testResults$V4)
      
    # loop over years
    for (year in years) {
    
      # get results of one year only
      testResultsYear <- testResults[testResults$V4 == year,]
        
      # calculate and sum up AUCs of individual years
      AUCs[lscape + 1, "base"] <- AUCs[lscape + 1, "base"] + auc(V1 ~ predBase, data = testResultsYear) * sum(testResultsYear$V1)/sum(testResults$V1)
      AUCs[lscape + 1, "rand"] <- AUCs[lscape + 1, "rand"] + auc(V1 ~ predRand, data = testResultsYear) * sum(testResultsYear$V1)/sum(testResults$V1)
    }
    
    # get the set of estimated slopes of the models and store them in data frame
    pars <- resultsNoautocov[[3]]
    parsAverage <- apply(pars, 2, FUN = function(x) mean(unlist(x)))
    parameters[lscape + 1, 1:8] <- parsAverage
    
    # get the residuals on training data
    trainResults <- resultsNoautocov[[1]]
    trainResults$baseRes <- trainResults$V1 - trainResults$predBase
    trainResults$randRes <- trainResults$V1 - trainResults$predRand
    residuals <- aggregate(cbind(baseRes, randRes) ~ V2 + V3, data = trainResults, FUN = mean)
    
    # merge all residuals together, average over sampling events
    residuals <- merge(residuals, autoRes, by = c("V2", "V3"))
    residuals <- merge(residuals, combRes, by = c("V2", "V3"))
    residuals <- merge(residuals, compRes, by = c("V2", "V3"))
    
    # create directory for residuals and store them
    dir.create(path = paste(c(directoryResults, subpath, lscape), collapse = ""), showWarnings = FALSE)
    filename <- paste(c(directoryResults, subpath, lscape, "/residuals.rdata"), collapse = "")
    save(residuals, file = filename)
    
    # load models for one year evaluation
    filename <- paste(c(directoryModels, subpath, lscape, "/resultsOneyearAutocov.rdata"), collapse = "")
    load(filename)
    
    # get the predictions on testing data
    testResults <- resultsAutocov[[2]]
    
    # get years/sampling events in data
    years <- unique(testResults$V4)
    
    # loop over years
    for (year in years) {
      
      # get results of one year only
      testResultsYear <- testResults[testResults$V4 == year,]
      
      # calculate and sum up AUCs of individual years
      AUCs[lscape + 1, "yearly"] <- AUCs[lscape + 1, "yearly"] + auc(V1 ~ predComb, data = testResultsYear) * sum(testResultsYear$V1)/sum(testResults$V1)
    }
    
    # calculate averages of parameters and store them
    pars <- resultsAutocov[[3]]
    parsAverage <- apply(pars, 2, FUN = function(x) mean(unlist(x)))
    parameters[lscape + 1, 29:34] <- parsAverage
  }
  
  # store the AUCs
  filename <- paste(c(directoryResults, subpath, "AUC.rdata"), collapse = "")
  save(AUCs, file = filename)
  filename <- paste(c(directoryResults, subpath, "AUC.txt"), collapse = "")
  write.table(AUCs, filename, row.names=FALSE)
  
  # store the estimated slopes
  filename <- paste(c(directoryResults, subpath, "slopes.rdata"), collapse = "")
  save(parameters, file = filename)
  filename <- paste(c(directoryResults, subpath, "slopes.txt"), collapse = "")
  write.table(parameters, filename, row.names=FALSE)
}
