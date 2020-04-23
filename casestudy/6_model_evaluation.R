
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
directoryModels <- settings$directory$models
directoryResults <- settings$directory$results

# values of parameter p_c
p_cs <- c(0.1,0.1*1.584893,0.1*1.584893^2,0.1*1.584893^3,0.1*1.584893^4,
          1,1*1.584893,1*1.584893^2,1*1.584893^3,1*1.584893^4,
          10)

# create directories if not already existing
dir.create(path = directoryResults, showWarnings = FALSE)
dir.create(path = paste(c(directoryResults, 'casestudy'), collapse = ""), showWarnings = FALSE)

# empty data frame for AUCs
AUCs <- data.frame(matrix(nrow = length(p_cs), ncol = 8, 0))
names(AUCs) <- c("base", "auto", "rand", "comb", "sd_base", "sd_auto", "sd_rand", "sd_comb")
  
# empty data frame for model slopes
parameters <- data.frame(matrix(nrow = length(p_cs), ncol = 20))
names(parameters) <- c("base_mean(b2)", "base_sd(b2)", "base_mean(b3)", "base_sd(b3)",
                       "rand_mean(b2)", "rand_sd(b2)", "rand_mean(b3)", "rand_sd(b3)",
                       "auto_mean(b2)", "auto_sd(b2)", "auto_mean(b3)", "auto_sd(b3)", "auto_mean(bc)", "auto_sd(bc)",
                       "comb_mean(b2)", "comb_sd(b2)", "comb_mean(b3)", "comb_sd(b3)", "comb_mean(bc)", "comb_sd(bc)")

# load models with autocovariate component
filename <- paste(c(directoryModels, 'casestudy/resultsAutocov.rdata'), collapse = "")
load(filename)
    
counter <- 1
for (p_c in 1:length(p_cs)) {
    
  # get the predictions on unknown data
  testResults <- resultsAutocov[[p_c]][[2]]
      
  # get years/sampling events in data
  years <- unique(testResults$year)
      
  # loop over years
  for (year in years) {
    
    # get results of one year only
    testResultsYear <- testResults[testResults$year == year,]
        
    # calculate and sum up AUCs of individual years
    AUCs[counter, "auto"] <- AUCs[counter, "auto"] + auc(Y ~ predAuto, data = testResultsYear) * sum(testResultsYear$Y)/sum(testResults$Y)
    AUCs[counter, "comb"] <- AUCs[counter, "comb"] + auc(Y ~ predComb, data = testResultsYear) * sum(testResultsYear$Y)/sum(testResults$Y)

  }
      
  # loop over years again
  for (year in years) {
        
    # get results of one year only (again)
    testResultsYear <- testResults[testResults$year == year,]
        
    # calculate standard deviation for the AUCs
    AUCs[counter, "sd_auto"] <- AUCs[counter, "sd_auto"] + ((auc(Y ~ predAuto, data = testResultsYear) - AUCs[counter, "auto"]) * sum(testResultsYear$Y)/sum(testResults$Y))^2
    AUCs[counter, "sd_comb"] <- AUCs[counter, "sd_comb"] + ((auc(Y ~ predComb, data = testResultsYear) - AUCs[counter, "comb"]) * sum(testResultsYear$Y)/sum(testResults$Y))^2
  }
      
  # get square root of AUCs
  AUCs[counter, "sd_auto"] <- sqrt(AUCs[counter, "sd_auto"])
  AUCs[counter, "sd_comb"] <- sqrt(AUCs[counter, "sd_comb"])

  # get the parameters
  pars <- resultsAutocov[[p_c]][[3]]
  parsAverage <- apply(pars, 2, FUN = function(x) mean(unlist(x)))
  
  # store the parameters
  parameters[counter, 9:20] <- parsAverage
      
  # counter to next
  counter <- counter + 1
}

# load models with no autocovariate component
filename <- paste(c(directoryModels, 'casestudy/resultsNoautocov.rdata'), collapse = "")
load(filename)
    
# get the predictions on unknown data
testResults <- resultsNoautocov[[2]]
    
# get years/sampling events in data
years <- unique(testResults$year)
      
# loop over years
for (year in years) {
    
  # get results of one year only
  testResultsYear <- testResults[testResults$year == year,]
        
  # calculate and sum up AUCs of individual years
  AUCs[1,"base"] <- AUCs[1,"base"] + auc(Y ~ predBase, data = testResultsYear) * sum(testResultsYear$Y)/sum(testResults$Y)
  AUCs[1,"rand"] <- AUCs[1,"rand"] + auc(Y ~ predRand, data = testResultsYear) * sum(testResultsYear$Y)/sum(testResults$Y)
}

for (year in years) {
  
  # get results of one year only
  testResultsYear <- testResults[testResults$year == year,]
  
  # calculate standard deviation
  AUCs[1, "sd_base"] <- AUCs[1, "sd_base"] + ((auc(Y ~ predBase, data = testResultsYear) - AUCs[1, "base"]) * sum(testResultsYear$Y)/sum(testResults$Y))^2
  AUCs[1, "sd_rand"] <- AUCs[1, "sd_rand"] + ((auc(Y ~ predRand, data = testResultsYear) - AUCs[1, "rand"]) * sum(testResultsYear$Y)/sum(testResults$Y))^2
}

# take square root of sum
AUCs[1, "sd_base"] <- sqrt(AUCs[1, "sd_base"])
AUCs[1, "sd_rand"] <- sqrt(AUCs[1, "sd_rand"])

# get the slope estimates
pars <- resultsNoautocov[[3]]
parsAverage <- apply(pars, 2, FUN = function(x) mean(unlist(x)))

# store the slope estimates
parameters[1, 1:8] <- rep(parsAverage, 3)

# save the AUCs
filename <- paste(c(directoryResults, "casestudy/AUC.rdata"), collapse = "")
save(AUCs, file = filename)
filename <- paste(c(directoryResults, "casestudy/AUC.txt"), collapse = "")
write.table(AUCs, filename, row.names = FALSE)
  
# save the estimated slopes
filename <- paste(c(directoryResults, "casestudy/slopes.rdata"), collapse = "")
save(parameters, file = filename)
filename <- paste(c(directoryResults, "casestudy/slopes.txt"), collapse = "")
write.table(parameters, filename, row.names = FALSE)