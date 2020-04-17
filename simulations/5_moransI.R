
# import packages

library(rjson)
library(pgirmess)

# read settings file
f <- "../settings.json"
settings <- fromJSON(file = f)

# read scenario file
f <- "../scenarios.json"
scenarios <- fromJSON(file = f)

# read the number of distance classes for Moran's I calculation
distanceClasses <- settings$distanceClasses

# read the model list
modelList <- read.table("../modelList.csv", sep = ',')

# get lists of all scenarios and landscapes
scenarioList <- unique(modelList[,1])
lscapeList <- unique(modelList[,2])

# read directories
directoryResults <- settings$directory$results

# generate empty structure to store results
moransI <- data.frame(matrix(nrow = 10, ncol = 10, 0))
names(moransI) <- c("baseCoeff", "baseP", "randCoeff", "randP", "autoCoeff", "autoP", "combCoeff", "combP", "compCoeff", "compP")

# loop over all scenarios
for (row in 1:nrow(modelList)) {
  
  # print message (progress keeping)
  print(row)
  
  # get scenario and grid from list
  scenario <- modelList[row, 1]
  lscape <- modelList[row, 2]
  
  # get the subpath for scenario
  subpath <- scenarios$scenarios[[scenario]]$path
  
  # load the residuals
  filename <- paste(c(directoryResults, subpath, lscape, "/residuals.rdata"), collapse = "")
  load(file = filename)
  
  # calculate correlograms for all models and store them
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"baseRes"], method = 'Moran', nbclass = distanceClasses)
  moransI[, c("baseCoeff", "baseP")] <- moransI[, c("baseCoeff", "baseP")] + res[,2:3]/length(lscapeList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"randRes"], method = 'Moran', nbclass = distanceClasses)
  moransI[, c("randCoeff", "randP")] <- moransI[, c("randCoeff", "randP")] + res[,2:3]/length(lscapeList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"autoRes"], method = 'Moran', nbclass = distanceClasses)
  moransI[, c("autoCoeff", "autoP")] <- moransI[, c("autoCoeff", "autoP")] + res[,2:3]/length(lscapeList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"combRes"], method = 'Moran', nbclass = distanceClasses)
  moransI[, c("combCoeff", "combP")] <- moransI[, c("combCoeff", "combP")] + res[,2:3]/length(lscapeList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"compRes"], method = 'Moran', nbclass = distanceClasses)
  moransI[, c("compCoeff", "compP")] <- moransI[, c("compCoeff", "compP")] + res[,2:3]/length(lscapeList)
}
    
# save the results
filename <- paste(c(directoryResults, subpath, "moransI.rdata"), collapse = "")
save(moransI, file = filename)
