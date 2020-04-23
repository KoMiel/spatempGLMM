
# import packages

library(rjson)
library(pgirmess)
library(foreach)
library(doParallel)

registerDoParallel(cores=30)

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

# read directories
directoryResults <- settings$directory$results

# loop over all scenarios
moransI <- foreach (row = 1:nrow(modelList)) %dopar% {
  
  # print message (progress keeping)
  print(row)
  
  # generate empty structure to store results
  moran <- data.frame(matrix(nrow = 10, ncol = 10, 0))
  names(moran) <- c("baseCoeff", "baseP", "randCoeff", "randP", "autoCoeff", "autoP", "combCoeff", "combP", "compCoeff", "compP")

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
  moran[, c("baseCoeff", "baseP")] <- res[,2:3]/nrow(modelList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"randRes"], method = 'Moran', nbclass = distanceClasses)
  moran[, c("randCoeff", "randP")] <- res[,2:3]/nrow(modelList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"autoRes"], method = 'Moran', nbclass = distanceClasses)
  moran[, c("autoCoeff", "autoP")] <- res[,2:3]/nrow(modelList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"combRes"], method = 'Moran', nbclass = distanceClasses)
  moran[, c("combCoeff", "combP")] <- res[,2:3]/nrow(modelList)
        
  res <- correlog(coords = residuals[,c("V2", "V3")], z = residuals[,"compRes"], method = 'Moran', nbclass = distanceClasses)
  moran[, c("compCoeff", "compP")] <- res[,2:3]/nrow(modelList)
    
  moran
}
    
# generate empty structure to store results
moransICombined <- data.frame(matrix(nrow = 10, ncol = 10, 0))
names(moransICombined) <- c("baseCoeff", "baseP", "randCoeff", "randP", "autoCoeff", "autoP", "combCoeff", "combP", "compCoeff", "compP")

# combine results
for (i in 1:length(moransI)){
  moransICombined <- moransICombined + moransI[[i]]
}

# save the results
filename <- paste(c(directoryResults, "moransI.rdata"), collapse = "")
save(moransICombined, file = filename)
filename <- paste(c(directoryResults, "moransI.csv"), collapse = "")
write.csv(moransICombined, file = filename)
