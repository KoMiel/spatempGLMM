
# import packages

library(dplyr)
library(rjson)
library(foreach)
library(doParallel)

registerDoParallel(cores=30)


# functions

# this function generates a list of matrices (one for each distribution/sampling event)
# with distances of all locations (presence or absence) to all presence locations
# this is then later used as a base in the calculation of the autocovariate

dist_list <- function(Y, coords, year){
  
  # generate an empty list for the distance matrices
  distances <- list()

  # counter to append matrices later on
  counter <- 1
  
  # loop over sampling events (called year here, but could be anything)
  for (k in (min(year)):(max(year))) {
    
    # filter all observations of one sampling event
    Y_Year <- Y[year == k] # the presence/absence encoding
    coordsYear <- coords[year == k, ] # the position
    presCoordsYear <- coordsYear[Y_Year == 1, ] # positions, but only those of presences
    
    # generate an empty matrix
    dist <- matrix(0, nrow = nrow(coordsYear), ncol = nrow(presCoordsYear))
    
    # loop over all positions (presences and absences)
    for (i in 1:nrow(coordsYear)){
      
      # loop over all presence positions
      for (j in 1:nrow(presCoordsYear)){
        
        # calculate all the distances from the i-th position to all presences
        dist[i,j] <- sqrt((coordsYear[i,1] - presCoordsYear[j,1])*(coordsYear[i,1] - presCoordsYear[j,1]) + (coordsYear[i,2] - presCoordsYear[j,2])*(coordsYear[i,2] - presCoordsYear[j,2]))
        
        #if distance is 0 (same position), set value to NA
        if (dist[i,j] == 0) {
          dist[i,j] <- NA
        }
      }
    }
    
    # add the matrix to list and continue with counter
    distances[[counter]] <- dist
    counter <- counter+1
  }
  
  # return the list
  return(distances)
}

# this function calculates the autocovariate from the distance matrices (output of dist_list)

autocovFunc <- function(distances, p, Y, year){
  
  # generate a vector
  autocov <- vector('numeric')
  
  # counter to select the distances matrix
  counter <- 1
  
  # loop over sampling events (called year here, but can be anything)
  for (k in (min(year)):(max(year))) {
    
    # get the data that belongs to the sampling event
    distYear <- distances[[counter]] # distances matrix
    Y_Year <- Y[year == k] # the presence/absence encoding
    
    # empty matrix for autocovariate
    autocovYear <- matrix(nrow = nrow(distYear), ncol = 1, 0)
    
    # loop over all presences
    for (presence in 1:ncol(distYear)){
      distTemp <- distYear
      
      # set the influence of that presence to 0
      distTemp[Y_Year == 0, presence] <- NA
      
      # calculate the autocovariate that is given from each of the presences
      weightedDistYear <- exp( - distTemp * p)

      # for each presence or absence, add up the autocovariate of all presences
      autocovYear <- autocovYear + rowSums(weightedDistYear, na.rm = TRUE)/ncol(distYear)
    }
    
    # append the result to the vector
    autocov <- append(autocov, autocovYear)
    
    # continue with counter
    counter <- counter + 1
  }
  
  # return the autocovariate
  return(autocov)
}


# read settings file
f <- "../settings.json"
settings <- fromJSON(file = f)

# read scenario file
f <- "../scenarios.json"
scenarios <- fromJSON(file = f)

# reading number of absences
absNum <- settings$absNum

# read the model list
modelList <- read.table("../modelList.csv", sep = ',')

# read directories
directoryDatasets <- settings$directory$datasets
directoryModels <- settings$directory$models

# generate a random seed and use it
seed <- runif(n = 1) * 1000000
set.seed(seed)

# save the random seed 
randomSeed <- c(seed, '\n')
filename <- paste(c(directoryDatasets, "prep_seed.txt"), collapse = "")
dir.create(path = directoryModels, showWarnings = FALSE) # generate directories if not already existing
sink(file = filename, append = TRUE)
cat(randomSeed)
sink()

# considered options for p_c
p_c <- c(0.64,0.8,1,1.25,1.5625)

# loop over all scenarios (parallel computing)
foreach (line = 1:nrow(modelList)) %dopar% {

  # get scenario and grid from list
  scenario <- modelList[line, 1]
  lscape <- modelList[line, 2]
  
  # get the subpath for scenario
  subpath <- scenarios$scenarios[[scenario]]$path
    
  # read the data
  filename <- paste(c(directoryDatasets, subpath, "dataset", lscape, ".txt"), collapse = "")
  dataFrame <- read.table(filename, head=F)
  
  # split in presences and absences
  presFrame <- dataFrame[dataFrame$V1 == 1,]
  absFrame <- dataFrame[dataFrame$V1 == 0,]
  
  # sample a number of absences per sampling event
  absFrame <- absFrame %>% group_by(V4) %>% sample_n(size = absNum)
  absFrame <- as.data.frame(absFrame)
  
  # combine absences and presences again
  dataFrame <- rbind(presFrame, absFrame)
  
  # order by sampling event
  dataFrame <- dataFrame[order(dataFrame$V4),]
  dataFrame <- dataFrame %>% mutate(id = row_number())

  # calculate distances object using dist_list function
  distances <- dist_list(Y = dataFrame[,"V1"], coords = dataFrame[,c("V2", "V3")], year = dataFrame[,"V4"])
  
  # calculate autocovariate for different p_c
  for (value in p_c) {
    dataFrame <- cbind(dataFrame, autocovFunc(distances = distances, p = value, Y = dataFrame[,"V1"], year = dataFrame[,"V4"]))
    names(dataFrame)[ncol(dataFrame)] <- paste0("p_c_", toString(value))
  }
  
  # save the resulting (sampled) dataframe
  dir.create(path = paste(c(directoryDatasets, subpath), collapse = ""), showWarnings = FALSE)
  filename <- paste(c(directoryDatasets, subpath, "sampled_dataset_", lscape, ".rdata"), collapse = "")
  save(dataFrame, file = filename)

}
