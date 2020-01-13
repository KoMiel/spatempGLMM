

### import packages ###

library("dplyr")
library("rjson")
library("spaMM")


### functions ###

# this function generates a list of matrices (one for each distribution/sampling event)
# with distances of all locations (presence or absence) to all presence locations
# this is then later used as a base in the calculation of the autocovariate

dist_list <- function(y, coords, year){
  
  # generate an empty list for the distance matrices
  distances <- list()
  # loop over sampling events (called year here, but can be anything)
  for (k in 0:(max(year))) {
    # filter all observations of one sampling event
    yYear <- y[year == k] #the presence/absence encoding
    coordsYear <- coords[year == k, ] #the position
    presCoordsYear <- coordsYear[yYear == 1, ] #positions, but only presences
    # generate an empty matrix
    dist <- matrix(0, nrow = nrow(coordsYear), ncol = nrow(presCoordsYear))
    # loop over all positions (presences and absences)
    for (i in 1:nrow(coordsYear)){
      # loop over all presence positions, calculate all the distances from the i-th position to all presences
      for (j in 1:nrow(presCoordsYear)){
        dist[i,j] <- sqrt((coordsYear[i,1] - presCoordsYear[j,1])*(coordsYear[i,1] - presCoordsYear[j,1]) + (coordsYear[i,2] - presCoordsYear[j,2])*(coordsYear[i,2] - presCoordsYear[j,2]))
        #if distance is 0 (same position), set value to NA
        if (dist[i,j] == 0) {
          dist[i,j] <- NA
        }
      }
    }
    distances[[k+1]] <- dist #add the matrix to the list
  }
  # return the list
  return(distances)
}


# this function calculates the autocovariate from the distanc matrices (output of dist_list)

autocov <- function(distances, year, param){
  # generate a vector
  autocov <- vector('numeric')
  # loop over sampling events (called year here, but can be anything)
  for (k in 0:(max(year))) {
    # get the matrix that belongs to the sampling event
    distYear <- distances[[k+1]]
    # calculate the autocovariate that is given from each of the presences
    weightedDistYear <- exp(- distYear * param)
    # for each presence or absence, add up the conspecific interaction of all presences
    autocovYear <- rowSums(weightedDistYear, na.rm = TRUE)
    # append the result to the vector
    autocov <- append(autocov, autocovYear)
  }
  # return the vector
  return(autocov)
}


### import settings ### 

# command line input to select the scenario/grid combination to work on

line <- commandArgs(trailingOnly = TRUE)
line <- strtoi(line)

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

# find the scenario and grid from the lit

scenario <- modelList[line, 1]
grid <- modelList[line, 2]

# read directories

directoryDatasets <- settings$directory$datasets
directoryModels <- settings$directory$models
subpath <- scenarios$scenarios[[scenario]]$path

# read autocovariate information

autocov <- settings$autocovariate

# generate a random seed and use it

seed <- runif(n = 1, min = 0, max = 1000000)
set.seed(seed)

# save seed 
randomSeed <- c(grid, seed)
filename <- paste(c(directoryModels, subpath, "seed.txt"), collapse = "")
dir.create(path = directoryModels, showWarnings = FALSE) # generate directories if not already existing
dir.create(path = paste(c(directoryModels, subpath), collapse = ""), showWarnings = FALSE)
capture.output(randomSeed, file = filename, append = TRUE)


### main program ###

# read the data
filename <- paste(c(directoryDatasets, subpath, "dataset", grid, ".txt"), collapse = "")
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

# save the resulting (sampled) dataframe
filename <- paste(c(directoryDatasets, subpath, "sampled_dataset_", grid, ".rdata"), collapse = "")
save(dataFrame, file = filename)

# calculate the standard deviations of the environmental variables (for re-calculation after the model fitting)
norm <- c(grid, sd(dataFrame$V6), sd(dataFrame$V7))

# standardize the environmental variables
for (var in names(dataFrame)) {
  if (var != "V1" && var != "V2" && var != "V3" && var != "V4") {
    dataFrame[var] <- (dataFrame[,var] - mean(dataFrame[,var]))/sd(dataFrame[,var])
  }
}

# calculate distances object using dist_list function
distances <- dist_list(y = dataFrame[,1], coords = dataFrame[,2:3], year = dataFrame[,4])

# calculate the autocovariate
dataFrame[,autocov$label] <- autocov(distance = distances, year = dataFrame[,4], param = autocov$parameter)

# store the standard deviation of the autocovariate
norm <- c(norm, sd(dataFrame[,autocov$label]))
filename <- paste(c(directoryModels, subpath,"sd.txt"), collapse = "")
capture.output(norm, file = filename, append = TRUE)
  
# standardize the autocovariate
dataFrame[,autocov$label] <- (dataFrame[,autocov$label] - mean(dataFrame[,autocov$label]))/sd(dataFrame[,autocov$label])
  
# fit the baseline model
temporalGLMM <- fitme(V1 ~ V6 + V7 + (1|V4) + dataFrame[,autocov$label], data = dataFrame, family = binomial)

# save the resulting model as well as its parameters (parameters are all saved to one file)
estimationTemporalGLMM <- c(grid, temporalGLMM$fixef[[2]], temporalGLMM$fixef[[3]], temporalGLMM$fixef[[4]], sd(dataFrame[,autocov$label]), sqrt(temporalGLMM$beta_cov[[2,2]]), sqrt(temporalGLMM$beta_cov[[3,3]]), sqrt(temporalGLMM$beta_cov[[4,4]]))
filename <- paste(c(directoryModels, subpath, "temporalGLMM.txt"), collapse = "")
capture.output(estimationTemporalGLMM, file = filename, append = TRUE)
filename <- paste(c(directoryModels, subpath, "temporalGLMMmodel_", grid, ".rdata"), collapse = "")
save(temporalGLMM, file = filename)

# fit the spatiotemporal model
spatiotemporalGLMM <- fitme(V1 ~ V6 + V7 + (1|V4) + dataFrame[,autocov$label] + Matern(1|V2+V3), data = dataFrame, family = binomial, fixed = list(nu = 0.5))
  
# save the resulting model as well as its parameters (parameters are all saved to one file, again)
estimationSpatiotemporalGLMM <- c(grid, spatiotemporalGLMM$fixef[[2]], spatiotemporalGLMM$fixef[[3]], spatiotemporalGLMM$fixef[[4]], sd(dataFrame[,autocov$label]), sqrt(spatiotemporalGLMM$beta_cov[[2,2]]), sqrt(spatiotemporalGLMM$beta_cov[[3,3]]), sqrt(spatiotemporalGLMM$beta_cov[[4,4]]))   
filename <- paste(c(directoryModels, subpath, "spatiotemporalGLMM.txt"), collapse = "")
capture.output(estimationSpatiotemporalGLMM, file = filename, append = TRUE)
filename <- paste(c(directoryModels, subpath, "spatiotemporalGLMMmodel_", grid, ".rdata"), collapse = "")
save(spatiotemporalGLMM, file = filename)