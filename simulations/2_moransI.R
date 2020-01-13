

### import packages ###

library("rjson")
library("pgirmess")
library("spaMM")


### import settings and scenarios ###

# command line input to select the scenario/grid combination to work on

line <- commandArgs(trailingOnly = TRUE)
line <- strtoi(line)

# read settings file

f <- "../../settings.json"
settings <- fromJSON(file = f)

# read scenario file

f <- "../../scenarios.json"
scenarios <- fromJSON(file = f)

# read the number of distance classes for Moran's I calculation

distanceClasses <- settings$distanceClasses

# read the model list

modelList <- read.table("../../modelList.csv", sep = ',')

# find the scenario and grid from the lit

scenario <- modelList[line, 1]
grid <- modelList[line, 2]

# read directories

directoryDatasets <- settings$directory$datasets
directoryModels <- settings$directory$models
directoryPlots <- settings$directory$plots
directoryResults <- settings$directory$results
subpath <- scenarios$scenarios[[scenario]]$path

# read autocovariate information

autocov <- settings$autocovariate


### main program ###

# load the fitted models
filename <- paste(c(directoryModels, subpath, "temporalGLMMmodel_", grid, ".rdata"), collapse = "")
load(file = filename)
  
filename <- paste(c(directoryModels, subpath, "spatiotemporalGLMMmodel_", grid, ".rdata"), collapse = "")
load(file = filename)
    
# get the data the models were fitted with
dataSpatiotemporalGLMM <- spatiotemporalGLMM$data
dataTemporalGLMM <- temporalGLMM$data

# get the residuals of the models
dataSpatiotemporalGLMM$residuals <- residuals(spatiotemporalGLMM)
dataTemporalGLMM$residuals <- residuals(temporalGLMM)
      
# throw everything away besides the coordinates and the residuals
dataSpatiotemporalGLMM <- data.frame(dataSpatiotemporalGLMM$V2, dataSpatiotemporalGLMM$V3, dataSpatiotemporalGLMM$residuals)
dataTemporalGLMM <- data.frame(dataTemporalGLMM$V2, dataTemporalGLMM$V3, dataTemporalGLMM$residuals)
      
# calculate Moran's I
corTemp <- correlog(coords = dataTemporalGLMM[,1:2], z = dataTemporalGLMM[,3], method = 'Moran', nbclass = distanceClasses)
corSpa <- correlog(coords = dataSpatiotemporalGLMM[,1:2], z = dataSpatiotemporalGLMM[,3], method = 'Moran', nbclass = distanceClasses)

# save the calculated Moran's I (first the grid number to indicate the grid worked on, then the Moran's I values and then the p values)
spatiotemporalGLMM <- c(grid, corSpa[,2], corSpa[,3], '\n')   
filename <- paste(c(directoryModels, subpath, "moranCorrelogSpatiotemporal.txt"), collapse = "")
sink(file = filename, append = TRUE)
cat(spatiotemporalGLMM)
sink()

# save the object
filename <- paste(c(directoryModels, subpath, "moranCorrelogSpatiotemporal_", grid, ".rdata"), collapse = "")
save(corSpa, file = filename)
    
# save the calculated Moran's I (first the grid number to indicate the grid worked on, then the Moran's I values and then the p values)
temporalGLMM <- c(grid, corTemp[,2], corTemp[,3], '\n')   
filename <- paste(c(directoryModels, subpath, "moranCorrelogTemporal.txt"), collapse = "")
sink(file = filename, append = TRUE)
cat(temporalGLMM)
sink()

# save the object
filename <- paste(c(directoryModels, subpath, "moranCorrelogTemporal_", grid, ".rdata"), collapse = "")
save(corTemp, file = filename)