

### import packages ###

library(dplyr)
library(rjson)
library(spaMM)
library(nlme)
library(stats)
library(pROC)


### main program ###

# read settings file
f <- "../../settings.json"
settings <- fromJSON(file = f)

# get directories
directoryModels <- settings$directory$models
directoryCasestudy <- settings$directory$casestudy

# the options for the year
folds <- c(2001,2002,2003,2004,2005,2006,2007)

# load the data
filename <- paste(c(directoryCasestudy, "corncrakes/survey2Frame.csv"), collapse = "")
dataFrame <- read.csv(filename)

# generate empty dataframe for information on test observations
testFrames <- data.frame()

# generate empty lists for predictions of the two models
probsTemp <- list()
probsSpa <- list()

# loop over all folds
for (fold in folds) {
    
  # get the test dataset
  testIndexes <- which(dataFrame$year==fold,arr.ind=TRUE)
  testFrame <- dataFrame[testIndexes, ]

  # load the baseline GLMM
  filename <- paste(c(directoryModels, "casestudy/temporalGLMM_", fold, ".rdata"), collapse = "")
  load(file = filename)
  
  # calculate the predictions of the model
  probTemp <- as.vector(predict(temporalGLMM, newdata = testFrame))
  
  # append to the list of predictions
  probsTemp <- c(probsTemp, probTemp)
  
  # same for the true values
  testFrames <- rbind(testFrames, testFrame)

  # load the spatiotemporal GLMM
  filename <- paste(c(directoryModels, "casestudy/spatiotemporalGLMM_", fold, ".rdata"), collapse = "")
  load(file = filename)

  # calculate the predictions of the model
  probSpa <- as.vector(predict(spatiotemporalGLMM, newdata = testFrame))

  # append to the list of predictions
  probsSpa <- c(probsSpa, probSpa)

}

# calculate the roc for both models
rocTemp <- roc(Y ~ unlist(probsTemp), data = testFrames)
rocSpa <- roc(Y ~ unlist(probsSpa), data = testFrames)

# get the auc
aucTemp <- rocTemp$auc
aucSpa <- rocSpa$auc

# get the auc variance
errTemp <- sqrt(var(rocTemp))
errSpa <- sqrt(var(rocSpa))

# save the aucs and variances
aucVector <- c(aucTemp, errTemp, aucSpa, errSpa)
filename <- paste(c(directoryModels, "casestudy/auc.txt"), collapse = "")
sink(file = filename, append = TRUE)
cat(aucVector)
sink()
