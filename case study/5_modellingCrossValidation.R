

### import packages ###

library(dplyr)
library(rjson)
library(spaMM)
library(nlme)
library(stats)


### function definitions ###

# a function that is used to calculate the value of lambda (see the supporting material of the paper)
lambda_func <- function(lambda, dist, n, independency){
  dist <- exp(- lambda * dist) #apply exponential function to distances
  value <- colSums(dist) * n #sum the influences of all positions
  value <- abs(mean(value) - (independency + 1) * mean(n)) #calculate the function that is to be optimized
  return(value)
}


### main program ###

# command line input to select the fold of the cross-validation
line <- commandArgs(trailingOnly = TRUE)
line <- strtoi(line)

# read settings file
f <- "../../settings.json"
settings <- fromJSON(file = f)

# get directories
directoryModels <- settings$directory$models
directoryCasestudy <- settings$directory$casestudy

# get the year of the cross-validation fold
folds <- c(2001,2002,2003,2004,2005,2006,2007)
fold <- folds[line]

# load the corncrake data
filename <- paste(c(directoryCasestudy, "corncrakes/survey2Frame.csv"), collapse = "")
dataFrame <- read.csv(filename)

# generate a random seed and save it
seed <- runif(n = 1, min = 0, max = 1000000)
set.seed(seed)
seed_vector <- c(fold, seed)
filename <- paste(c(directoryModels, "casestudy/seed_cv.txt"), collapse = "")
dir.create(path = directoryModels, showWarnings = FALSE)
dir.create(path = paste(c(directoryModels, "casestudy/"), collapse = ""), showWarnings = FALSE)
write.table(seed_vector, file = filename, row.names = FALSE, col.names = FALSE, append = TRUE)

# split the corncrake data in testing and training dataset
testIndexes <- which(dataFrame$year==fold,arr.ind=TRUE)
testFrame <- dataFrame[testIndexes, ]
trainFrame <- dataFrame[-testIndexes, ]

# data frame with only the presence observations for the lambda calculation
presenceFrame <- testFrame[testFrame$Y == 1,]

# generate an array with all unique presence positions
uniqueFrame <- data.frame(x = numeric(), y = numeric(), num = numeric())

for (i in unique(presenceFrame$x)) {
  insideFrame <- presenceFrame[presenceFrame$x == i,]
  for (j in unique(insideFrame$y)) {
    insiderFrame <- insideFrame[insideFrame$y == j,]
    newFrame <- data.frame(x = i, y = j, num = nrow(insiderFrame))
    uniqueFrame <- rbind(uniqueFrame, newFrame)
  }
}

# generate a matrix with distances between all unique presence positions
distancesMatrix <- matrix(0, nrow = nrow(uniqueFrame), ncol = nrow(uniqueFrame))

for (i in 1:nrow(uniqueFrame)) {
  distances <- sqrt((uniqueFrame[i,1] - uniqueFrame[,1])*(uniqueFrame[i,1] - uniqueFrame[,1]) + (uniqueFrame[i,2] - uniqueFrame[,2])*(uniqueFrame[i,2] - uniqueFrame[,2]))
  distancesMatrix[i,] <- distances
}

# calculate lambda
lambdaOpt <- optim(1, lambda_func, dist = distancesMatrix, n = uniqueFrame$num, independency = nrow(presenceFrame)/2)$par

# fit the models with trainFrame
temporalGLMM <- fitme(Y ~ suitableHabitat + averageElevation + (1|year) + Autocovariate, data = trainFrame, family = binomial)
spatiotemporalGLMM <- fitme(Y ~ suitableHabitat + averageElevation + (1|year) + Autocovariate + Matern(1|x+y), data = trainFrame, family = binomial(), fixed = list(nu = 0.5, rho = lambdaOpt))
  
# save the estimated slope values for the baseline GLMM
estimationTemporalGLMM <- c(fold, temporalGLMM$fixef[[2]], temporalGLMM$fixef[[3]], temporalGLMM$fixef[[4]], vcov(temporalGLMM)[2,2], vcov(temporalGLMM)[3,3], vcov(temporalGLMM)[4,4], '\n')
filename <- paste(c(directoryModels, "casestudy/temporalGLMM.txt"), collapse = "")
sink(file = filename, append = TRUE)
cat(estimationTemporalGLMM)
sink()

# save the parameter values for the spatiotemporal GLMM
estimationSpatiotemporalGLMM <- c(fold, spatiotemporalGLMM$fixef[[2]], spatiotemporalGLMM$fixef[[3]], spatiotemporalGLMM$fixef[[4]], vcov(spatiotemporalGLMM)[2,2], vcov(spatiotemporalGLMM)[3,3], vcov(spatiotemporalGLMM)[4,4], '\n')
filename <- paste(c(directoryModels, "casestudy/spatiotemporalGLMM.txt"), collapse = "")
sink(file = filename, append = TRUE)
cat(estimationSpatiotemporalGLMM)
sink()

# store the models
filename <- paste(c(directoryModels, "casestudy/temporalGLMM_", fold, ".rdata"), collapse = "")
save(temporalGLMM, file = filename)

filename <- paste(c(directoryModels, "casestudy/spatiotemporalGLMM_", fold, ".rdata"), collapse = "")
save(spatiotemporalGLMM, file = filename)
