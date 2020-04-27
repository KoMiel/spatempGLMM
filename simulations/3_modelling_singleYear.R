
# import packages

library(rjson)
library(INLA)
library(boot)
library(foreach)
library(doParallel)

INLA:::inla.dynload.workaround()
registerDoParallel(cores=10)

# functions

# this function performs a cross-validation for the Comb model

CV1 <- function(dataFrame, folds, lscape, subpath, directoryModels){

  # generate directories for models
  dir.create(path = paste(c(directoryModels, subpath), collapse = ""), showWarnings = FALSE)
  dir.create(path = paste(c(directoryModels, subpath, lscape, "/"), collapse = ""), showWarnings = FALSE)

  # generate lists for results
  predictionTrain <- list()
  predictionTest <- list()
  parameterFrame <- list()

  # generate INLA mesh, spde objects 
  boundary <- cbind(c(min(dataFrame$V2),min(dataFrame$V2),max(dataFrame$V2),max(dataFrame$V2)), c(min(dataFrame$V3),max(dataFrame$V3),max(dataFrame$V3),min(dataFrame$V3)))
  mesh <- inla.mesh.2d(boundary, max.edge = 5, offset = c(1.5,3), cutoff = 3)

  # alpha = 1.5 -> exponential structure
  spde <- inla.spde2.matern(mesh = mesh, alpha = 1.5, constr = TRUE)
  indexs <- inla.spde.make.index("s", spde$n.spde)

  # loop over all folds (parallel computing)
  res <- foreach(nFold = 1:length(folds)) %dopar% {
    
    # get fold
    fold <- folds[nFold][[1]]

    # split the corncrake data in testing and training dataset
    testIndexes <- which(dataFrame$V4 %in% fold,arr.ind=TRUE)
    testFrame <- dataFrame[testIndexes, ]
    trainFrame <- dataFrame[-testIndexes, ]
    
    # dummy variable to not use the testing data in fitting
    testFrame$FakeY <- NA

    # generate A for training data and test data
    coordsTrain <- cbind(trainFrame$V2, trainFrame$V3)
    A <- inla.spde.make.A(mesh = mesh, loc = coordsTrain)
    coordsTest <- cbind(testFrame$V2, testFrame$V3)
    ATest <- inla.spde.make.A(mesh = mesh, loc = coordsTest)
    
    # generate stack for training data with data and effects
    trainStack <- inla.stack(
      tag = "est",
      data = list(y = trainFrame$V1),
      A = list(1, A),
      effects = list(data.frame(b0 = 1, V5 = trainFrame$V5, V6 = trainFrame$V6, V7 = trainFrame$V7, autocovariate = trainFrame$Autocovariate), s = indexs)
    )
      
    # for test data, tag is pred and y is FakeY (NA for all points)
    testStack <- inla.stack(
      tag = "pred",
      data = list(y = testFrame$FakeY),
      A = list(1, ATest),
      effects = list(data.frame(b0 = 1, V5 = testFrame$V5, V6 = testFrame$V6, V7 = testFrame$V7, autocovariate = testFrame$Autocovariate), s = indexs)
    )
      
    # stk.full has stk.e and stk.p
    stk.full <- inla.stack(trainStack, testStack)
    
    # formulas of the Comb model (note: Only one sampling event, so no random effect for it here)
    formulaComb <- y ~ 0 + b0 + V6 + V7 + autocovariate + f(s, model = spde)

    # fit the model
    # note: for one of the simulations, the model did not converge. Therefore, we used the simpler gaussian approximation in that one case
    comb <- inla(formulaComb,
                 family = "binomial",
                 control.family = list(link = "logit"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)
                 ),
#                 control.inla = list(strategy = 'gaussian'),
                 verbose = TRUE
    )

    
    # combine the slopes into one vector
    slopes <- cbind(comb$summary.fixed['V6', 'mean'],
                    comb$summary.fixed['V6', 'sd'],
                    comb$summary.fixed['V7', 'mean'],
                    comb$summary.fixed['V7', 'sd'],
                    comb$summary.fixed['autocovariate', 'mean'],
                    comb$summary.fixed['autocovariate', 'sd'])
    
    
    # get the predictions for training data
    index <- inla.stack.index(stack = stk.full, tag = "est")$data
    trainFrame$predComb <- comb$summary.fitted.values[index, "mean"]

    # get the predictions for test data
    index <- inla.stack.index(stack = stk.full, tag = "pred")$data
    testFrame$predComb <- comb$summary.fitted.values[index, "mean"]

    # combine frames and parameters into one list to return in parallel computing
    list(trainFrame, testFrame, slopes)
  }
  
  # take the results of all folds and combine them into seperate objects for training predictions, testing predictions and slopes
  for (list in res) {
    predictionTrain <- rbind(predictionTrain, list[[1]])
    predictionTest <- rbind(predictionTest, list[[2]])
    parameterFrame <- rbind(parameterFrame, list[[3]])
  }

  # return the lists
  return(list(predictionTrain, predictionTest, parameterFrame))
}


# read settings file
f <- "../settings.json"
settings <- fromJSON(file = f)

# read scenario file
f <- "../scenarios.json"
scenarios <- fromJSON(file = f)

# read the model list
modelList <- read.table("../modelList.csv", sep = ',')

# read directories
directoryDatasets <- settings$directory$datasets
directoryModels <- settings$directory$models

# 10-fold cross-validation (all but one year for testing, each)
folds <- list(c(0,1,2,3,4,5,6,7,8),
              c(0,1,2,3,4,5,6,7,9),
              c(0,1,2,3,4,5,6,8,9),
              c(0,1,2,3,4,5,7,8,9),
              c(0,1,2,3,4,6,7,8,9),
              c(0,1,2,3,5,6,7,8,9),
              c(0,1,2,4,5,6,7,8,9),
              c(0,1,3,4,5,6,7,8,9),
              c(0,2,3,4,5,6,7,8,9),
              c(1,2,3,4,5,6,7,8,9))

# random seeds
seeds <- runif(n = nrow(modelList)) * 1000000

# loop over all scenarios
for (row in 1:nrow(modelList)) {
  
  # print message (progress keeping)
  print(row)
  
  # get scenario and grid from list
  scenario <- modelList[row, 1]
  lscape <- modelList[row, 2]

  # get the subpath for scenario
  subpath <- scenarios$scenarios[[scenario]]$path

  # use the random seed
  seed <- seeds[row]
  set.seed(seed)
  
  # save the random seed 
  randomSeed <- c(seed, lscape, '\n')
  filename <- paste(c(directoryModels, subpath, "/SingleYearSeed.txt"), collapse = "")
  dir.create(path = directoryModels, showWarnings = FALSE) # generate directories if not already existing
  dir.create(path = paste(c(directoryModels, subpath), collapse = ""), showWarnings = FALSE)
  sink(file = filename, append = TRUE)
  cat(randomSeed)
  sink()
  
  # read the data
  filename <- paste(c(directoryDatasets, subpath, "sampled_dataset_", lscape, ".rdata"), collapse = "")
  load(file = filename) #dataFrame

  # set range parameter to 1
  value <- 1
  dataFrame$Autocovariate <- dataFrame[,paste0("p_c_", toString(value))]
  
  # perform cross-validation
  resultsAutocov <- CV1(dataFrame, folds, lscape, subpath, directoryModels)

  # save the results  
  filename <- paste(c(directoryModels, subpath, lscape, "/resultsOneyearAutocov.rdata"), collapse = "")
  save(resultsAutocov, file = filename)
}
