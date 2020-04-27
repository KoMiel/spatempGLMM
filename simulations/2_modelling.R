
# import packages

library(rjson)
library(INLA)
library(boot)
library(foreach)
library(doParallel)

INLA:::inla.dynload.workaround() 
registerDoParallel(cores=25)


# functions

# this function performs a cross-validation for the models that rely on an autocovariate

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

    # split the data in testing and training dataset
    testIndexes <- which(dataFrame$V4 %in% fold,arr.ind=TRUE)
    testFrame <- dataFrame[testIndexes, ]
    trainFrame <- dataFrame[-testIndexes, ]
    
    # dummy variable to ignore the testing data in fitting
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
      effects = list(data.frame(b0 = 1, V5 = trainFrame$V5, V6 = trainFrame$V6, V7 = trainFrame$V7, autocovariate = trainFrame$Autocovariate, year = trainFrame$V4), s = indexs)
    )
      
    # for test data, tag is pred and y is FakeY (NA for all points)
    testStack <- inla.stack(
      tag = "pred",
      data = list(y = testFrame$FakeY),
      A = list(1, ATest),
      effects = list(data.frame(b0 = 1, V5 = testFrame$V5, V6 = testFrame$V6, V7 = testFrame$V7, autocovariate = testFrame$Autocovariate, year = testFrame$V4), s = indexs)
    )
      
    # stk.full has stk.e and stk.p
    stk.full <- inla.stack(trainStack, testStack)
      
    # formulas for all three models
    formulaAuto <- y ~ 0 + b0 + V6 + V7 + f(year, model = "iid") + autocovariate
    formulaComb <- y ~ 0 + b0 + V6 + V7 + f(year, model = "iid") + autocovariate + f(s, model = spde)
    formulaComp <- y ~ 0 + b0 + V6 + V7 + f(year, model = "iid") + autocovariate + V5
    
    # fit all three models
    auto <- inla(formulaAuto,
                 family = "binomial",
                 control.family = list(link = "logit"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)
                 ),
#                 control.inla = list(strategy = 'gaussian'),
                 verbose = FALSE
    )
    
    comb <- inla(formulaComb,
                 family = "binomial",
                 control.family = list(link = "logit"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)
                 ),
#                 control.inla = list(strategy = 'gaussian'),
                 verbose = FALSE
    )

    comp <- inla(formulaComp,
                 family = "binomial",
                 control.family = list(link = "logit"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)
                 ),
#                 control.inla = list(strategy = 'gaussian'),
                 verbose = FALSE
    )
    
    # combine the slopes into one vector
    slopes <- cbind(auto$summary.fixed['V6', 'mean'],
                    auto$summary.fixed['V6', 'sd'],
                    auto$summary.fixed['V7', 'mean'],
                    auto$summary.fixed['V7', 'sd'],
                    auto$summary.fixed['autocovariate', 'mean'],
                    auto$summary.fixed['autocovariate', 'sd'],
                    comb$summary.fixed['V6', 'mean'],
                    comb$summary.fixed['V6', 'sd'],
                    comb$summary.fixed['V7', 'mean'],
                    comb$summary.fixed['V7', 'sd'],
                    comb$summary.fixed['autocovariate', 'mean'],
                    comb$summary.fixed['autocovariate', 'sd'],
                    comp$summary.fixed['V6', 'mean'],
                    comp$summary.fixed['V6', 'sd'],
                    comp$summary.fixed['V7', 'mean'],
                    comp$summary.fixed['V7', 'sd'],
                    comp$summary.fixed['autocovariate', 'mean'],
                    comp$summary.fixed['autocovariate', 'sd'],
                    comp$summary.fixed['V5', 'mean'],
                    comp$summary.fixed['V5', 'sd'])

    # get the predictions for training data
    index <- inla.stack.index(stack = stk.full, tag = "est")$data
    trainFrame$predAuto <- auto$summary.fitted.values[index, "mean"]
    trainFrame$predComb <- comb$summary.fitted.values[index, "mean"]
    trainFrame$predComp <- comp$summary.fitted.values[index, "mean"]
    
    # get the predictions for test data
    index <- inla.stack.index(stack = stk.full, tag = "pred")$data
    testFrame$predAuto <- auto$summary.fitted.values[index, "mean"]
    testFrame$predComb <- comb$summary.fitted.values[index, "mean"]
    testFrame$predComp <- comp$summary.fitted.values[index, "mean"]
    
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

# this function performs a cross-validation for the models that DO NOT rely on an autocovariate
# we split them because these models do not depend on a range parameter: we only need to fit them once

CV2 <- function(dataFrame, folds, lscape, subpath, directoryModels){
  
  # generate directories for model
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
  
  # loop over all folds
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
      effects = list(data.frame(b0 = 1, V6 = trainFrame$V6, V7 = trainFrame$V7, year = trainFrame$V4), s = indexs)
    )
    
    # for test data, tag is pred and y is FakeY (NA for all points)
    testStack <- inla.stack(
      tag = "pred",
      data = list(y = testFrame$FakeY),
      A = list(1, ATest),
      effects = list(data.frame(b0 = 1, V6 = testFrame$V6, V7 = testFrame$V7, year = testFrame$V4), s = indexs)
    )
    
    # stk.full has stk.e and stk.p
    stk.full <- inla.stack(trainStack, testStack)
    
    # formulas for both models
    formulaBase <- y ~ 0 + b0 + V6 + V7 + f(year, model = "iid")
    formulaRand <- y ~ 0 + b0 + V6 + V7 + f(year, model = "iid") + f(s, model = spde)
    
    # fit both models
    base <- inla(formulaBase,
                 family = "binomial",
                 control.family = list(link = "logit"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)
                 ),
#                 control.inla = list(strategy = 'gaussian'),
                 verbose = FALSE
    )
    
    rand <- inla(formulaRand,
                 family = "binomial",
                 control.family = list(link = "logit"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)
                 ),
#                 control.inla = list(strategy = 'gaussian'),
                 verbose = FALSE
    )
    
    # combine the slopes into one vector
    slopes <- cbind(base$summary.fixed['V6', 'mean'],
                    base$summary.fixed['V6', 'sd'],
                    base$summary.fixed['V7', 'mean'],
                    base$summary.fixed['V7', 'sd'],
                    rand$summary.fixed['V6', 'mean'],
                    rand$summary.fixed['V6', 'sd'],
                    rand$summary.fixed['V7', 'mean'],
                    rand$summary.fixed['V7', 'sd'])
    
    # get the predictions for training data
    index <- inla.stack.index(stack = stk.full, tag = "est")$data
    trainFrame$predBase <- base$summary.fitted.values[index, "mean"]
    trainFrame$predRand <- rand$summary.fitted.values[index, "mean"]
    
    # get the predictions for testing data
    index <- inla.stack.index(stack = stk.full, tag = "pred")$data
    testFrame$predBase <- base$summary.fitted.values[index, "mean"]
    testFrame$predRand <- rand$summary.fitted.values[index, "mean"]
    
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

# considered options for p_c
p_c <- c(0.64,0.8,1,1.25,1.5625)

# folds of the cross-validation
folds <- list(c(0,1), c(2,3), c(4,5), c(6,7), c(8,9))

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
  filename <- paste(c(directoryModels, subpath, "/seed.txt"), collapse = "")
  dir.create(path = directoryModels, showWarnings = FALSE) # generate directories if not already existing
  dir.create(path = paste(c(directoryModels, subpath), collapse = ""), showWarnings = FALSE)
  sink(file = filename, append = TRUE)
  cat(randomSeed)
  sink()
  
  # read the data
  filename <- paste(c(directoryDatasets, subpath, "sampled_dataset_", lscape, ".rdata"), collapse = "")
  load(file = filename) #dataFrame

  # perform the cross-validation for different values of p_c (parallel computing)
  resultsAutocov <- foreach(nValue = 1:length(p_c)) %dopar% {
    value <- p_c[nValue]
    dataFrame$Autocovariate <- dataFrame[,paste0("p_c_", toString(value))]
    singleResultAutocov <- CV1(dataFrame, folds, lscape, subpath, directoryModels)
  }
  
  # save the results
  filename <- paste(c(directoryModels, subpath, lscape, "/resultsAutocov.rdata"), collapse = "")
  save(resultsAutocov, file = filename)
  
  # perform the cross-validation for the non-autocovariate models
  resultsNoautocov <- CV2(dataFrame, folds, lscape, subpath, directoryModels)
  
  # save the results
  filename <- paste(c(directoryModels, subpath, lscape, "/resultsNoautocov.rdata"), collapse = "")
  save(resultsNoautocov, file = filename)
}
