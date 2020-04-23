
# import packages

library(rjson)
library(INLA)
library(boot)
library(foreach)
library(doParallel)

INLA:::inla.dynload.workaround() 
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
  
  # loop over years
  for (k in (min(year)):(max(year))) {

    # filter all observations of one sampling event
    Y_Year <- Y[year == k] # the presence/absence encoding
    coordsYear <- coords[year == k, ] # the position
    presCoordsYear <- coordsYear[Y_Year == 1, ] # positions, but only presences
    
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
  
  # loop over year
  for (k in (min(year)):(max(year))) {
    
    # get the matrix that belongs to the sampling event
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

# this function performs a cross-validation for the models that rely on an autocovariate

CV1 <- function(dataFrame, folds, distances, directoryModels, parameter){
  
  # generate directory for the models
  dir.create(path = paste(c(directoryModels, 'casestudy/'), collapse = ""), showWarnings = FALSE)

  # generate lists for results
  predictionTrain <- list()
  predictionTest <- list()
  parameterFrame <- list()
  
  # generate INLA mesh, spde objects 
  boundary <- cbind(c(min(dataFrame$x),min(dataFrame$x),max(dataFrame$x),max(dataFrame$x)), c(min(dataFrame$y),max(dataFrame$y),max(dataFrame$y),min(dataFrame$y)))
  mesh <- inla.mesh.2d(boundary, max.edge = 5, offset = c(1.5,3), cutoff = 3)
  
  # alpha = 1.5 -> exponential structure
  spde <- inla.spde2.matern(mesh = mesh, alpha = 1.5, constr = TRUE)
  indexs <- inla.spde.make.index("s", spde$n.spde)
  
  # calculate the autocovariate
  dataFrame[,'Autocovariate'] <- autocovFunc(distances = distances, p = parameter, Y = dataFrame['Y'], year = dataFrame['year'])

  # loop over all folds (parallel computing)
  res <- foreach(nFold = 1:length(folds)) %dopar% {
    
    # get fold
    fold <- folds[nFold][[1]]
    
    # split the corncrake data in testing and training dataset
    testIndexes <- which(dataFrame$year %in% fold,arr.ind=TRUE)
    testFrame <- dataFrame[testIndexes, ]
    trainFrame <- dataFrame[-testIndexes, ]
    
    # dummy variable to not use the testing data in fitting
    testFrame$FakeY <- NA
    
    # generate A for training data and test data
    coordsTrain <- cbind(trainFrame$x, trainFrame$y)
    A <- inla.spde.make.A(mesh = mesh, loc = coordsTrain)
    coordsTest <- cbind(testFrame$x, testFrame$y)
    ATest <- inla.spde.make.A(mesh = mesh, loc = coordsTest)
    
    # generate stack for training data with data and effects
    trainStack <- inla.stack(
      tag = "est",
      data = list(y = trainFrame$Y),
      A = list(1, A),
      effects = list(data.frame(b0 = 1, habitat = trainFrame$suitableHabitat, elevation = trainFrame$averageElevation, autocovariate = trainFrame$Autocovariate, year = trainFrame$year), s = indexs)
    )
    
    # for test data, tag is pred and y is FakeY (NA for all points)
    testStack <- inla.stack(
      tag = "pred",
      data = list(y = testFrame$FakeY),
      A = list(1, ATest),
      effects = list(data.frame(b0 = 1, habitat = testFrame$suitableHabitat, elevation = testFrame$averageElevation, autocovariate = testFrame$Autocovariate, year = testFrame$year), s = indexs)
    )
    
    # stk.full has stk.e and stk.p
    stk.full <- inla.stack(trainStack, testStack)

    # formulas for both models
    formulaAuto <- y ~ 0 + b0 + habitat + elevation + f(year, model = "iid") + autocovariate
    formulaComb <- y ~ 0 + b0 + habitat + elevation + f(year, model = "iid") + autocovariate + f(s, model = spde)

    # fit the models
    auto <- inla(formulaAuto,
                 family = "binomial",
                 control.family = list(link = "logit"),
                 data = inla.stack.data(stk.full),
                 control.predictor = list(
                   compute = TRUE, link = 1,
                   A = inla.stack.A(stk.full)
                 ),
#                control.inla = list(strategy = 'gaussian'),
#                 verbose = TRUE
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
#                 verbose = TRUE
    )
    
    # combine the slopes into one vector
    slopes <- cbind(auto$summary.fixed['habitat', 'mean'] * sd(trainFrame[,'suitableHabitat']),
                    auto$summary.fixed['habitat', 'sd'] * sd(trainFrame[,'suitableHabitat']),
                    auto$summary.fixed['elevation', 'mean'] * sd(trainFrame[,'averageElevation']),
                    auto$summary.fixed['elevation', 'sd'] * sd(trainFrame[,'averageElevation']),
                    auto$summary.fixed['autocovariate', 'mean'] * sd(trainFrame[,'Autocovariate']),
                    auto$summary.fixed['autocovariate', 'sd'] * sd(trainFrame[,'Autocovariate']),
                    comb$summary.fixed['habitat', 'mean'] * sd(trainFrame[,'suitableHabitat']),
                    comb$summary.fixed['habitat', 'sd'] * sd(trainFrame[,'suitableHabitat']),
                    comb$summary.fixed['elevation', 'mean'] * sd(trainFrame[,'averageElevation']),
                    comb$summary.fixed['elevation', 'sd'] * sd(trainFrame[,'averageElevation']),
                    comb$summary.fixed['autocovariate', 'mean'] * sd(trainFrame[,'Autocovariate']),
                    comb$summary.fixed['autocovariate', 'sd'] * sd(trainFrame[,'Autocovariate']))

    # get the predictions for training data
    index <- inla.stack.index(stack = stk.full, tag = "est")$data
    trainFrame$predAuto <- auto$summary.fitted.values[index, "mean"]
    trainFrame$predComb <- comb$summary.fitted.values[index, "mean"]

    # get the predictions for test data
    index <- inla.stack.index(stack = stk.full, tag = "pred")$data
    testFrame$predAuto <- auto$summary.fitted.values[index, "mean"]
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

# this function performs a cross-validation for the models that DO NOT rely on an autocovariate
# we split them because these models do not depend on a range parameter: we only need to fit them once

CV2 <- function(dataFrame, folds, directoryModels){
  
  # generate directory for models
  dir.create(path = paste(c(directoryModels, 'casestudy/'), collapse = ""), showWarnings = FALSE)
  
  # generate lists for results
  predictionTrain <- list()
  predictionTest <- list()
  parameterFrame <- list()
  
  # generate INLA mesh, spde objects 
  boundary <- cbind(c(min(dataFrame$x),min(dataFrame$x),max(dataFrame$x),max(dataFrame$x)), c(min(dataFrame$y),max(dataFrame$y),max(dataFrame$y),min(dataFrame$y)))
  mesh <- inla.mesh.2d(boundary, max.edge = 5, offset = c(1.5,3), cutoff = 3)
  
  # alpha = 1.5 -> exponential structure
  spde <- inla.spde2.matern(mesh = mesh, alpha = 1.5, constr = TRUE)
  indexs <- inla.spde.make.index("s", spde$n.spde)
  
  # loop over all folds
  res <- foreach(nFold = 1:length(folds)) %dopar% {
    
    # get fold
    fold <- folds[nFold][[1]]
    
    # split the corncrake data in testing and training dataset
    testIndexes <- which(dataFrame$year %in% fold,arr.ind=TRUE)
    testFrame <- dataFrame[testIndexes, ]
    trainFrame <- dataFrame[-testIndexes, ]
    
    # dummy variable to not use the testing data in fitting
    testFrame$FakeY <- NA
    
    # generate A for training data and test data
    coordsTrain <- cbind(trainFrame$x, trainFrame$y)
    A <- inla.spde.make.A(mesh = mesh, loc = coordsTrain)
    coordsTest <- cbind(testFrame$x, testFrame$y)
    ATest <- inla.spde.make.A(mesh = mesh, loc = coordsTest)
    
    # generate stack for training data with data and effects
    trainStack <- inla.stack(
      tag = "est",
      data = list(y = trainFrame$Y),
      A = list(1, A),
      effects = list(data.frame(b0 = 1, habitat = trainFrame$suitableHabitat, elevation = trainFrame$averageElevation, year = trainFrame$year), s = indexs)
    )
    
    # for test data, tag is pred and y is FakeY (NA for all points)
    testStack <- inla.stack(
      tag = "pred",
      data = list(y = testFrame$FakeY),
      A = list(1, ATest),
      effects = list(data.frame(b0 = 1, habitat = testFrame$suitableHabitat, elevation = testFrame$averageElevation, year = testFrame$year), s = indexs)
    )
    
    # stk.full has stk.e and stk.p
    stk.full <- inla.stack(trainStack, testStack)
    
    # formulas for both models
    formulaBase <- y ~ 0 + b0 + habitat + elevation + f(year, model = "iid")
    formulaRand <- y ~ 0 + b0 + habitat + elevation + f(year, model = "iid") + f(s, model = spde)
    
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
#                 verbose = TRUE
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
#                 verbose = TRUE
    )
    
    # combine the slopes into one vector
    slopes <- cbind(base$summary.fixed['habitat', 'mean'] * sd(trainFrame[,'suitableHabitat']),
                    base$summary.fixed['habitat', 'sd'] * sd(trainFrame[,'suitableHabitat']),
                    base$summary.fixed['elevation', 'mean'] * sd(trainFrame[,'averageElevation']),
                    base$summary.fixed['elevation', 'sd'] * sd(trainFrame[,'averageElevation']),
                    rand$summary.fixed['habitat', 'mean'] * sd(trainFrame[,'suitableHabitat']),
                    rand$summary.fixed['habitat', 'sd'] * sd(trainFrame[,'suitableHabitat']),
                    rand$summary.fixed['elevation', 'mean'] * sd(trainFrame[,'averageElevation']),
                    rand$summary.fixed['elevation', 'sd'] * sd(trainFrame[,'averageElevation']))
    
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

# get directories
directoryModels <- settings$directory$models
directoryCasestudy <- settings$directory$casestudy

# folds of the cross-validation
folds <- list(2001, 2002, 2003, 2004, 2005, 2006, 2007)

# considered options for p_c
p_cs <- c(0.1,0.1*1.584893,0.1*1.584893^2,0.1*1.584893^3,0.1*1.584893^4,
          1,1*1.584893,1*1.584893^2,1*1.584893^3,1*1.584893^4,
          10)

# load the corncrake data
filename <- paste(c(directoryCasestudy, "corncrakes/survey2.csv"), collapse = "")
dataFrame <- read.csv(filename)

# filter out all incomplete cases
dataFrame <- dataFrame[complete.cases(dataFrame),]

# generate and use a random seed
seed <- runif(n = 1) * 1000000
set.seed(seed)

# save the random seed 
randomSeed <- c(seed, '\n')
filename <- paste(c(directoryModels, "casestudy/seed.txt"), collapse = "")
dir.create(path = directoryModels, showWarnings = FALSE) # generate directories if not already existing
dir.create(path = paste(c(directoryModels, "casestudy/"), collapse = ""), showWarnings = FALSE)
sink(file = filename, append = TRUE)
cat(randomSeed)
sink()

# order by year
dataFrame <- dataFrame[order(dataFrame$year),]

# calculate distances object
distances <- dist_list(Y = dataFrame[,"Y"], coords = dataFrame[,c("x", "y")], year = dataFrame[,"year"])

# perform cross-validation for different values of p_c
resultsAutocov <- foreach(nValue = 1:length(p_cs)) %dopar% {
  value <- p_cs[nValue]
  singleResultAutocov <- CV1(dataFrame, folds = folds, distances = distances, directoryModels = directoryModels, parameter = value)
}

# save the results
filename <- paste(c(directoryModels, "casestudy/resultsAutocov.rdata"), collapse = "")
save(resultsAutocov, file = filename)

# perform cross-validation for models without autocovariate
resultsNoautocov <- CV2(dataFrame, folds, directoryModels)

# save the results
filename <- paste(c(directoryModels, "casestudy/resultsNoautocov.rdata"), collapse = "")
save(resultsNoautocov, file = filename)
