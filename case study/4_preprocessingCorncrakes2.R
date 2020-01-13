

### In this file, we calculate the autocovariate


### import packages ###

library(dplyr)
library(rjson)
library(spaMM)
library(rjson)


### function definitions ###

# this function rounds to a specified base (e.g. to 5)
mround <- function(x,base){
  base*round(x/base)
}


# this function generates a list of matrices (one for each distribution/sampling event)
# with distances of all locations (presence or absence) to all presence locations
# this is then later used as a base in the calculation of the autocovariate

dist_list <- function(y, coords, year){
  # generate an empty list for the distance matrices
  distances <- list()
  # loop over years
  for (k in min(year):max(year)) {
    # filter all observations of one year
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
    distances[[k]] <- dist #add the matrix to the list
  }
  # return the list
  return(distances)
}

# this function calculates the autocovariate from the distanc matrices (output of dist_list)

autocov <- function(distances, year, param = 1){
  # generate a vector
  autocov <- vector('numeric')
  # loop over all years
  for (k in min(year):max(year)) {
    # get the matrix that belongs to the year
    distYear <- distances[[k]]
    # calculate the autocovariate that is given from each of the presences
    weightedDistYear <- exp(- distYear * param)
    # for each presence or absence, add up the conspecific interaction of all presences
    autocovYear <- rowSums(weightedDistYear, na.rm = TRUE)
    # append the result to the vector
    autocov <- append(autocov, autocovYear)
  }
  # standardize
  autocov <- (autocov - mean(autocov))/sd(autocov)
  # return the vector
  return(autocov)
}


### main program ###

# read settings file
f <- "../../settings.json"
settings <- fromJSON(file = f)

# get directory
directoryCasestudy <- settings$directory$casestudy

# load corncrake data
filename <- paste(c(directoryCasestudy, "corncrakes/survey2.csv"), collapse = "")
dataFrame <- read.csv(filename, head=T)

# round the coordinates (because the presences are not exact, but rounded to 5 meters)
dataFrame$x <- mround(dataFrame$x, 5)/1000
dataFrame$y <- mround(dataFrame$y, 5)/1000

# calculate distances object
distances <- dist_list(y = dataFrame['Y'], coords = dataFrame[,3:4], year = dataFrame['year'])

# calculate the autocovariate
dataFrame[,'Autocovariate'] <- autocov(distance = distances, year = dataFrame['year'])

# save the dataframe
filename <- paste(c(directoryCasestudy, "corncrakes/survey2Frame.csv"), collapse = "")
write.csv(dataFrame, filename)