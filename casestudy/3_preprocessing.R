
# this file performs random sampling of absences and calculation of environmental site properties

# import packages

library(raster)
library(sf)
library(rjson)
library(dplyr)


# functions

# this function rounds to a specified base (e.g. to 5)
# necessary because the corncrake observations are only given to 5 meters accuracy

mround <- function(x,base){
  base*round(x/base)
}

# list of all years with observations
years <- c(2001,2002,2003,2004,2005,2006,2007)

# open and read settings
f <- "../settings.json"
settings <- fromJSON(file = f)

# get the directories
directoryCasestudyExternal <- settings$directory$casestudy_external
directoryCasestudy <- settings$directory$casestudy

# get the home range diameter
homerangeDiameter <- settings$homerangeDiameter

# get the number of pseudo-absences
numAbsences <- settings$absNum

# generate a random seed and use it
seed <- runif(n = 1) * 1000000
set.seed(seed)

# save the random seed 
randomSeed <- c(seed, '\n')
filename <- paste(c(directoryCasestudy, "corncrakes/prep_seed.txt"), collapse = "")
dir.create(path = directoryCasestudy, showWarnings = FALSE) # generate directories if not already existing
dir.create(path = paste(c(directoryCasestudy, "corncrakes/"), collapse = ""), showWarnings = FALSE)
sink(file = filename, append = TRUE)
cat(randomSeed)
sink()

# open the corncrake observations
dsn <- paste(c(directoryCasestudyExternal, "corncrakes/"), collapse = "")
presences <- read_sf(dsn = dsn, layer = "alle roepplekken 2001-2007")

# some observations don't have a date, replace it with a date outside of the surveys (handling of NA)
presences$DATUM[is.na(presences$DATUM)] <- as.Date(1, origin="1970-01-01")

# only take those that are during the second survey (either because they were during the surveys or because they were flagged as survey observations)
survey2_1 <- presences[presences$DATUM == '2001-06-29' |
                       presences$DATUM == '2001-06-30' |
                       presences$DATUM == '2001-07-01' |
                       presences$DATUM == '2001-07-02' |
                       presences$DATUM == '2002-06-21' |
                       presences$DATUM == '2002-06-22' |
                       presences$DATUM == '2002-06-23' |
                       presences$DATUM == '2002-06-24' |
                       presences$DATUM == '2003-06-20' |
                       presences$DATUM == '2003-06-21' |
                       presences$DATUM == '2003-06-22' |
                       presences$DATUM == '2003-06-23' |
                       presences$DATUM == '2005-06-24' |
                       presences$DATUM == '2005-06-25' |
                       presences$DATUM == '2005-06-26' |
                       presences$DATUM == '2005-06-27' |
                       presences$DATUM == '2006-06-23' |
                       presences$DATUM == '2006-06-24' |
                       presences$DATUM == '2006-06-25' |
                       presences$DATUM == '2006-06-26' |
                       presences$DATUM == '2007-06-22' |
                       presences$DATUM == '2007-06-23' |
                       presences$DATUM == '2007-06-24' |
                       presences$DATUM == '2007-06-25',]
survey2_2 <- presences[presences$SIMULTAAN == 2,]

# combine both
survey2 <- unique(rbind(survey2_1, survey2_2))
survey2 <- survey2[1]

# open the converted ecotopes
dsn <- paste(c(directoryCasestudy, "ecotopes/"), collapse = "")
sfdata <- read_sf(dsn = dsn, layer = "cycle2_converted")

# manual selection of the case study region
IJssel <- st_crop(sfdata[sfdata$GEBIED == 'IJssel', ], c(xmin = 1, xmax = 1000000, ymin = 1, ymax = 508190))
Lek1 <- sfdata[sfdata$GEBIED == 'Lek', ]
Lek2 <- st_crop(sfdata[sfdata$GEBIED == 'Dordsche Kil, Oude Maas, Spui, Noord, Lek', ], c(xmin = 1, xmax = 1000000, ymin = 433578, ymax = 1000000))
Nederrijn <- sfdata[sfdata$GEBIED == 'Neder-Rijn', ]
Waal1 <- sfdata[sfdata$GEBIED == 'Waal', ]
Waal2 <- st_crop(sfdata[sfdata$GEBIED == 'Boven Merwede, Waal', ], c(xmin = 124275, xmax = 1000000, ymin = 1, ymax = 1000000))

# throw away everything else
data <- do.call(rbind, list(IJssel, Lek1, Lek2, Nederrijn, Waal1, Waal2))
ecotopes <- data['ECOTOOP']

# generate empty dataframes for presences and absences
absencesFrame <- data.frame(matrix(nrow = round(numAbsences*length(years)*1.1), ncol = 6))
names(absencesFrame) <- c("Y", "x", "y", "year", "suitableHabitat", "averageElevation")
presencesFrame <- data.frame(matrix(nrow = nrow(survey2), ncol = 6))
names(presencesFrame) <- c("Y", "x", "y", "year", "suitableHabitat", "averageElevation")

# generate points at the observations
point <- survey2$geometry
  
# get the crs of the data
crs <- st_crs(ecotopes)

# generate the home ranges
homerange <- st_sfc(st_buffer(point, dist = homerangeDiameter))
st_crs(homerange) <- crs
  
# make the points compatible
point <- st_sfc(point)
st_crs(point) <- crs

# open the elevation map
filename <- paste(c(directoryCasestudyExternal, "elevation/demRhineDelta_25m.map"), collapse = "")
elevation_raster <- raster(filename)

# extract the average elevation at all presence positions
elevation <- extract(elevation_raster, survey2, buffer = homerangeDiameter, fun = mean)  

# initialize a counter
counter <- 1

# loop over presences
for (i in 1:length(homerange)) {
  
  # calculate the overlap of homerange with any ecotope
  intersection <- st_intersection(ecotopes, homerange[i])
  area <- sum(st_area(intersection))
  
  # calculate the overlap of homerange with suitable ecotope
  suitableIntersection <- intersection[intersection$ECOTOOP == 1,]
  areaSuitable <- sum(st_area(suitableIntersection))
  
  # if there is no overlap, we cannot divide the two, so we have to exclude that case
  if (as.numeric(area) > 0) {
    areapercentage <- as.numeric(areaSuitable/area)
  } else {
    areapercentage <- NA
  }
  
  # new row in data frame, counter to next
  presencesFrame[counter, ] <- c(1, st_coordinates(survey2)[i,1], st_coordinates(survey2)[i,2], year = survey2$JAAR[i], areapercentage, elevation[i])
  counter <- counter + 1
}


# initialize new counter
counter <- 1

# loop over years
for (nYear in 1:length(years)) {
   
  # get the year
  year <- years[nYear]
  
  # print message
  print(year)

  # sample absences
  sample <- st_sample(ecotopes, size = round(numAbsences**1.1))
  
  # round coordinates to 5 meters
  for (i in 1:length(sample)) {
    sample[[i]][1] <- mround(sample[[i]][1], 5)
    sample[[i]][2] <- mround(sample[[i]][2], 5)
  }

  # get coordinates and extract elevation from elevation raster
  x <- st_coordinates(sample)[,1]
  y <- st_coordinates(sample)[,2]
  elevation <- extract(elevation_raster, data.frame(x = x, y = y), buffer = homerangeDiameter, fun = mean)

  # set the home ranges
  homerange <- st_sfc(st_buffer(sample, dist = homerangeDiameter))
  st_crs(homerange) <- crs
  
  # loop over absences
  for (i in 1:length(homerange)) {
    
    # calculate the overlap of homerange with any ecotope
    intersection <- st_intersection(ecotopes, homerange[i])
    area <- sum(st_area(intersection))
    
    # calculate the overlap of homerange with suitable ecotope
    suitableIntersection <- intersection[intersection$ECOTOOP == 1,]
    areaSuitable <- sum(st_area(suitableIntersection))
    
    # divide
    areapercentage <- as.numeric(areaSuitable/area)
    
    # new row in data frame, counter to next
    absencesFrame[counter, ] <- c(0, x[i], y[i], year, areapercentage, elevation[i])
    counter <- counter + 1
  }
}

# exclude observations where elevation is NA
absencesFrame <- absencesFrame[complete.cases(absencesFrame),]

#subsample the required amount of absences per year
absencesFrame <- absencesFrame %>% group_by(year) %>% sample_n(size = numAbsences)

# combine data
dataFrameSurvey2 <- rbind(presencesFrame, absencesFrame)

# convert coordinates from meters to kilometers
dataFrameSurvey2$x <- dataFrameSurvey2$x/1000
dataFrameSurvey2$y <- dataFrameSurvey2$y/1000

# save the presence-absence data
filename <- paste(c(directoryCasestudy, "corncrakes/survey2.csv"), collapse = "")
write.csv(dataFrameSurvey2, file = filename)