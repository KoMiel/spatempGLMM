

### this file does random Sampling of pseudo-absences and calculation of abiotic site properties


### import packages ###

library(raster)
library(sf)
library(rjson)


# list of all years with observations
years <- c(2001,2002,2003,2004,2005,2006,2007)

# open and read settings
f <- "../../settings.json"
settings <- fromJSON(file = f)

# get the directories
directoryCasestudyExternal <- settings$directory$casestudy_external
directoryCasestudy <- settings$directory$casestudy

# get the home range diameter
homerangeDiameter <- settings$homerangeDiameter

# get the number of pseudo-absences
numAbsences <- settings$absNum

# generate a random seed and store it for reproducibility
seed <- runif(n = 1, min = 0, max = 1000000)
set.seed(seed)
filename <- paste(c(directoryCasestudy, "corncrakes/seed.txt"), collapse = "")
dir.create(path = paste(c(directoryCasestudy, "corncrakes/"), collapse = ""), showWarnings = FALSE)
write.table(seed, file = filename, row.names = FALSE, col.names = FALSE, append = TRUE)


### main program ###

# open the corncrake observations
dsn <- paste(c(directoryCasestudyExternal, "corncrakes/"), collapse = "")
presences <- read_sf(dsn = dsn, layer = "alle roepplekken 2001-2007")

# some observations don't have a date, replace it with a data outside of the surveys (handling of NA)
presences$DATUM[is.na(presences$DATUM)] <- as.Date(1, origin="1970-01-01")

# only take those that are either during the first survey or the second survey (either because they were during the surveys or because they were flagged as survey observations)
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

# generate an empty dataframe for the survey
dataFrameSurvey2 <- data.frame(x = numeric(), y = numeric(), suitableHabitat = numeric(), averageHeight = numeric(), distanceToRiver = numeric(), jaggedness = numeric(), suitableJaggedness = numeric(), year = numeric())

# get the crs of the data
crs <- st_crs(ecotopes)

# open the elevation map
filename <- paste(c(directoryCasestudyExternal, "elevation/demRhineDelta_25m.map"), collapse = "")
elevation_raster <- raster(filename)

# extract the average elevation at all presence positions
elevationFrame2 <- extract(elevation_raster, survey2, buffer = homerangeDiameter, fun = mean)  

#loop over all observations from survey 2
for (j in 1:nrow(survey2)) {
  
  # progress control
  print(j)
  
  # generate a point at the observation
  point <- survey2$geometry[j]
  
  # generate the home range
  homerange <- st_sfc(st_buffer(point, dist = homerangeDiameter))
  st_crs(homerange) <- crs
  
  # make the point compatible
  point <- st_sfc(point)
  st_crs(point) <- crs
  
  # calculate the suitable habitat proportion
  intersection <- st_intersection(ecotopes, homerange)
  area <- sum(st_area(intersection))
  suitableIntersection <- intersection[intersection$ECOTOOP == 1,]
  areaSuitable <- sum(st_area(suitableIntersection))
  
  # if there is no overlap, the division cannot be carried out, so we have to exclude that case
  if (as.numeric(area) > 0) {
    areapercentage <- as.numeric(areaSuitable/area)
  } else {
    areapercentage <- NA
  }

  # generate a row for the data frame with presence/absence-indicator, coordinates, suitable habitat proportion, elevation and year of observation
  newRow <- data.frame(Y = 1, x = st_coordinates(survey2)[j,1], y = st_coordinates(survey2)[j,2], suitableHabitat = areapercentage, averageElevation = elevationFrame2[j], year = survey2$JAAR[j])
  dataFrameSurvey2 <- rbind(dataFrameSurvey2, newRow)
}


#######################

# loop over years and the required number of absences (we sample the same number per year)
for (year in years) {
  for (j in 1:numAbsences) {
    # keep sampling until a valid one has been found (this is necessary because the sample function used here does not necessarily provides an obervation from the case study region)
    finished <- FALSE
    while (finished == FALSE) {
      sample <- st_sample(ecotopes, size = 1)
      while (nrow(st_coordinates(sample)) == 0) {
        sample <- st_sample(ecotopes, size = 1)
      }
      
      # generate the home range
      homerange <- st_sfc(st_buffer(sample[1], dist = homerangeDiameter))
      st_crs(homerange) <- crs
      
      # calculate the suitable habitat proportion
      intersection <- st_intersection(ecotopes, homerange)
      area <- sum(st_area(intersection))
      suitableIntersection <- intersection[intersection$ECOTOOP == 1,]
      areaSuitable <- sum(st_area(suitableIntersection))
      
      # if there is overlap, then we have found a valid position. If not, then try again
      if (as.numeric(areaSuitable) > 0) {
        areapercentage <- as.numeric(areaSuitable/area)
        finished <- TRUE
      } else {
        next
      }
      
      # convert coordinates to simple features to extract elevation
      x <- st_coordinates(sample)[1,1]
      y <- st_coordinates(sample)[1,2]
      elevation <- extract(elevation_raster, data.frame(x = x, y = y), buffer = homerangeDiameter, fun = mean)
      
      # if there is no elevation data, then sample again
      if(is.na(elevation)) {
        finished <- FALSE
        next
      }
    }
    
    # once we have ensured our position is valid, add a row to the data frame
    newRow <- data.frame(Y = 0, x = x, y = y, suitableHabitat = areapercentage, averageElevation = elevatio, year = year)
    dataFrameSurvey2 <- rbind(dataFrameSurvey2, newRow)
    
    # progress control
    print(j)
  }
}

# save the presence-(pseudo-)absence data
filename <- paste(c(directoryCasestudy, "corncrakes/survey2.csv"), collapse = "")
write.csv(dataFrameSurvey2, file = filename)