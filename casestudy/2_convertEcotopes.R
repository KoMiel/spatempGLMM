
# this file converts ecotope map to suitable and unsuitable

# import packages

library(sf)
library(plyr)
library(rjson)

# open and read settings
f <- "../../settings.json"
settings <- fromJSON(file = f)

# get the directories
directoryCasestudyExternal <- settings$directory$casestudy_external
directoryCasestudy <- settings$directory$casestudy

# read the ecotopes
dsn <- paste(c(directoryCasestudyExternal, "ecotopenKarteringCycli/ecotopen_herziene_tweede_cyclus/geogegevens/"), collapse = "")
sfdata <- read_sf(dsn = dsn, layer = "vegetatiestructuur_CC2")

# read the ecotope list with suitable/unsuitable-classification (1: suitable, 0: unsuitable)
filename <- paste(c(directoryCasestudyExternal, "classifiedUniqueEcotopeList.csv"), collapse = "")
classifiedUniqueEcotopeList <- read.csv(filename)

# only keep the needed variables
ecotopes <- sfdata[c('ECOTOOP', 'GEBIED')]

# get a list of all ecotopes
ecotopeList <- unique(ecotopes$ECOTOOP)

# loop over all ecotopes and replace them with 1 and 0 (suitable and unsuitable)
for (ecotope in ecotopeList) {
  ecotopes$ECOTOOP <- mapvalues(ecotopes$ECOTOOP, from = ecotope, to = classifiedUniqueEcotopeList[which(classifiedUniqueEcotopeList$x == ecotope, arr.ind = TRUE),3])
}

# save the converted ecotope maps
dsn <- paste(c(directoryCasestudy, "ecotopes/cycle2_converted.shp"), collapse = "")
st_write(obj = ecotopes, dsn = dsn, delete_layer = TRUE)