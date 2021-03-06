
# this file generates a list of all ecotopes which we then classified as suitable and unsuitable

# import packages

library(sf)
library(rjson)

# open and read settings
f <- "../../settings.json"
settings <- fromJSON(file = f)

# get the directories
directoryCasestudyExternal <- settings$directory$casestudy_external
directoryCasestudy <- settings$directory$casestudy

# read the ecotope maps
dsn <- paste(c(directoryCasestudyExternal, "ecotopenKarteringCycli/ecotopen_herziene_tweede_cyclus/geogegevens/"), collapse = "")
sfdata <- read_sf(dsn = dsn, layer = "vegetatiestructuur_CC2")
ecotopes <- sfdata['ECOTOOP']
ecotopeList <- unique(ecotopes$ECOTOOP)

# get a list of all ecotopes
uniqueEcotopeList <- unique(ecotopeList)

# generate directories if necessary
dir.create(path = paste(c("../data"), collapse = ""), showWarnings = FALSE)
dir.create(path = paste(c(directoryCasestudy), collapse = ""), showWarnings = FALSE)
dir.create(path = paste(c(directoryCasestudy, "ecotopes/"), collapse = ""), showWarnings = FALSE)

# save the list
filename <- paste(c(directoryCasestudy, "ecotopes/uniqueEcotopeList.csv"), collapse = "")
write.csv(uniqueEcotopeList, filename)