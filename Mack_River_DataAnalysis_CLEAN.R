library(data.table) 
library(dataRetrieval)
library(ggplot2)
library(scales)
library(ggpubr)
library(maps)
library(lubridate)
library(GGally)
library(raster)
library(rnaturalearth)
library(dplyr)
library(plyr)
library(zoo)
library(broom)
library(rstatix)
library(rcartocolor)
library(npreg)



### SET DIRECTORIES ####
## CHANGE THESE TO YOUR OWN COMPUTER
# Set root directory
wd_root <- "/Users/jordanfields/Documents/Graduate_School/Research/Arctic_Meander_Rate/Data_Analysis/October2023_Results/"

# Imports folder (store all import files here)
wd_imports <- paste0(wd_root,"GEE_Exports/")
# Exports folder (save all figures, tables here)
wd_exports <- paste0(wd_root,"R_Exports/")
# Study Sites 
wd_sites <- paste0(wd_root,"Study_Sites/")

wd_figures <- paste0(wd_exports, "Figures/")


# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_imports, wd_exports, wd_figures)
for(i in 1:length(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}

#############################################################################
##                           Add Themes                                    ##
#############################################################################


# Make a plotting theme
theme_pubr <- theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0))
# Make a plot label theme
fancy_scientific_modified <- function(l) { 
  # turn in to character string in scientific notation 
  if(abs(max(log10(l), na.rm = T) - min(log10(l), na.rm = T)) > 2 | 
     # min(l, na.rm = T) < 0.01 | 
     max(l, na.rm = T) > 1e4){ 
    l <- log10(l)
    label <- parse(text = paste("10^",as.character(l),sep = ""))
  }else{
    label <- parse(text = paste(as.character(l), sep = ""))
  }
  # print(label)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(label)
}

my_label <- function(x) {
  # Wrap var names
  names(x)[[1]] <- stringr::str_wrap(gsub("_", " ", names(x)[[1]]), 10)
  # Wrap value labels
  x[[1]] <- stringr::str_wrap(gsub("_", " ", x[[1]]), 10)
  # Call label both with sep "\n"
  label_value(x, sep = "\n")
}

#############################################################################
##                            Data Organization                            ##
#############################################################################

############ FILE AGGREGATION #############
### COMBINE ALL DATA FROM VARIOUS THIESSEN BLOCKS (groups necessary for GEE exports to proceed within Google's user memory limits)


## AREA CHANGED
wd_AreaChange <- paste0(wd_imports, "AreaChanged/")
setwd(wd_AreaChange)
# Import Area Changed (Land Lost/Gained) data for each batch of Thiessens (study sites)
# Get filenames of each batch of Thiessens
AreaChange_data_files <- list.files(pattern = 'AreaChanged_Thiessens_', wd_AreaChange) ## Pattern: string to look for and working directory to search within
AreaChange_import <- rbindlist(  		## TABLE NAME - Update
   lapply(paste0(wd_AreaChange, AreaChange_data_files),  ## Folder name, list name
   fread),
   use.names = T, fill = T)[,':='(.geo = NULL)]


## TOTAL AREA
wd_TotalArea <- paste0(wd_imports, "TotalArea/")
setwd(wd_TotalArea)
# Import Area Changed (Land Lost/Gained) data for each batch of Thiessens (study sites)
# Get filenames of each batch of Thiessens
TotalArea_data_files <- list.files(pattern = 'TotalArea_Thiessens_', wd_TotalArea) ## Pattern: string to look for and working directory to search within
TotalArea_import <- rbindlist(  		## TABLE NAME - Update
  lapply(paste0(wd_TotalArea, TotalArea_data_files),  ## Folder name, list name
         fread),
  use.names = T, fill = T)[,':='(.geo = NULL)]


## NDVI
wd_NDVI <- paste0(wd_imports, "NDVI/")
setwd(wd_NDVI)
# Import Area Changed (Land Lost/Gained) data for each batch of Thiessens (study sites)
NDVI_all <- fread('NDVI_all_1.csv')
#mean NDVI over 5 year windows, like all other variables
NDVI_5_roll <- NDVI_all[, list(NDVI_5yr = frollmean(NDVI, 5, na.rm =T)),
                               by = Thiessen_ID
                        ]
# order both tables by year and Thiessen ID 
NDVI_all<- NDVI_all[order(Thiessen_ID, Year),] 
NDVI_5_roll<- NDVI_5_roll[order(Thiessen_ID),]
#join year column from NDVI_all to NDVI_5_roll 
NDVI_5_roll <-cbind(NDVI_5_roll,NDVI_all$Year)
#remove NA rows: all those before 1990
NDVI_5_roll <- na.omit(NDVI_5_roll)
#actually have 1989 data, so remove that year
NDVI_5_roll <- NDVI_5_roll[V2>1989,]
NDVI_5_roll <- setnames(NDVI_5_roll, 'V2', 'Year')
#create baseline data
# Filter out for 1990-2006 baseline calculation
NDVI_90_06 <- NDVI_5_roll[Year<=2006,]
# Calculate mean/sd from 1975-1995 
# remember that these metrics are for 5year periods, so to get rough idea of annual avg, divide by 5
Baseline_NDVI <- NDVI_90_06[, list(mean_NDVI_90_06 = mean(NDVI_5yr),
                                   sd_NDVI_90_06 = sd(NDVI_5yr)),
                            by = Thiessen_ID]

# Add mean/sd for each segment buffer into the master table
NDVI_5_roll <- NDVI_5_roll[Baseline_NDVI,on=.(Thiessen_ID = Thiessen_ID)]
#calc anomaly
NDVI_5_roll <- NDVI_5_roll[, NDVI_Anomaly := (NDVI_5yr - mean_NDVI_90_06)/sd_NDVI_90_06]
hist(NDVI_5_roll$NDVI_Anomaly) #just making sure numbers look reasonable. They do. 

NDVI_Summary <- NDVI_5_roll[, list(mean_NDVI_Anomaly = mean(NDVI_Anomaly), 
                                   mean_NDVI = mean(NDVI_5yr),
                                   se_NDVI_Anomaly = sd(NDVI_Anomaly)/sqrt(length(NDVI_Anomaly)),
                                   se_NDVI = sd(NDVI_5yr)/sqrt(length(NDVI_5yr))
                                   ),
                               by = Year
                            ]


## TEMPERATURE
wd_Temp <- paste0(wd_imports, "Temperature/")
setwd(wd_Temp)
# Import Area Changed (Land Lost/Gained) data for each batch of Thiessens (study sites)
# Get filenames of each batch of Thiessens
Temperature_data_files <- list.files(pattern = 'Temp_Thiessens_', wd_Temp) ## Pattern: string to look for and working directory to search within
Temperature_import <- rbindlist(  		## TABLE NAME - Update
  lapply(paste0(wd_Temp, Temperature_data_files),  ## Folder name, list name
         fread),
  use.names = T, fill = T)[,':='(.geo = NULL)]
Temperature_import <- unique(Temperature_import) # there are only 109806 entries, should be 41yrs*2717 Thiesses = 111397
Temp_Thiessens <- data.table(unique(Temperature_import$Thiessen_ID)) #have all the thiessens, 2717. Must be missing some years for certain Thiessens, oh well 

## FUTURE TEMPERATURE
# CMIP6 Ensemble Mean 
wd_TempFuture <- paste0(wd_imports, "Future_Temperature/Ensemble_Mean/")
setwd(wd_TempFuture)
# Import Future Temp data (CMIP6 Ensemble Mean) for each batch of Thiessens (study sites)
# Get filenames of each batch of Thiessens
CMIP6_data_files <- list.files(pattern = 'FutureTemp_CMIP6_', wd_TempFuture) ## Pattern: string to look for and working directory to search within
CMIP6_import <- rbindlist(  		## TABLE NAME - Update
  lapply(paste0(wd_TempFuture, CMIP6_data_files),  ## Folder name, list name
         fread),
  use.names = TRUE, fill = TRUE)[,':='(.geo = NULL)]
CMIP6_import <- unique(CMIP6_import) # there are only 109806 entries, should be 41yrs*2717 Thiesses = 111397
FutureTemp <- data.table(unique(CMIP6_import$Thiessen_ID))
fwrite(FutureTemp, paste0(wd_exports, "Future_Temperature.csv"))


#CanESM5 -- model with little bias and tuned for Canada specifically
wd_TempFuture1 <- paste0(wd_imports, "Future_Temperature/CanESM5/")
setwd(wd_TempFuture1)
# if you use this, need to update the  references in the creation of 'Future_Summary" tables later 
FutureTemp_data_files <- fread('FutureTemp_CMIP6_all.csv') 
fwrite(FutureTemp_data_files, paste0(wd_exports, "Future_Temperature.csv"))

## Discharge data -- averaged from three stations across the Mackenzie River Basin and normalized by the pre-1996 Q average, like temp 
#raw Q data is in Carl's data sheet from 08_02_23, in the same root directory as the other source files used in this script
setwd(wd_root)
#Q_Mack <- fread(paste0(wd_root, "Mack_5yr_Mean_Normalized_Q.csv"))
Q_Mack <- fread(paste0(wd_root, "Discharge/MackenzieRiver_Discharge.csv"))
Q_Mack_max <- fread(paste0(wd_root, "Discharge/MackenzieRiver_DischargeMAX.csv"))
Q_Mack_max <- Q_Mack_max[, c('Year', 'STATION_NUMBER', 'STATION_NAME', 'Q_cms')] # select columns
#merge tables 
Q_Mack <- Q_Mack[Q_Mack_max, on= .(Year = Year, STATION_NAME = STATION_NAME, STATION_NUMBER = STATION_NUMBER)]
#set column names 
Q_Mack <- setnames(Q_Mack, c('Q_cms', "i.Q_cms"), c('Q_MEAN_cms', 'Q_MAX_cms'))
#remove unnecessary columns
Q_Mack[, c("V1", "Parameter", "Sum_stat") := NULL]

Q_RedRiver <- fread(paste0(wd_root, "Discharge/Mackenzie_ArcticRedRiver_2020-2023.csv"))
Q_RedRiver <- Q_RedRiver[, Q_MEAN_cms := mean(discharge, na.rm = T), by = 'Year']
Q_RedRiver <- Q_RedRiver[, Q_MAX_cms := max(discharge, na.rm = T), by = 'Year'] 
Q_RedRiver <- Q_RedRiver[, c('STATION_NUMBER', "STATION_NAME", "LATITUDE", "LONGITUDE", "Year", "Q_MEAN_cms", "Q_MAX_cms")]
Q_RedRiver <- unique(Q_RedRiver)
Q_RedRiver <- Q_RedRiver[STATION_NUMBER != 'NA'] # get rid of an extra row which appears because Station Name is missing in a few rows 
Q_Mack <- setcolorder(Q_Mack, c('STATION_NUMBER', "STATION_NAME", "LATITUDE", "LONGITUDE", "Year", "Q_MEAN_cms", "Q_MAX_cms"))
Q_Mack <- rbind(Q_Mack, Q_RedRiver)

#Calculate rolling average over 5 year window for each year at each station
#normalize the Q at each station by that station's long term mean (using full record at each station)
Q_Summary1 <- Q_Mack[, Q_MEAN_norm := Q_MEAN_cms/(mean(Q_MEAN_cms, na.rm=T)), by = STATION_NAME]
Q_Summary1 <- Q_Summary1[, Q_MAX_norm := Q_MAX_cms/(mean(Q_MAX_cms, na.rm=T)), by = STATION_NAME] 
#take mean in each year now across all stations in each year
Q_Summary <- Q_Summary1[, list(Q_Mean_norm_Annual = mean(Q_MEAN_norm, na.rm=T), 
                               Q_Max_norm_Annual = mean(Q_MAX_norm, na.rm=T)), 
                          by = Year
                        ]

Q_Summary <- Q_Summary[order(Year),]
Q_yrs <- as.data.table(list(1939:2023)) #make a full list of years
setnames(Q_yrs, c("V1"), c("Year"))
Q_Summary <- left_join(Q_yrs, Q_Summary, by = "Year") # make Q_Summary have rows for every year, even those missing Q data originally
Q_5yr <- data.table(frollmean(Q_Summary[, Q_Max_norm_Annual], 5, na.rm=T)) # taking rolling average of annual max Q over 5 year window
Q_Mean_5yr <- data.table(frollmean(Q_Summary[, Q_Mean_norm_Annual], 5, na.rm=T))
setnames(Q_5yr,  c("V1"), c("Q_MAX_5yr"))
setnames(Q_Mean_5yr,  c("V1"), c("Q_MEAN_5yr")) 
Q_5yr <- Q_5yr[, Year := list(1939:2023)] # add a year column
Q_Mean_5yr <- Q_Mean_5yr[, Year := list(1939:2023)] # add a year column
Q_Summary <- Q_5yr[Q_Summary, on=.(Year = Year)] # join tables
Q_Summary <- Q_Summary[Q_Mean_5yr, on=.(Year = Year)] 


#calculate baseline data for anomalies, using 1995 and before -- same ref period as the temperature data
Q_ref <- data.table(Q_Summary[Year <= 1995]) 
Q_ref <- Q_ref[, Q_MEAN_5yr := ifelse(Q_MEAN_5yr == "NaN", NA, Q_MEAN_5yr)]
#mean and standard deviation at each station in ref period
Q_baseline <- Q_ref[, Q_MAX_baseline := mean(Q_MAX_5yr, na.rm=T)]
Q_baseline <- Q_ref[, Q_MAX_baseline_SD := sd(Q_MAX_5yr, na.rm=T)]
Q_baseline <- Q_ref[, Q_MEAN_baseline := mean(Q_MEAN_5yr, na.rm=T)]
Q_baseline <- Q_ref[, Q_MEAN_baseline_SD := sd(Q_MEAN_5yr, na.rm=T)]
#grab the columns we want and make the table simple since all years have same value
Q_ref <- Q_ref[1,c('Q_MAX_baseline', 'Q_MAX_baseline_SD', 'Q_MEAN_baseline', 'Q_MEAN_baseline_SD')]
#add baseline values to Q_Summary
Q_Summary <- Q_Summary[, Q_MAX_baseline := Q_ref$Q_MAX_baseline]
Q_Summary <- Q_Summary[, Q_MAX_baseline_SD := Q_ref$Q_MAX_baseline_SD]
Q_Summary <- Q_Summary[, Q_MEAN_baseline := Q_ref$Q_MEAN_baseline]
Q_Summary <- Q_Summary[, Q_MEAN_baseline_SD := Q_ref$Q_MEAN_baseline_SD]
#calculate anomalies 
Q_Summary <- Q_Summary[, Q_MAX_Anomaly := (Q_MAX_5yr-Q_MAX_baseline)/Q_MAX_baseline_SD] 
Q_Summary <- Q_Summary[, Q_MEAN_Anomaly := (Q_MEAN_5yr-Q_MEAN_baseline)/Q_MEAN_baseline_SD] 


## SSC Data from Evan Dethier - Prodduced December 2023
SSC_data <- fread(paste0(wd_root,"SSC/mackenzie_river_ssc_warm_months_summary.csv"))

#mean SSC over 5 year windows, like all other variables
SSC_5_roll <- SSC_data[, list(SSC_5yr = frollmean(SSC_mgL, 5, na.rm =T)),
                        by = site_no #site_no is same as Thiessen_ID but with "st_" prefix 
                        ]

# order both tables by year and Thiessen ID 
SSC_data <- SSC_data[order(site_no, year),] 
SSC_5_roll<- SSC_5_roll[order(site_no)]
#join year column from NDVI_all to NDVI_5_roll 
SSC_5_roll <-cbind(SSC_5_roll,SSC_data$year)
#remove NA rows: all those before 1990
SSC_5_roll <- na.omit(SSC_5_roll)
#actually have 1989 data, so remove that year
SSC_5_roll <- SSC_5_roll[V2>1989,]
SSC_5_roll <- setnames(SSC_5_roll, 'V2', 'Year')
#create baseline data
# Filter out for 1990-2006 baseline calculation
SSC_90_06 <- SSC_5_roll[Year<=2006,]
# Calculate mean/sd from 1975-1995 
# remember that these metrics are for 5year periods, so to get rough idea of annual avg, divide by 5
Baseline_SSC <- SSC_90_06[, list(mean_SSC_90_06 = mean(SSC_5yr),
                                   sd_SSC_90_06 = sd(SSC_5yr)),
                            by = site_no]

# Add mean/sd for each segment buffer into the master table
SSC_5_roll <-SSC_5_roll[Baseline_SSC,on=.(site_no = site_no)]
#calc anomaly
SSC_5_roll <- SSC_5_roll[, SSC_Anomaly := (SSC_5yr - mean_SSC_90_06)/sd_SSC_90_06]
hist(SSC_5_roll$SSC_Anomaly) #just making sure numbers look reasonable. They do. 

SSC_Summary <- SSC_5_roll[, list(mean_SSC_Anomaly = mean(SSC_Anomaly), 
                                   mean_SSC = mean(SSC_5yr),
                                   se_SSC_Anomaly = sd(SSC_Anomaly)/sqrt(length(SSC_Anomaly)),
                                   se_SSC = sd(SSC_5yr)/sqrt(length(SSC_5yr)),
                                   median_SSC_Anomaly = median(SSC_Anomaly)
                                  ),
                            by = Year
                          ]

                 
# Read in Permafrost data
wd_Permafrost <- paste0(wd_root, "Permafrost/")
setwd(wd_Permafrost)
Permafrost_Data <- fread(paste0(wd_Permafrost, "Thiessen_Permafrost_Data_0926.csv")) 

# Clean up the tables
## Join the Permafrost data to the areaChanged data 
#Keep only relevant columns
PF <- Permafrost_Data[,c('point_ID', 'NUM_CODE', 'COMBO', 'RELICT', 'EXTENT', 'CONTENT', 'LANDFORM')]
PF <- setnames(PF, c('NUM_CODE', 'COMBO', 'RELICT', 'EXTENT', 'CONTENT', 'LANDFORM'), c('PF_NUM_CODE', 'PF_COMBO', 'PF_RELICT', 'PF_EXTENT', 'ICE_CONTENT', 'PF_LANDFORM'))
areaChanged_1 <- AreaChange_import[PF, on=.(Thiessen_ID = point_ID)]

###TIDYING: 
#Get rid of upstream area and upstream pix because they are whack. Also the .geo column and the system index column b/c not needed
Water <- data.table::copy(TotalArea_import)
Water <- Water[, c('Upstream_Area_km2', 'Upstream_pix', 'system:index') := NULL]
Water <- setcolorder(Water, c('Thiessen_ID','start_year', 'latitude', 'longitude', 'Polygon_Area_km2', 
                           'Stream_Order', 'Width', 'Elv', 'Hydro_Adj_Elv', 'Water_Area_km2'))

Water <- Water[, join_year := start_year + 5] # add this column for joining Water and Change tables later

Change <- data.table::copy(areaChanged_1) #dublicate table
Change <- Change[, c('Upstream_Area_km2', 'Upstream_pix', 'system:index') := NULL] #delete unnecessary columns
Change <- setcolorder(Change, c('Thiessen_ID','start_year', 'latitude', 'longitude', 'Polygon_Area_km2',
                                'Stream_Order', 'Width', 'Elv', 'Hydro_Adj_Elv', 'PF_COMBO', 'PF_RELICT', 
                                'PF_EXTENT', 'ICE_CONTENT', 'PF_LANDFORM', 'area_gained_km2', 'area_lost_km2')) #reorder 
colnames(Change)[2] <- 'End_Year' # rename column so that it makes sense with respect to period we are measuring erosion over 

### CHECK FOR ERRORS: from Earth Engine (i.e. negative land area)
# Calculate land area in each Thiessen Polygon for each year
Water <- Water[, ":=" (Land_Area_km2 = Polygon_Area_km2 - Water_Area_km2)] #Looks good 10/25/23, no negatives, smallest "land" value is about 2km^2
fwrite(Water, paste0(wd_exports, "Total_Area_Oct25_23.csv"))

# JOIN total area values into area lost table; 
## What we want here is a table of length # of Thiessens x 34 years (1990-2023, inclusive) that, for each year, has an entry for 
## water/land area at the start year (e.g. 1985), water/land area at the end (e.g. 1990), and, of course, area lost/gained in each year.
## These data are necessary for the discharge corrections. 

#create a table that stores water in the start year
Water1 <- Water[, c('Thiessen_ID', 'start_year', 'Water_Area_km2', 'Land_Area_km2')]
colnames(Water1)[3] <- 'Water_Area_Start_km2'
colnames(Water1)[4] <- 'Land_Area_Start_km2'
Water1 <- Water1[, join_year := start_year + 5]
Water1 <- Water1[join_year<=2023,] #get rid of water area for years 2016-2020 bc no start years after that; yes it should be 2023 here
# due to the offset in 'join_year' 

#Master_Area <- Master1[WaterArea_dt, on=.(Thiessen_ID = Thiessen_ID)]


#create a table that stores water in the end/current year
Water2 <- Water[, c('Thiessen_ID', 'start_year', 'Water_Area_km2', 'Land_Area_km2')]
colnames(Water2)[3] <- 'Water_Area_End_km2'
colnames(Water2)[4] <- 'Land_Area_End_km2'
colnames(Water2)[2] <- 'Year'
Water2 <- Water2[Year>=1990,] #there are no 'end' years before 1990, so get rid of those entries 

#PERFORM JOINS
Master1 <- Change[Water1, on=.(Thiessen_ID = Thiessen_ID, End_Year = join_year)]
Master1 <- Master1[Water2, on=.(Thiessen_ID = Thiessen_ID, End_Year = Year)]

#NOT USING CURRENTLY: this makes a table with NA values for area lost/gained first 5 years
#Master <- Change[Water, on=.(Thiessen_ID = Thiessen_ID, End_Year = start_year, latitude = latitude, 
                             # longitude = longitude, Polygon_Area_km2 = Polygon_Area_km2, Stream_Order = Stream_Order, 
                             # Width = Width, Elv = Elv, Hydro_Adj_Elv = Hydro_Adj_Elv ), nomatch = NA] #join


################   END FILE AGGREGATION #############   



########### DISCHARGE CORRECTIONS #############

# Begin Section To Adjust Area Calculations For Higher Flows 
# Logic for this is outlined in Fields et al., Online Methods Figure 2

#correction for scenario in which increase in Q causes apparent erosion (scenario 3 in figure)
# scenario 4 also accounted for here
Master1 <- Master1[, Real_Area_Lost_km2 := ifelse(Water_Area_Start_km2 == Water_Area_End_km2, area_lost_km2, 
                                                  ifelse(Water_Area_Start_km2 < Water_Area_End_km2, 
                                                         area_lost_km2-(Water_Area_End_km2-Water_Area_Start_km2), area_lost_km2))]
#correct for scenario in which change in water area is greater than recorded area lost and we might get negative area lost 
# basically, enforce 0 as lower bound
Master1 <- Master1[, Real_Area_Lost_km2 := ifelse(Real_Area_Lost_km2 < 0, 0, Real_Area_Lost_km2)]

#correction for scenario in which deecrease in Q causes apparent land gained (scenario 5 in figure)
# scenario 6 also accounted for here
Master1 <- Master1[, Real_Area_Gained_km2 := ifelse(Water_Area_Start_km2 == Water_Area_End_km2, area_gained_km2, 
                                                    ifelse(Water_Area_Start_km2 > Water_Area_End_km2, 
                                                           area_gained_km2-(Water_Area_Start_km2-Water_Area_End_km2), area_gained_km2))]
#enforce 0 as lower bound
Master1 <- Master1[, Real_Area_Gained_km2 := ifelse(Real_Area_Gained_km2 < 0, 0, Real_Area_Gained_km2)]

### NORMALIZATION: 
Master1 <- Master1[, Normalized_Area_Lost := Real_Area_Lost_km2/Water_Area_Start_km2]
Master1 <- Master1[, Normalized_Area_Gained := Real_Area_Gained_km2/Water_Area_Start_km2] #using water area here bc I want to normalize 
                                                                                                              #relative to river size
#deal with the few (343) cases of division by zero by just setting those NAs to 0
Master1[, Normalized_Area_Lost := ifelse(is.na(Normalized_Area_Lost), 0, Normalized_Area_Lost)]
Master1[, Normalized_Area_Lost := ifelse(is.infinite(Normalized_Area_Lost), 1.0, Normalized_Area_Lost)]
Master1[, Normalized_Area_Gained := ifelse(is.na(Normalized_Area_Gained), 0, Normalized_Area_Gained)]

WaterArea <- data.table(unique(Master1$Thiessen_ID))

### TEMPERATURE ###

#Above0 <- fread(paste0(wd_imports, 'Above0_Sep23_23_ALL_YEARS.csv'))
#slim down the table to just what's needed, i.e. don't need all the Thiessen variables
Above0 <- Temperature_import[, list(Thiessen_ID, year, latitude, longitude, degreeDays, above0Days)]

# Filter out for 1979-1995 baseline calculation
temp_79_95 <- Above0[year<=1995,]
# Calculate mean/sd from 1975-1995 
# remember that these metrics are for 5year periods, so to get rough idea of annual avg, divide by 5
Baseline_temp <- temp_79_95[, list(mean_degreeDays_79_95 = mean(degreeDays), 
                            mean_above0_79_95 = mean(above0Days),
                            sd_above0_79_95 = sd(above0Days),
                            sd_degreeDays_79_95 = sd(degreeDays)), 
                           by = Thiessen_ID]

# Add mean/sd for each segment buffer into the master table
Master1 <- Master1[Baseline_temp,on=.(Thiessen_ID = Thiessen_ID)] 

# Filter out for period of area change record (1990-2020)
temp_studyPeriod <- Above0[year>=1990,]
temp_studyPeriod <- temp_studyPeriod[, list(Thiessen_ID, year, degreeDays, above0Days)] #get rid of redundant variables that already exist in Master table
# add these data to the master area change table
Master1 <- Master1[temp_studyPeriod, on=.(Thiessen_ID = Thiessen_ID, End_Year = year)]


## Generate latitude breaks
Master1[, north_south := cut(latitude, 
                              breaks = c(min(Master1$latitude-1, na.rm = T), 
                                          53,54,55, 56, 57, 58, 59, 60, 61,
                                          62, 63, 64, 65, 66, 67, 68, 
                                          max(Master1$latitude, na.rm = T), include.lowest = T))]
Master1[, N_S_2deg := cut(latitude, 
                          breaks = c(
                                     52,54,56,58,60, 
                                     62,64,66,68,
                                     Inf),
                          labels = rev(c("(68,70]","(66,68]", "(64,66]", "(62,64]", "(60,62]",
                                     "(58,60]", "(56,58]", "(54,56]", "(52,54]")))]
Master1 <- Master1[, Percent_Area_Lost := (Normalized_Area_Lost*100)]

### Normalize Percent Area Lost (Each Lat Bin, Stream Order, and year) by the pre-2006 average of Percent Area Lost
Perc_Lost_Pre2006 <- Master1[End_Year<=2006,]
# Calculate mean/sd percent area lost from 1990-2006 for each Thiessen
Baseline_Erosion <- Perc_Lost_Pre2006[, list(mean_Area_Lost_1990_2006 = mean(Percent_Area_Lost),
                                             SD_Area_Lost_1990_2006 = sd(Percent_Area_Lost)), 
                                      by = .(Thiessen_ID)] #Split up by order and latitude later
Baseline_Erosion[, SD_Area_Lost_1990_2006 := ifelse(SD_Area_Lost_1990_2006 == 0, 0.01, SD_Area_Lost_1990_2006)] # doing this to avoid division by zero in a few cases

# Join this table to the Master1 Table 
Master1 <- Master1[Baseline_Erosion, on=.(Thiessen_ID = Thiessen_ID)] #, N_S_2deg = N_S_2deg, Stream_Order = Stream_Order
Master1 <- Master1[Thiessen_ID != 1179839] # delete a problem thiessen - this one produces wild anomalies because it is a tiny stream. 

#Master_test <- Master1[, Perc_Lost_Anomaly := (Percent_Area_Lost - mean_Area_Lost_1990_2006) / SD_Area_Lost_1990_2006]

Master1 <- Master1[, PF_EXTENT := ifelse(PF_EXTENT == "", "None", 
                                         ifelse(PF_EXTENT == "C", "Continuous", 
                                                ifelse(PF_EXTENT == "D", "Discontinuous", 
                                                       ifelse(PF_EXTENT == "S", "Sporadic", 
                                                              ifelse(PF_EXTENT == "I", "Isolated", PF_EXTENT)
                                                              )
                                                       )
                                                )
                                         )
                   ]

#join to NDVI data master 
Master1 <- NDVI_5_roll[Master1,on=.(Thiessen_ID = Thiessen_ID, Year = End_Year)]

fwrite(Master1, paste0(wd_exports, "Master1_ArcticErosion.csv"))





############## GROUNDTRUTH DATA #####################
#Groundtruth site selection:
# sites: 2163953, 578989, 4644595, 5610933, 2704674, 746629, 5241127, 5989641, 8489055, 8840305
# listed in order: 1-10; The highest latitude site (site 1; Thiessen 2163953) was ultimately excluded because we excluded the delta region data from
  # our final analysis

GT_sel <- Master1[Thiessen_ID == 8840305, ]
#Groundtruth discharge correction:
setwd(wd_sites)
GT_data<- fread("Groundtruth_Data.csv")

### DISCHARGE CORRECTIONS: 
# Begin Section To Adjust Area Calculations For Higher Flows 
# Logic for this is outlined in a very nice figure

#correction for scenario in which increase in Q causes apparent erosion (scenario 3 in figure)
# scenario 4 also accounted for here
GT_data <- GT_data[, GT_Real_Area_Lost_km2 := ifelse(GT_Start_Area_sq_km == GT_End_Area_sq_km, GT_Area_Lost_sq_km, 
                                                  ifelse(GT_Start_Area_sq_km < GT_End_Area_sq_km, 
                                                         GT_Area_Lost_sq_km-(GT_End_Area_sq_km-GT_Start_Area_sq_km), GT_Area_Lost_sq_km))]

#correct for scenario in which change in water area is greater than recorded area lost and we might get negative area lost 
# basically, enforce 0 as lower bound
GT_data <- GT_data[, GT_Real_Area_Lost_km2 := ifelse(GT_Real_Area_Lost_km2 < 0, 0, GT_Real_Area_Lost_km2)]
## not accounting for land gained, because doesn't matter in this analysis

### NORMALIZATION: 
GT_data <- GT_data[, GT_Perc_Area_Lost := (GT_Real_Area_Lost_km2/GT_Start_Area_sq_km)*100]
GT_data <- GT_data[, GEE_Perc_Area_Lost := (REAL_Area_Lost_GEE_sq_km/GEE_Start_Area_sq_km)*100]
GT_data <- GT_data[Site != 1,]
GT_data <- GT_data[, Site_new := seq(1,9,1)]

#mean % error in water detection
GT_mean_err <- GT_data[, list(mean_start_err = mean(((GT_Start_Area_sq_km - GEE_Start_Area_sq_km)/GT_Start_Area_sq_km)*100), 
                              mean_end_err = mean(((GT_End_Area_sq_km - GEE_End_Area_sq_km)/GT_End_Area_sq_km)*100)
                              )
                              ]

## Plots of Groundtruth data 
GT_plot <- ggplot(GT_data, aes(x=Site, y = GT_Perc_Area_Lost )) +
  geom_point(size = 5, pch = 21, stroke = 2, color = "#8480F2") +
  geom_point(aes(x=Site, y = GEE_Perc_Area_Lost), size = 5, pch = 3, stroke = 2, color = "#F2A03D") +
  geom_text(aes(x = Site, y = GT_Perc_Area_Lost
                , label = paste0(round(GT_Perc_Area_Lost,2),"%")
                , hjust = -0.25
                
  ), size = 5, fontface = "bold") +
  #geom_text(aes(label = round(GEE_Perc_Area_Lost, digits = 2)), nudge_x = 0.5, nudge_y = -2) +
  geom_text(aes(x = Site, y = GEE_Perc_Area_Lost
                , label = paste0(round(GEE_Perc_Area_Lost,2),"%")
                , hjust = -0.25
                
  ), size = 5, fontface = "bold") +
  #ylim(0, 30) +
  geom_errorbar(aes(ymin = (GEE_Perc_Area_Lost), ymax = (GT_Perc_Area_Lost)), color = 'black', width = 0.2, alpha = 0.5) +
  theme_pubr() +
  theme(legend.position = 'right',
        axis.text = element_text(size=16, face = 'bold')) +
  ggtitle("Hand-Delineated vs. Algorithmic Area Lost \n") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  labs(
    x = '\n Site',
    y = 'Percent Area Lost \n'
  ) +
  font("ylab", size = 18, face = "bold") +
  font("xlab", size =18, face = "bold") 

print(GT_plot)

GT_water_plot <- ggplot(GT_data, aes(x=Site_new, y = GT_Start_Area_sq_km, label = paste0(round(GEE_Start_Area_sq_km,2)," km2"))) +
  geom_point(size = 5, pch = 21, stroke = 2, color = "#033F73") +
  geom_text(aes(x = Site_new, y =  GT_Start_Area_sq_km
                , label = paste0(round(GT_Start_Area_sq_km,2)," km2")
                , hjust = -0.25
                
  ), size = 5, fontface = "bold") +
  #geom_text_repel(size = 5, fontface = "bold", nudge_x = -0.5)+
  geom_point(aes(x=Site_new, y = GEE_Start_Area_sq_km), size = 5, pch = 3, stroke = 2, color = "#0597F2") +
  geom_text(aes(x = Site_new, y = GEE_Start_Area_sq_km
                , label = paste0(round(GEE_Start_Area_sq_km,2)," km2")
                , hjust = -0.25
                
  ), size = 5, fontface = "bold") +
  #geom_text_repel(label = paste0(round(GEE_Start_Area_sq_km,2),"km2"), size = 5, fontface = "bold")+
  #geom_text(aes(label = round(GEE_Perc_Area_Lost, digits = 2)), nudge_x = 0.5, nudge_y = -2) +
  #ylim(0, 30) +
  geom_errorbar(aes(ymin = (GEE_Start_Area_sq_km), ymax = (GT_Start_Area_sq_km )), color = 'black', width = 0.2, alpha = 0.5) +
  scale_x_continuous(limits = c(1, 9), breaks = seq(1,9,1)) +
  theme_pubr() +
  theme(legend.position = 'right',
        axis.text = element_text(size=16, face = 'bold')) +
  ggtitle("Hand-Delineated vs. Algorithmic Area \n") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  labs(
    x = '\n Site',
    y = 'Area (km2)\n'
  ) +
  font("ylab", size = 18, face = "bold") +
  font("xlab", size =18, face = "bold") 

print(GT_water_plot)

GT_end_water_plot <- ggplot(GT_data, aes(x=Site_new, y = GT_End_Area_sq_km)) +
  geom_point(size = 5, pch = 21, stroke = 2, color = "#033F73") +
  geom_text(aes(x = Site_new, y =  GT_End_Area_sq_km
                , label = paste0(round(GT_End_Area_sq_km,2)," km2")
                , hjust = -0.25
                
  ), size = 5, fontface = "bold") +
  #geom_text_repel(size = 5, fontface = "bold", nudge_x = -0.5)+
  geom_point(aes(x=Site_new, y = GEE_End_Area_sq_km), size = 5, pch = 3, stroke = 2, color = "#0597F2") +
  geom_text(aes(x = Site_new, y = GEE_End_Area_sq_km
                , label = paste0(round(GEE_End_Area_sq_km,2)," km2")
                , hjust = -0.25
                
  ), size = 5, fontface = "bold") +
  #geom_text_repel(label = paste0(round(GEE_Start_Area_sq_km,2),"km2"), size = 5, fontface = "bold")+
  #geom_text(aes(label = round(GEE_Perc_Area_Lost, digits = 2)), nudge_x = 0.5, nudge_y = -2) +
  scale_x_continuous(limits = c(1, 9), breaks = seq(1,9,1)) +
  geom_errorbar(aes(ymin = (GEE_End_Area_sq_km), ymax = (GT_End_Area_sq_km )), color = 'black', width = 0.2, alpha = 0.5) +
  theme_pubr() +
  theme(legend.position = 'right',
        axis.text = element_text(size=16, face = 'bold')) +
  ggtitle("Hand-Delineated vs. Algorithmic Area \n") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  labs(
    x = '\n Site',
    y = 'Area (km2)\n'
  ) +
  font("ylab", size = 18, face = "bold") +
  font("xlab", size =18, face = "bold") 

print(GT_end_water_plot)




############ EXCLUDE THE NORTHERNMOST (Mackenzie River Delta) STUDY SITES #############

#Look at results without the highest latitude bin, which is affected by lakes
Master1 <- Master1[N_S_2deg != '(68,70]',]
setnames(Master1, 'Year', 'End_Year')

Master_Thiessens <- length(unique(Master1$Thiessen_ID)) #seeing how many Thiessens are left after taking out the northernmost points

# Only use lines below to generate first figure
Master1 <- Master1[, Perc_Lost_Anomaly := 
                            (Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006] 

Master1 <- Master1[, Perc_Lost_Anomaly_SE := 
                            sd((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006)/sqrt(length(Percent_Area_Lost))]


############ Simple model of thaw depth to see if that predicts erosion -- not used in the end #############

# Mean Thaw Deoth (Active Layer Calculation)
## Using Farquharson et al.,  2019 Eq. 2 and getting A from the slope of the line of best fit for Degree Days vs. their observed thaw depths
### I digitized their data from Figure 2a and 2b and then did the regression myself, the equation of that line is: y = 0.1033x + 21.306; R2 = 0.78
Master1 <- Master1[, Mean_Thaw_Depth := sqrt(degreeDays/5)*0.1] #use slope of 0.018 if taken from multiple regression in their Table S5
plot(Master1$degreeDays/5, Master1$Mean_Thaw_Depth)
#create bins 
Master1[, MTD_bins := cut(Mean_Thaw_Depth,
                          breaks = c(
                            -Inf,3.0,3.5,4.0,Inf),
                          labels = c("(2.5,3.0]", "(3.0,3.5]", "(3.5,4.0]", "(4.0,4.5]"))]
# Master1[, MTD_bins := cut(Mean_Thaw_Depth, 
#                           breaks = c(
#                             0.43,0.50,0.60,0.70,0.76),
#                           labels = c("(0.43,0.50]","(0.50,0.60]", "(0.60,0.70]", "(0.70,0.76]"))]






################ SUMMARIZE DATA IN MASTER1 by YEAR, LATITUDE BIN, STREAM ORDER, PERMAFROST EXTENT, ETC ###############
  ## These are the tables used to make most of the plots presented in Fields et al.

#Summarize by year only 
Summary1 <- Master1[, list(Perc_Lost_Anomaly = 
                             mean((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006), 
                          Perc_Lost_Anomaly_SE = 
                             (sd((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006))/sqrt(length(Percent_Area_Lost))
                          #meanAnnual_NDVI_ratio = 
                            # mean(NDVI_Ratio)
                          ),
                    by = c('End_Year')] 

#sumamrize by stream order for each year
Summary2 <- Master1[, list(Perc_Lost_Anomaly = 
                             mean((Percent_Area_Lost - mean_Area_Lost_1990_2006) / SD_Area_Lost_1990_2006),
                           Perc_Lost_Anomaly_SE = 
                             (sd((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006))/sqrt(length(Percent_Area_Lost))
                          ),
                             by = .(End_Year, Stream_Order)
                    ]

#summarize by latitude bin for each year 
Summary3 <- Master1[, list(Perc_Lost_Anomaly = 
                                         mean((Percent_Area_Lost - mean_Area_Lost_1990_2006) / SD_Area_Lost_1990_2006),
                           Perc_Lost_Anomaly_SE = 
                             (sd((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006))/sqrt(length(Percent_Area_Lost))
                           ), 
                                by = .(End_Year, N_S_2deg)
                    ]

#summarize by Permafrost Extent\ for each year
Summary4 <- Master1[, list(Perc_Lost_Anomaly = 
                             mean((Percent_Area_Lost - mean_Area_Lost_1990_2006) / SD_Area_Lost_1990_2006),
                           Perc_Lost_Anomaly_SE = 
                             (sd((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006))/sqrt(length(Percent_Area_Lost))
                           ), 
                             by = .(End_Year, PF_EXTENT)
                    ]

#summarize by mean thaw depth bin for each year 
Summary5 <- Master1[, list(Perc_Lost_Anomaly = 
                             mean((Percent_Area_Lost - mean_Area_Lost_1990_2006) / SD_Area_Lost_1990_2006),
                           Perc_Lost_Anomaly_SE = 
                             (sd((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006))/sqrt(length(Percent_Area_Lost))
                            ), 
                              by = .(End_Year, MTD_bins)
                    ]


#add the discharge data to the summary table
Summary1 <- Q_Summary[Summary1, on=.(Year = End_Year)]
Summary2 <- Q_Summary[Summary2, on=.(Year = End_Year)]
Summary3 <- Q_Summary[Summary3, on=.(Year = End_Year)]
Summary4 <- Q_Summary[Summary4, on=.(Year = End_Year)]
Summary5 <- Q_Summary[Summary5, on=.(Year = End_Year)]

#add the NDVI data to the summary table
Summary1 <- NDVI_Summary[Summary1, on=.(Year = Year)]
Summary2 <- NDVI_Summary[Summary2, on=.(Year = Year)]
Summary3 <- NDVI_Summary[Summary3, on=.(Year = Year)]
Summary4 <- NDVI_Summary[Summary4, on=.(Year = Year)]
Summary5 <- NDVI_Summary[Summary5, on=.(Year = Year)]

#add the NDVI data to the summary table
Summary1 <- SSC_Summary[Summary1, on=.(Year = Year)]
Summary2 <- SSC_Summary[Summary2, on=.(Year = Year)]
Summary3 <- SSC_Summary[Summary3, on=.(Year = Year)]
Summary4 <- SSC_Summary[Summary4, on=.(Year = Year)]
Summary5 <- SSC_Summary[Summary5, on=.(Year = Year)]


#Create a Temperature Summary Table - average in each year
  #Calculate anomaly from 1979-1995 averages

Master1[, DegDay_Anomaly := (degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95] 

Master1[, Above0_Anomaly := (above0Days-mean_above0_79_95)/sd_above0_79_95]


#average in each year
DegDay_Summary <- Master1[, list(Annual_DegD_Anomaly = mean((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95),
                                 Annual_DegD_SD = sd(degreeDays),
                                 Annual_DegD_Mean = mean(degreeDays),
                                 Annual_DegD_Anomaly_SE = sd((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95)/sqrt(length(degreeDays))
                                 ), 
                            by = .(End_Year)
                          ]

#average in each year for each STREAM ORDER
Temp_Summary <- data.table(Master1[, list(End_Year, Stream_Order, N_S_2deg, Thiessen_ID, DegDay_Anomaly, Above0_Anomaly)])
DegDay_Summary1 <- Master1[, list(Annual_DegD_Anomaly = mean(DegDay_Anomaly),
                                 Annual_DegD_SD = sd(degreeDays),
                                 Annual_DegD_Mean = mean(degreeDays),
                                 Annual_DegD_Anomaly_SE = sd((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95)/sqrt(length(degreeDays))
                                 ), 
                              by = .(End_Year, Stream_Order)
                          ]
#average in each year for each LATITUDE BIN
DegDay_Summary2 <- Master1[, list(Annual_DegD_Anomaly = mean(DegDay_Anomaly),
                                  Annual_DegD_SD = sd(degreeDays),
                                  Annual_DegD_Mean = mean(degreeDays),
                                  Annual_DegD_Anomaly_SE = sd((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95)/sqrt(length(degreeDays))
                                  ), 
                              by = .(End_Year, N_S_2deg)
                            ]

#average in each year for each PERMAFROST EXTENT 
DegDay_Summary3 <- Master1[, list(Annual_DegD_Anomaly = mean(DegDay_Anomaly),
                                  Annual_DegD_SD = sd(degreeDays),
                                  Annual_DegD_Mean = mean(degreeDays),
                                  Annual_DegD_Anomaly_SE = sd((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95)/sqrt(length(degreeDays))
                                  ), 
                             by = .(End_Year, PF_EXTENT)
                            ]

#average in each year for each MEAN THAW DEPTH BIN
DegDay_Summary4 <- Master1[, list(Annual_DegD_Anomaly = mean(DegDay_Anomaly),
                                  Annual_DegD_SD = sd(degreeDays),
                                  Annual_DegD_Mean = mean(degreeDays),
                                  Annual_DegD_Anomaly_SE = sd((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95)/sqrt(length(degreeDays))
                                  ), 
                              by = .(End_Year, MTD_bins)
                            ]

#Add Degree Day Summary Tables to Larger Summary Tables
Summary1 <- DegDay_Summary[Summary1, on=.(End_Year = Year)]
Summary2 <- DegDay_Summary1[Summary2, on=.(End_Year = Year, Stream_Order = Stream_Order)]
Summary3 <- DegDay_Summary2[Summary3, on=.(End_Year = Year, N_S_2deg = N_S_2deg)]
Summary4 <- DegDay_Summary3[Summary4, on=.(End_Year = Year, PF_EXTENT = PF_EXTENT)]
Summary5 <- DegDay_Summary4[Summary5, on=.(End_Year = Year, MTD_bins = MTD_bins)]

#Save all the summary tables 
fwrite(Summary1, paste0(wd_exports, "Summary1_ArcticErosion_Year.csv")) # table summarized by END_YEAR
fwrite(Summary2, paste0(wd_exports, "Summary2_ArcticErosion_StrmOrd.csv")) # table summarized by STREAM ORDER
fwrite(Summary3, paste0(wd_exports, "Summary3_ArcticErosion_LatBins.csv")) # table summarized by LATITUDE BIN
fwrite(Summary4, paste0(wd_exports, "Summary4_ArcticErosion_Permafrost_Extent.csv")) # table summarized by PERMAFROST EXTENT 
fwrite(Summary5, paste0(wd_exports, "Summary5_ArcticErosion_Thaw_Depth.csv")) # table by summarized MEAN THAW DEPTH



###### FUTURE TEMP/EROSION PROJECTIONS ###########
Future <- FutureTemp_data_files[Baseline_temp, on=.(Thiessen_ID = Thiessen_ID)] # Future temperature -- using CanESM5 only, change to CMIP6_Import if you want that
#clean up the table, keeping only columns we want
Future <- Future[, .(Thiessen_ID, year, above0Days, degreeDays, mean_degreeDays_79_95, mean_above0_79_95, sd_degreeDays_79_95, sd_above0_79_95)] 


Future_Summary <- Future[, list(F_Annual_DegD_Anomaly = mean((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95),
                                 F_Annual_DegD_SD = sd(degreeDays),
                                 F_Annual_DegD_Mean = mean(degreeDays),
                                 F_Annual_DegD_Anomaly_SE = sd((degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95)/sqrt(length(degreeDays))
                                ), 
                            by = .(year)
                         ]
#calculate erosion in each year based on regression equation: Annual Erosion Anomaly = -0.051 + 0.12(Annual_Temp_Anomaly); R2 = 0.29; P = 0.00098
Future_Summary <- Future_Summary[, predicted_Eros_Anomaly := -0.051 + 0.12*F_Annual_DegD_Anomaly]

#Add these data to summary1 table
Summary_F1 <- data.table(rbind(Summary1, Future_Summary, fill = T))
Summary_F1 <- Summary_F1[, End_Year := ifelse(is.na(End_Year), year, End_Year)]
        
Summary_F1 <- Summary_F1[, Perc_Area_Lost := ifelse(is.na(predicted_Eros_Anomaly), NA, ((predicted_Eros_Anomaly*3.35)+5.31))] #mean and SD from baseline erosion period 

Summary_F1a <- Summary_F1[, .(End_Year, Perc_Area_Lost)] # select only columns we want to add to big table
setnames(Summary_F1a, 'Perc_Area_Lost', 'mean_Annual_Perc_Lost') # rename before combining 


#use this for estimate of total area eroded  
WaterArea_dt <- Water1[, list(mean_Water_km2 = mean(Water_Area_Start_km2)),
                       by = Thiessen_ID
]

#create a table with all years, all thiessens, each with mean water area 
Yr <- c(1990:2050)

#adding water area to future years - imperfect as this may change in the future, but the best we can do and need to do it to get anomalies
start_table <- WaterArea_dt

for(i in 1:61){
  Year_Select <- Yr[i]
  output <- start_table[, Year := Year_Select]
  WaterArea_dt <- rbind(WaterArea_dt, output)
}

#get rid of duplicated 1990 rows
# Set keys - this sorts the data based on these values
setkeyv(WaterArea_dt, c('Thiessen_ID', 'Year'))
# keep unique observations (I also remove the variable x)
uniqdat <- subset(unique(WaterArea_dt))
WaterArea_dt <- uniqdat
WaterArea_dt <- WaterArea_dt[Summary_F1a, on=.(Year=End_Year)]

#calc area eroded in each thiessen in each year; this is not as good as the 'real_Area_lost_km2' data, which is directly 
# measured but it is the best we can do for the future years. 
WaterArea_dt <- WaterArea_dt[!is.na(mean_Annual_Perc_Lost), Annual_Km_Eroded := mean_Water_km2*(mean_Annual_Perc_Lost/100)] #ignore but retain NAs

#for the real data:
Eroded_Area_REAL <- Master1[, list(Km2_eroded = sum(Real_Area_Lost_km2)),
                                  by = .(End_Year)
                            ]

Eroded_Area <- WaterArea_dt[, list(Km2_eroded = sum(Annual_Km_Eroded)),
                             by = .(Year)
                           ]

Eroded_Area <- Eroded_Area[, Km2_eroded := ifelse(is.na(Km2_eroded), Eroded_Area_REAL$Km2_eroded, Km2_eroded)] #add in the real data

#Quick sanity check: recall that our discharge corrections are so conservative that the area eroded in km2 should not be taken as anything close to 
  #the truth, it is an extreme underestimate. Just making sure the future erosion is approximatley in line with the real data from the past. 
ErodedArea_REAL_gg <- ggplot(Eroded_Area, aes(x=Year, y=Km2_eroded)) +
  geom_line(size = 2, color = "red", linetype= 'dotted', alpha = 0.8) +
  geom_line(data = Eroded_Area_REAL, aes(x=End_Year, y=Km2_eroded), size = 2, color = "#bd925a")+
  theme_pubr

print(ErodedArea_REAL_gg)




###############################################
##                PLOTTTING                  ##
###############################################

#############################################################################
##                      Normalized Area Loss Plots                         ##
#############################################################################

# Z-SCORES: Start by calculating the mean and std dev for each latitude and stream order combination 
## For normalized Data
Annual_mean <- Master1[, list(mean_Norm_Area_Lost = mean(Percent_Area_Lost, na.rm = TRUE)), by=c('Stream_Order', 'N_S_2deg')]
Annual_sd <- Master1[, list(sd_Norm_Area_Lost = sd(Percent_Area_Lost, na.rm = TRUE)), by=c('Stream_Order', 'N_S_2deg')]
# Now join these parameters to main table in order to calc z-scores: dt_a[dt_b, on = .(id = id, date = date), roll = TRUE]
Annual_Stats <- Annual_mean[Annual_sd, on = .(Stream_Order = Stream_Order, N_S_2deg = N_S_2deg)]

## For Raw area lost data
Annual_mean_Raw <- Master1[, list(mean_Area_Lost = mean(Real_Area_Lost_km2, na.rm = TRUE)), by=c('Stream_Order', 'N_S_2deg')]
Annual_sd_Raw <- Master1[, list(sd_Area_Lost = sd(Real_Area_Lost_km2, na.rm = TRUE)), by=c('Stream_Order', 'N_S_2deg')]
# Now join these parameters to main table in order to calc z-scores: dt_a[dt_b, on = .(id = id, date = date), roll = TRUE]
Annual_Stats_Raw <- Annual_mean_Raw[Annual_sd_Raw, on = .(Stream_Order = Stream_Order, N_S_2deg = N_S_2deg)]

# Add mean and sd for each Stream Order and Latitude Bin to the corresponding rows in Master table
Master1 <- merge(Master1, Annual_Stats, by = c("Stream_Order", "N_S_2deg"), all.x = TRUE)
Master1 <- merge(Master1, Annual_Stats_Raw, by = c("Stream_Order", "N_S_2deg"), all.x = TRUE)

# Calculate the z-score of percent area lost
Master1[,Z_Score_Percent := (Percent_Area_Lost - mean_Norm_Area_Lost)/sd_Norm_Area_Lost]
#deal with the few (343) cases of division by zero by just setting those NAs to 0
Master1[, Z_Score_Percent := ifelse(is.na(Z_Score_Percent), 0, Z_Score_Percent)] # is setting to zero fair? If not, what should I do? Somehow use NAs in lm(); Ask Evan. 

# Calculate the z-score of RAW area lost (not normalized)
Master1[,Z_Score_AreaLost := (Real_Area_Lost_km2 - mean_Area_Lost)/sd_Area_Lost]
#deal with the few (343) cases of division by zero by just setting those NAs to 0
Master1[, Z_Score_AreaLost := ifelse(is.na(Z_Score_AreaLost), 0, Z_Score_AreaLost)] # is setting to zero fair? If not, what should I do? Somehow use NAs in lm(); Ask Evan. 
#Master1[, Z_Score := ifelse(is.na(Normalized_Area_Gained), 0, Normalized_Area_Gained)]



# Reordering group factor levels so I can plot North to South 
Master1$N_S_2deg <- factor(Master1$N_S_2deg,               
                                    levels = c("(66,68]", "(64,66]", "(62,64]", "(60,62]",
                                               "(58,60]", "(56,58]", "(54,56]", "(52,54]")) #"(68,70]",



####
## Now, for the master plot of Normalized Area Loss over time by river order and latitude ##
###
Strm_Ord <- c(7:11)
Lat_Bins <- c("(66,68]", "(64,66]", "(62,64]", "(60,62]",
              "(58,60]", "(56,58]", "(54,56]", "(52,54]") #"(68,70]"

#IMPORTANT: Here is how to extract the stats for slopes and p-values of each stream order and latitude bin
for(i in 1:5){
    for(j in 1:9){
      Lat_Bin_select <- Lat_Bins[j]
      Strm_ord_select <- Strm_Ord[i]
      riv_data <- Master1[Stream_Order == Strm_ord_select & N_S_2deg == Lat_Bin_select]
      if(nrow(riv_data)>3){
        mod <- lm(Perc_Lost_Anomaly ~ End_Year, data = riv_data)
        summary_table_select <- data.table(tidy(mod))[2,.(estimate, p.value)]
        summary_table_select[, ":=" (Stream_Order = Strm_ord_select, 
                                     N_S_2deg = Lat_Bin_select,
                                     significance = ifelse(p.value < 0.001, '**',
                                                      ifelse(p.value <0.05, "*", "")))]
        if(i==1 & j==1){
          summary_table <- summary_table_select
        }else{
          summary_table <- rbind(summary_table, summary_table_select, use.names=T)
        }
      }
    }
}


coeff1 = 5 # this is. for scaling secondary axes in later plots
#plot itself -- Note that the regression lines here are still based on all the data, not the annual means, which are plotted as points
## Change input to percent area lost or real area lost, depending on what you want 
BigBoy <- ggplot(Master1, aes(x = End_Year, y = Perc_Lost_Anomaly)) +
  #geom_point(aes(color = Stream_Order), pch = 21, alpha = 0.2) +
  stat_summary(fun = "mean", geom = "point", aes(color = Stream_Order), pch = 21, alpha = 0.6) +
  stat_summary(fun = "mean", geom = "line", aes(color = Stream_Order), color = "#bd925a", linewidth = 2) +
  #geom_line(color = "#bd925a", linewidth = 4) +
  geom_smooth(method = 'lm', lty = 'dashed', linewidth = 1, color = '#0554F2') + 
  #geom_smooth(lty = 'dashed', linewidth = 2, color = '#0554F2')+
  #geom_smooth(aes(x=End_Year, y=DegDay_Anomaly/coeff1), color = "#F21905", linewidth = 2) +
  #facet_wrap(.~north_south + RIV_ORD) +
  #scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5)) +
  facet_grid(rows = vars(N_S_2deg), cols=vars(Stream_Order)) +
  #coord_cartesian(ylim = c(- 0, 1.5)) +
  #coord_cartesian(ylim = c(0.0, 0.6)) +
  scale_y_continuous(
    #name = "Erosion Anomaly \n",
    #sec.axis = sec_axis(~.*coeff, name="Degree Day Anomaly")
  ) +
  theme_pubr +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    strip.text.y = element_text(size = 15, face = "bold"),
    legend.position = "none",
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="bold"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=16)
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 20, face = "bold")) +
  #stat_regline_equation(label.y = .75, size=2) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  # the next line adds the slope and p-values to each facet
  geom_text(data = summary_table, aes(x=1991, y=2.5, label = paste0(round(estimate, 3), '%/yr ', significance)), hjust = 0, size = 5, fontface= 'bold') + 
  labs(
    title = "Percent Area Lost vs. Time",
    subtitle = "by Stream Order and Latitude",
    x = '',
    y = 'Percent Area Lost Anomaly \n' #\n (Land Lost / Total Water at Start) \n'
  ) 

print(BigBoy)

## Scatterplot and regression of Temp Anomaly and Erosion Anomaly
T_Erosion_Regression <- ggplot(Summary1, aes(x = Annual_DegD_Anomaly, y = Perc_Lost_Anomaly)) + # 
  theme_pubr() +
  geom_point(color = "#025949", stroke =2, size = 6, pch =21) +
  geom_smooth(method = 'lm', formula= (y ~ x), lty = 'dashed', lwd = 1, color = 'black') + #formula= (y ~exp(x)),
  xlab("\n Annual Temp Anomaly") +
  ylab("Annual Erosion Anomaly \n")+
  font("ylab", size = 18) +
  font("xlab", size =18) +
  #stat_regline_equation(label.y =20000, aes(label = ..rr.label..)) 
  stat_regline_equation(
    formula = (y ~ x), 
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~`,`~")), #label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    label.x = -1, 
    label.y = 0.5, 
    size = 6
  ) +
  stat_cor(aes(label = paste(..p.label.., sep = "~`,`~")), label.x = 0.2, label.y = 0.4, size = 6)


print(T_Erosion_Regression)

#Look at regression report
mod <- lm(Summary1$Perc_Lost_Anomaly ~ Summary1$Annual_DegD_Anomaly, df = (33/5))
print(summary(mod))
#extract RMSE of regression model: RMSE = 0.1679565
sqrt(mean(mod$residuals^2))

#Multiple Linear Regression
Mult_Reg <- lm(formula = Perc_Lost_Anomaly ~ Annual_DegD_Anomaly + Q_MEAN_Anomaly, data = Summary1)
print(summary(Mult_Reg))




####
## Residuals of PERCENT AREA LOST OVER TIME BY RIVER ORDER AND LATITUDE ##
###

mod <- lm(Percent_Area_Lost ~ End_Year, data = Master1)
Master1 <- Master1[, Residuals := mod$residuals] #343 rows are removed in the regression for "non-finite" values; they are NAs, b/c of division by zero earlier, fixed now

for(i in 1:5){
  for(j in 1:9){
    Lat_Bin_select <- Lat_Bins[j]
    Strm_ord_select <- Strm_Ord[i]
    riv_data <- Master1[Stream_Order == Strm_ord_select & N_S_2deg == Lat_Bin_select]
    if(nrow(riv_data)>3){
      #mod <- lm(Z_Score ~ End_Year, data = riv_data)
      summary_table_select <- data.table(tidy(mod))[2,.(estimate, p.value)]
      summary_table_select[, ":=" (Stream_Order = Strm_ord_select, 
                                   N_S_2deg = Lat_Bin_select,
                                   significance = ifelse(p.value < 0.01, '**',
                                                         ifelse(p.value <0.05, "*", "")))]
      if(i==1 & j==1){
        summary_table <- summary_table_select
      }else{
        summary_table <- rbind(summary_table, summary_table_select, use.names=T)
      }
    }
  }
}

#plot itself -- Note that the regression lines here are still based on all the data, not the annual means, which are plotted as points
Residuals_Plot <- ggplot(Master1, aes(x = End_Year, y = Residuals)) +
  #geom_point(aes(color = Stream_Order), pch = 21, alpha = 0.2) +
  stat_summary(fun = "mean", geom = "point", color = '#2B65D9', pch = 19, alpha = 0.6) +
  #geom_smooth(method = 'lm', lty = 'dashed', linewidth = 1, color = '#0554F2') + 
  #facet_wrap(.~N_S_2deg) +
  facet_grid(rows = vars(N_S_2deg), cols=vars(Stream_Order)) +
  #coord_cartesian(ylim = c(0.0, 0.6)) +
  theme_pubr +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    strip.text.y = element_text(size = 15, face = "bold"),
    legend.position = "none",
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="plain"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=14)
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  #stat_regline_equation(label.y = .75, size=2) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  # the next line adds the slope and p-values to each facet
  #geom_text(data = summary_table_Z, aes(x=1991, y=3, label = paste0(round(estimate, 3), significance)), hjust = 0) + 
  labs(
    title = "Residuals of Percent Area Loss vs. Time",
    subtitle = "by Latitude and Stream Order",
    x = '',
    y = 'Residual (annual mean)'
  ) 

print(Residuals_Plot)


#################################################
#                                               #
###  CORRELATION ANALYSIS & RESIDUAL TESTS   ####
#             By Stream Order                   #
#################################################

## Determine if the temperature data is serially correlated:
x <- DegDay_Summary$End_Year
y <- DegDay_Summary$Annual_DegD_Anomaly
# fit spline 
mod.ss <- ss(x, y, nknots = 10) 
# plot results
plot(x, y) +
lines(x, mod.ss$y, lty = 2, col = 2, lwd = 2)
#model summary to get residuals
mod.sum <- summary(mod.ss)
Temperature_Residuals <- data.table(mod.sum$residuals)
Temperature_Residuals <- Temperature_Residuals[, Year := list(1990:2023)]
SummaryResd <- Summary1[Temperature_Residuals, on =. (End_Year = Year)] # add to Summary table
setnames(SummaryResd, 'V1', 'DegDay_Anomaly_Residuals')
#plot residuals
zeroes <- data.table(matrix(0,34,1))
plot(Temperature_Residuals$Year, Temperature_Residuals$V1) +
  legend("bottomright", legend = "Temperature Anonmaly Residuals", bty = "n") 
  lines(x, zeroes$V1) 
  
x
Temperature_Residuals

## Determine if the erosion data is serially correlated:
#calculate the residuals relative to a spline fit
mod_ssA <- ss(Summary2$End_Year, Summary2$Perc_Lost_Anomaly, nknots = 10) ##ss does the spline fitting with a max of 10 nodes 
mod_sumA <- summary(mod_ssA)
mod_sumA
resdA <- mod_sumA$residuals
resdA
Summary2 <- Summary2[, Perc_Lost_Residuals := mod_sumA$residuals]


#create lag column
Summary2 <- Summary2[, Perc_Lost_Resd_Lags := c(NA, Perc_Lost_Residuals[-.N]), by = Stream_Order]
#test for serial correlation 
mod2 = lm(Perc_Lost_Residuals ~ Perc_Lost_Resd_Lags, data = Summary2)
summary(mod2)

#Plot of straight residuals 
ErosionResiduals <- ggplot(Summary2, aes(x= End_Year, y=Perc_Lost_Residuals)) +
  geom_point(size = 2, pch = 21, fill = "#63589f") +
  facet_wrap(~Stream_Order) +
  theme_pubr +
  ggtitle("Erosion Residuals vs. Time \n") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = 'Year',
    y = 'Residual (spline fit)'
  )

print(ErosionResiduals)

#determine significance of residuals versus their lags
for(i in 1:5){
    Strm_ord_select <- Strm_Ord[i]
    riv_data <- Summary2[Stream_Order == Strm_ord_select]
    #print(riv_data)
    if(nrow(riv_data)>3){
      mod2 = lm(Perc_Lost_Residuals ~ Perc_Lost_Resd_Lags, data = riv_data)
      summary_table_select <- data.table(tidy(mod2))[2,.(estimate, p.value)]
      summary_table_select[, ":=" (Adj.R = summary(mod2)$r.squared)]
      summary_table_select[, ":=" (Stream_Order = Strm_ord_select, 
                                   significance = ifelse(p.value < 0.01, '**',
                                                         ifelse(p.value <0.05, "*", "")))]
      if(i==1){
        summary_table_Resd <- summary_table_select
      }else{
        summary_table_Resd <- rbind(summary_table_Resd, summary_table_select, use.names=T)
      }
    }
  }

lags <- ggplot(Summary2, aes(x= Perc_Lost_Residuals, y=Perc_Lost_Resd_Lags)) +
  geom_smooth(method = 'lm', lty = 'dashed', linewidth = 1, color = 'black') +
  geom_point(size = 2, pch = 21, color = "#d1afe8") +
  facet_wrap(~Stream_Order) +
  theme_pubr +
  geom_text(data = summary_table_Resd, aes(x=-0.3, y=0.5, label = paste0('R Sq: ', round(Adj.R, 3), significance)), hjust = 0) +
  ggtitle("Erosion Residuals vs. First Lags") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = 'Residual',
    y = 'Lag'
  )

print(lags)

hist(Summary1$Annual_DegD_Anomaly) #approximately Normal 
hist(Summary1$Perc_Lost_Anomaly) #normal
hist(Summary1$Q_MEAN_Anomaly) #normal 

#Pearson Correlation between Area Lost (by Stream Order) and Degree Days
#first, on the raw data - which is normally distributed (see plots above) and not serially correlated (see plots above)
Erosion_DegD_Cor_Strm <- Summary2[, cor.test(Annual_DegD_Anomaly, Perc_Lost_Anomaly, method = "pearson"), by = 'Stream_Order']
Erosion_DegD_Cor_Strm[, conf.int := NULL]
Erosion_DegD_Cor_Strm <- unique(Erosion_DegD_Cor_Strm)

#Pearson Correlation between Area Lost (by Stream Order) and Discharge
#first, on the raw data - which is normally distributed (see plots above) and not serially correlated (see plots above)
Erosion_Q_Cor_StrmOrd <- Summary2[, cor.test(Q_MEAN_Anomaly, Perc_Lost_Anomaly, method = "pearson"), by = 'Stream_Order']
Erosion_Q_Cor_StrmOrd[, conf.int := NULL]
Erosion_Q_Cor_StrmOrd <- unique(Erosion_Q_Cor_StrmOrd)


#################################################
#                                               #
###             ANOMALY PLOTS                ####
#                                               #
#################################################

coeff <- 10 # for scaling the secondary axis values

CarlPlot <- ggplot(Summary2, aes(x = End_Year, y = Perc_Lost_Anomaly, color = "Stream_Order")) +
  geom_point(aes(color = Stream_Order), pch = 20, size = 4) +
  geom_line(aes(group = Stream_Order, color = Stream_Order),  linewidth = 2) +
  geom_errorbar(aes(ymin = (Perc_Lost_Anomaly - Perc_Lost_Anomaly_SE), ymax = (Perc_Lost_Anomaly + Perc_Lost_Anomaly_SE)), color = 'black', width = 0.2, alpha = 0.5) +
  facet_wrap(~Stream_Order)+
  #geom_line(aes(x=End_Year, y=Q_MEAN_Anomaly/coeff1), color = "#049DD9", linewidth = 2, alpha = 0.4) + #scaling by 10 for second axis
  geom_line(aes(x=End_Year, y=Annual_DegD_Anomaly/coeff1, group = Stream_Order), color = "#F21905", linewidth = 2, alpha = 0.4) + # scaling by 10 for second axis
  geom_errorbar(aes(ymin = (Annual_DegD_Anomaly - Annual_DegD_Anomaly_SE)/coeff1, ymax = (Annual_DegD_Anomaly + Annual_DegD_Anomaly_SE)/coeff1), color = "black", width = 0.2, alpha = 0.5) +
  scale_color_carto_c(name = "Stream order: ",
                      type = "Sequential", palette = "TealGrn", direction = -1) +
  scale_y_continuous(
    name = "Erosion Anomaly \n",
    sec.axis = sec_axis(~.*coeff1, name="Discharge and Degree Day Anomalies \n")
    ) +
  theme_pubr +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1995, y=0.5, label = paste0('Pearson Coeff: ', round(estimate, 3))), hjust = 0) +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1991, y=0.4, label = paste0('p-value: ', round(p.value,3))), hjust = 0) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = 'black'), 
    strip.text.y = element_text(size = 15, face = "bold"),
    legend.position = ,
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="bold"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=20, color = 'black')
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  labs(
    #title = "Normalized Area Lost vs. Time",
    #subtitle = "by Stream Order and Latitude",
    x = '',
    y = '' #Area Lost Anomaly \n' #\n (Land Lost / Total Water at Start) \n'
  )

print(CarlPlot)

#################################################
#                                               #
###         CORRELATION ANALYSIS             ####
#               By Latitude                     #
#################################################

#Pearson Correlation between Area Lost (by Latitude Bins) and Degree Days
#first, on the raw data - which is normally distributed (see plots above) and not serially correlated (see plots above)
Erosion_DegD_Cor_LAT <- Summary3[, cor.test(Annual_DegD_Anomaly, Perc_Lost_Anomaly, method = "pearson"), by = 'N_S_2deg']
Erosion_DegD_Cor_LAT[, conf.int := NULL]
Erosion_DegD_Cor_LAT <- unique(Erosion_DegD_Cor_LAT)

#Pearson Correlation between Area Lost (by Latitude Bins) and DISCHARGE
#first, on the raw data - which is normally distributed (see plots above) and not serially correlated (see plots above)
Erosion_Q_Cor_LAT <- Summary3[, cor.test(Q_MEAN_Anomaly, Perc_Lost_Anomaly, method = "pearson"), by = 'N_S_2deg']
Erosion_Q_Cor_LAT[, conf.int := NULL]
Erosion_Q_Cor_LAT <- unique(Erosion_Q_Cor_LAT)

##### ANOMALY PLOT BY LATITUDE BIN

coeff <- 10 # for scaling the secondary axis values

CarlPlot_Lat <- ggplot(Summary3, aes(x = End_Year, y = Perc_Lost_Anomaly, color = "N_S_2deg")) +
  geom_point(aes(color = N_S_2deg), pch = 20, size = 3) +
  geom_line(aes(group = N_S_2deg, color = N_S_2deg),  linewidth = 1) +
  geom_errorbar(aes(ymin = (Perc_Lost_Anomaly - Perc_Lost_Anomaly_SE), ymax = (Perc_Lost_Anomaly + Perc_Lost_Anomaly_SE)), color = 'black', width = 0.2, alpha = 0.5) +
  facet_wrap(~N_S_2deg)+
  #geom_line(aes(x=End_Year, y=Q_MEAN_Anomaly/coeff1), color = "#049DD9", linewidth = 2, alpha = 0.4) + #scaling by 10 for second axis
  geom_line(aes(x=End_Year, y=Annual_DegD_Anomaly/coeff1, group = N_S_2deg), color = "#F21905", linewidth = 2, alpha = 0.4) + # scaling by 10 for second axis
  geom_errorbar(aes(ymin = (Annual_DegD_Anomaly - Annual_DegD_Anomaly_SE)/coeff1, ymax = (Annual_DegD_Anomaly + Annual_DegD_Anomaly_SE)/coeff1), color = "black", width = 0.2, alpha = 0.5) +
  #scale_color_carto_c(name = "Latitude Bin: ",
                      #type = "diverging", palette = "Earth", direction = -1) +
  scale_y_continuous(
    name = "Erosion Anomaly \n",
    sec.axis = sec_axis(~.*coeff1, name="Discharge and Degree Day Anomalies \n")
    #limits = c(-1.2,1.2)
  ) +
  theme_pubr +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1995, y=0.5, label = paste0('Pearson Coeff: ', round(estimate, 3))), hjust = 0) +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1991, y=0.4, label = paste0('p-value: ', round(p.value,3))), hjust = 0) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"), 
    strip.text.y = element_text(size = 15, face = "bold", color = "black"),
    legend.position = ,
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="plain"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=20, face = "plain", color = "black")
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  labs(
    #title = "Normalized Area Lost vs. Time",
    #subtitle = "by Stream Order and Latitude",
    x = '',
    y = '' #Area Lost Anomaly \n' #\n (Land Lost / Total Water at Start) \n'
  )

print(CarlPlot_Lat)


#################################################
#                                               #
###         CORRELATION ANALYSIS             ####
#           By Permafrost Extent                #
#################################################

#Pearson Correlation between Area Lost (by Permafrost Extent) and Degree Days
#first, on the raw data - which is normally distributed (see plots above) and not serially correlated (see plots above)
Erosion_DegD_Cor_PF <- Summary4[, cor.test(Annual_DegD_Anomaly, Perc_Lost_Anomaly, method = "pearson"), by = 'PF_EXTENT']
Erosion_DegD_Cor_PF[, conf.int := NULL]
Erosion_DegD_Cor_PF <- unique(Erosion_DegD_Cor_PF)

#Pearson Correlation between Area Lost (by Permafrost Extent) and DISCHARGE
#first, on the raw data - which is normally distributed (see plots above) and not serially correlated (see plots above)
Erosion_Q_Cor_PF <- Summary4[, cor.test(Q_MEAN_Anomaly, Perc_Lost_Anomaly, method = "pearson"), by = 'PF_EXTENT']
Erosion_Q_Cor_PF[, conf.int := NULL]
Erosion_Q_Cor_PF <- unique(Erosion_Q_Cor_PF)

coeff <- 10 # for scaling the secondary axis values
Summary4$PF_EXTENT <- factor(Summary4$PF_EXTENT,levels=c("Continuous","Discontinuous","Sporadic", "Isolated", "None"))

#plot
CarlPlot_PF <- ggplot(Summary4, aes(x = End_Year, y = Perc_Lost_Anomaly, color = "PF_EXTENT")) +
  geom_point(aes(color = PF_EXTENT), pch = 20, size = 3) +
  geom_line(aes(group = PF_EXTENT, color = PF_EXTENT),  linewidth = 1) +
  geom_errorbar(aes(ymin = (Perc_Lost_Anomaly - Perc_Lost_Anomaly_SE), ymax = (Perc_Lost_Anomaly + Perc_Lost_Anomaly_SE)), color = 'black', width = 0.2, alpha = 0.5) +
  facet_wrap(~PF_EXTENT)+
  #geom_line(aes(x=End_Year, y=Q_MEAN_Anomaly/coeff1), color = "#049DD9", linewidth = 2, alpha = 0.4) + #scaling by 10 for second axis
  geom_line(aes(x=End_Year, y=Annual_DegD_Anomaly/coeff1, group = PF_EXTENT), color = "#F21905", linewidth = 2, alpha = 0.4) + # scaling by 10 for second axis 
  geom_errorbar(aes(ymin = (Annual_DegD_Anomaly - Annual_DegD_Anomaly_SE)/coeff1, ymax = (Annual_DegD_Anomaly + Annual_DegD_Anomaly_SE)/coeff1), color = "black", width = 0.2, alpha = 0.5) +
  #scale_color_carto_c(name = "Latitude Bin: ",
  #type = "diverging", palette = "Earth", direction = -1) +
  scale_y_continuous(
    name = "Erosion Anomaly \n",
    sec.axis = sec_axis(~.*coeff1, name="Degree Day Anomaly \n")
    #limits = c(-1.2,1.2)
  ) +
  theme_pubr +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1995, y=0.5, label = paste0('Pearson Coeff: ', round(estimate, 3))), hjust = 0) +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1991, y=0.4, label = paste0('p-value: ', round(p.value,3))), hjust = 0) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = 'black'), 
    strip.text.y = element_text(size = 15, face = "bold", color = 'black'),
    legend.position = ,
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="plain"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=20, color = "black")
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  labs(
    #title = "Normalized Area Lost vs. Time",
    #subtitle = "by Stream Order and Latitude",
    x = '',
    y = '' #Area Lost Anomaly \n' #\n (Land Lost / Total Water at Start) \n'
  )

print(CarlPlot_PF)


#################################################
#                                               #
###         CORRELATION ANALYSIS             ####
#           By Mean Thaw Depth                  #
#################################################

#Pearson Correlation between Area Lost (by Permafrost Extent) and Degree Days
#first, on the raw data - which is normally distributed (see plots above) and not serially correlated (see plots above)
Erosion_DegD_Cor_MTD <- Summary5[, cor.test(Annual_DegD_Anomaly, Perc_Lost_Anomaly, method = "pearson"), by = 'MTD_bins']
Erosion_DegD_Cor_MTD[, conf.int := NULL]
Erosion_DegD_Cor_MTD <- unique(Erosion_DegD_Cor_MTD)


##### ANOMALIES 

coeff <- 10 # for scaling the secondary axis values
Summary5$MTD_bins <- factor(Summary5$MTD_bins,levels=c("(2.5,3.0]","(3.0,3.5]", "(3.5,4.0]", "(4.0,4.5]"))

CarlPlot_MTD <- ggplot(Summary5, aes(x = End_Year, y = Perc_Lost_Anomaly, color = "MTD_bins")) +
  geom_point(aes(color = MTD_bins), pch = 20, size = 3) +
  geom_line(aes(group = MTD_bins, color = MTD_bins),  linewidth = 1) +
  geom_errorbar(aes(ymin = (Perc_Lost_Anomaly - Perc_Lost_Anomaly_SE), ymax = (Perc_Lost_Anomaly + Perc_Lost_Anomaly_SE)), color = 'black', width = 0.2, alpha = 0.5) +
  facet_wrap(~MTD_bins)+
  #geom_line(aes(x=End_Year, y=Q_MEAN_Anomaly/coeff1), color = "#049DD9", linewidth = 2, alpha = 0.4) + #scaling by 10 for second axis
  geom_line(aes(x=End_Year, y=Annual_DegD_Anomaly/coeff1, group = MTD_bins), color = "#F21905", linewidth = 2, alpha = 0.4) + # scaling by 10 for second axis 
  geom_errorbar(aes(ymin = (Annual_DegD_Anomaly - Annual_DegD_Anomaly_SE)/coeff1, ymax = (Annual_DegD_Anomaly + Annual_DegD_Anomaly_SE)/coeff1), color = "black", width = 0.2, alpha = 0.5) +
  #scale_color_carto_c(name = "Latitude Bin: ",
  #type = "diverging", palette = "Earth", direction = -1) +
  scale_y_continuous(
    name = "Erosion Anomaly \n",
    sec.axis = sec_axis(~.*coeff1, name="Degree Day Anomaly \n"),
    #xlimits = c(-1,1.5)
  ) +
  theme_pubr +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1995, y=0.5, label = paste0('Pearson Coeff: ', round(estimate, 3))), hjust = 0) +
  #geom_text(data = Erosion_DegD_Cor, aes(x=1991, y=0.4, label = paste0('p-value: ', round(p.value,3))), hjust = 0) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = 'black'), 
    strip.text.y = element_text(size = 15, face = "bold", color = 'black'),
    legend.position = ,
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="plain"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=20, color = "black")
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  labs(
    #title = "Normalized Area Lost vs. Time",
    #subtitle = "by Stream Order and Latitude",
    x = '',
    y = '' #Area Lost Anomaly \n' #\n (Land Lost / Total Water at Start) \n'
  )

print(CarlPlot_MTD)


######################## CRUCIAL STATISTICAL ANALYSES #############################################################################
## Correlation between Percent Lost Anomaly and Degree Day Anomaly across ALL SITES:
#data table with just annual mean: 
#b0f2bc,#89e8ac,#67dba5,#4cc8a3,#38b2a3,#2c98a0,#257d98 #TealGrn
#A16928,#bd925a,#d6bd8d,#edeac2,#b5c8b8,#79a7ac,#2887a1 #earth

#Pearson Correlation 
Erosion_DegD_Cor_ALL <- Summary1[, cor.test(Annual_DegD_Anomaly, Perc_Lost_Anomaly, method = "pearson")]
print(Erosion_DegD_Cor_ALL)
# Degrees of Freedom are not correctly accounted for as n-2 because we are doing a 5 year comparison window 
  # (we compare '85 to '90, '86 to '91, etc), so df should be 33/5
#first, calculate t statistic 
t = (0.54*((33/5)^0.5)) / (1- ((0.54)^2)) #divided by 5 -- for taking mean of window
t1 = (0.54*((34/2)^0.5)) / (1- ((0.54)^2)) #divided by 2 -- for comparing end points of window 

Erosion_Q_Cor_ALL <- Summary1[, cor.test(Q_MEAN_Anomaly, Perc_Lost_Anomaly, method = "pearson")]
Erosion_NDVI_Cor_ALL <- Summary1[, cor.test(mean_NDVI_Anomaly, Perc_Lost_Anomaly, method = "pearson")]
T_NDVI_Cor_ALL <- Summary1[, cor.test(Annual_DegD_Anomaly, mean_NDVI_Anomaly, method = "pearson")]
Erosion_SSC_Cor_ALL <- Summary1[, cor.test(mean_SSC_Anomaly, Perc_Lost_Anomaly, method = "pearson")]


## EVALUATING SIGNIFICANCE EVEN MORE RIGOROUSLY: MONTE CARLO SIMULATION
# Monte Carlo Simulation to get an actual sense of error that is ambiguous due to the degrees of freedom question
#find mean and SD of my annual Erosion anomlay data
mean_ErosAnom <- mean(Summary1$Perc_Lost_Anomaly)
sd_ErosAnom <- sd(Summary1$Perc_Lost_Anomaly)

x<-rnorm(34, mean = mean_ErosAnom, sd = sd_ErosAnom) #33 random values from a normal distribution
x # 34 numbers
w <- Summary1$Annual_DegD_Anomaly # also 34 numbers 
testy_test <- cor.test(w, x, method = "pearson")
testy_test$estimate


#first test - % of time that a random distribution with the mean and std dev of the annual erosion anomalies \
## is significant at p<0.05 for a two-sided t-test 
res<-replicate(1000, {
  y<-rnorm(34, mean = mean_ErosAnom, sd = sd_ErosAnom)
  R <- cor.test(w, y, method = "pearson")
  z<- R$p.value # find maximum correlation for two-sided t-test, hence abs()
  sum(z<0.05)#tests if the p value is less that 0.05
})
mean(res>0)#find the proportion of trials where at least some windows appeared to have a significant correlation

#Second Test, Monte Carlo 
res<-replicate(1000, {
  y<-rnorm(34, mean = mean_ErosAnom, sd = sd_ErosAnom)
  R <- cor.test(w, y, method = "pearson")
  z<-max(abs(R$estimate))# find maximum correlation for two-sided t-test, hence abs()
})
quantile(res, 0.975) #critical of maximum correlation at p=0.05 two-sided test.

#Third test - Monte Carlo CDF
res<-replicate(1000, {
  y<-rnorm(34, mean = mean_ErosAnom, sd = sd_ErosAnom)
  R <- cor.test(w, y, method = "pearson")
  z<-max(abs(R$estimate))# find maximum correlation for two-sided t-test, hence abs()
})
res
quantile(res, 0.975) #critical of maximum correlation at p=0.05 two-sided test.


### Wrangle Monte Carlo Data ###
res <- data.table(res)
MC <- res[order(res)]
MC[, Index := 1:.N] #index column w/ count for the CDF

#create the probability column for the cumulative distribution function
MC[, Divisor := max(Index)] 
MC[, Probability := Index/Divisor] 
setnames(MC, 'res', 'cor_coef')

#Plot a cumulative distribution function for the correlation coefficient
CDF_Monte_Carlo<-ggplot(MC,
                         aes(x = cor_coef, y = Probability)) +
  theme_pubr() +
  geom_point(color = "#d1afe8", size = 5, pch = 21) +
  #geom_line() +
  #scale_x_log10() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, color = 'black'), 
    strip.text.y = element_text(size = 20, face = "bold", color = 'black'),
    legend.position = "none",
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="bold"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=20, face= 'bold', color = 'black')
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  #stat_regline_equation(label.y = .75, size=2) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  # the next line adds the slope and p-values to each facet
  #geom_text(data = summary_table, aes(x=1991, y=0.4, label = paste0(round(estimate, 3), '%/yr ', significance)), hjust = 0) + 
  labs(
    title = "Monte Carlo Simulation CDF",
    subtitle = "1000 runs of a random normal distribution of annual erosion anomaly \n compared to annual degree day anomaly distribution",
    x = '\n Correlation Coefficient ',
    y = 'Probability \n'
  )

print(CDF_Monte_Carlo)


######################## END CRUCIAL STATISTICAL ANALYSES ##############################################################################
  

  

  
##### FUTURE PROJECTIONS AND HISTORICAL DATA PLOT -- FIGURE 2 in MAIN TEXT OF FIELDS ET AL. #####  

coeff1 = 5
#plot itself -- Note that the regression lines here are still based on all the data, not the annual means, which are plotted as points
## Change input to percent area lost or real area lost, depending on what you want 
Anomalies_AllYrs <- ggplot(Summary_F1, aes(x = End_Year, y = Perc_Lost_Anomaly)) +
  #geom_point(size = 6, pch = 21, alpha = 0.8) +
  geom_line(color = "#bd925a", linewidth = 4) +
  geom_errorbar(aes(ymin = (Perc_Lost_Anomaly - Perc_Lost_Anomaly_SE), ymax = (Perc_Lost_Anomaly + Perc_Lost_Anomaly_SE)), width = 0.2) +
  #geom_line(aes(x=End_Year, y= mean_NDVI_Anomaly/coeff1), color = "#6EA646" , linewidth = 4, alpha =0.8) +
  #geom_line(aes(x=End_Year, y=Q_MEAN_Anomaly/coeff1), color = "#049DD9", linewidth = 2, alpha = 0.4) + #scaling by 10 for second axis
  geom_line(aes(x=End_Year, y=Annual_DegD_Anomaly/coeff1), color = "#ee4d5a", linewidth = 4, alpha = 0.8) + # scaling by 10 for second axis,  #F21905"
  geom_errorbar(aes(ymin = (Annual_DegD_Anomaly - Annual_DegD_Anomaly_SE)/coeff1, ymax = (Annual_DegD_Anomaly + Annual_DegD_Anomaly_SE)/coeff1), width = 0.2) +
  geom_line(aes(x=End_Year, y=F_Annual_DegD_Anomaly/coeff1), color = "#ee4d5a", linewidth = 1.5, alpha = 0.6) +
  geom_errorbar(aes(ymin = (F_Annual_DegD_Anomaly - F_Annual_DegD_Anomaly_SE)/coeff1, ymax = (F_Annual_DegD_Anomaly + F_Annual_DegD_Anomaly_SE)/coeff1), width = 0.2) +
  geom_line(aes(x=End_Year, y=predicted_Eros_Anomaly), color = "#bd925a", linewidth = 1.5) +
  geom_ribbon(aes(ymin=predicted_Eros_Anomaly-0.1679565, ymax=predicted_Eros_Anomaly+0.1679565), fill = "#bd925a", linewidth = 1.5, alpha = 0.6) + #using  RMSE = 0.1679565 as error for future erosion
  #scale_color_carto_c(name = "Latitude Bin: ",
  scale_y_continuous(
    #name = "Erosion Anomaly \n",
    sec.axis = sec_axis(~.*coeff1, name="Degree Day Anomaly \n")
  ) +
  theme_pubr +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = 'black'), 
    strip.text.y = element_text(size = 20, face = "bold", color = 'black'),
    legend.position = "none",
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="bold"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=20, face= 'bold', color = 'black')
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  #stat_regline_equation(label.y = .75, size=2) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  # the next line adds the slope and p-values to each facet
  #geom_text(data = summary_table, aes(x=1991, y=0.4, label = paste0(round(estimate, 3), '%/yr ', significance)), hjust = 0) + 
  labs(
    title = "Percent Area Lost vs. Time",
    subtitle = "Mean annual values across all sites",
    x = '',
    y = 'Percent Area Lost Anomaly \n' #\n (Land Lost / Total Water at Start) \n'
  ) 

print(Anomalies_AllYrs)

## NDVI or Q VS EROSION ANOMALIES -- Figure 1b and 1c in Fields et al., 
  # toggle on/off various lines to plot different factors (NDVI, Q, etc)
NDVI_erosion_allYrs <- ggplot(Summary1, aes(x = End_Year, y = Perc_Lost_Anomaly)) +
  #geom_point(size = 6, pch = 21, alpha = 0.8) +
  geom_line(color = "#bd925a", linewidth = 4) +
  geom_errorbar(aes(ymin = (Perc_Lost_Anomaly - Perc_Lost_Anomaly_SE), ymax = (Perc_Lost_Anomaly + Perc_Lost_Anomaly_SE)), width = 0.2) +
  geom_line(aes(x=End_Year, y= Q_MEAN_Anomaly/coeff1), color = "#049DD9" , linewidth = 4, alpha =0.8) +
  #geom_line(aes(x=End_Year, y=mean_SSC_Anomaly), color = "grey", linewidth = 4, alpha = 0.8) + 
  #geom_errorbar(aes(ymin = (mean_SSC_Anomaly - se_SSC_Anomaly), ymax = (mean_SSC_Anomaly + se_SSC_Anomaly)), width = 0.2) +
  geom_line(aes(x=End_Year, y= mean_NDVI_Anomaly/coeff1), color = "#6EA646" , linewidth = 4, alpha =0.8) +
  geom_errorbar(aes(ymin = (mean_NDVI_Anomaly - se_NDVI_Anomaly)/coeff1, ymax = (mean_NDVI_Anomaly + se_NDVI_Anomaly)/coeff1), width = 0.2) +
  #geom_line(aes(x=End_Year, y=Q_MEAN_Anomaly/coeff1), color = "#049DD9", linewidth = 2, alpha = 0.4) + #scaling by 10 for second axis
  #geom_line(aes(x=End_Year, y=Annual_DegD_Anomaly/coeff1), color = "#ee4d5a", linewidth = 4, alpha = 0.8) + # scaling by 10 for second axis,  #F21905"
  #geom_errorbar(aes(ymin = (Annual_DegD_Anomaly - Annual_DegD_Anomaly_SE)/coeff1, ymax = (Annual_DegD_Anomaly + Annual_DegD_Anomaly_SE)/coeff1), width = 0.2) +
  #geom_line(aes(x=End_Year, y=F_Annual_DegD_Anomaly/coeff1), color = "#ee4d5a", linewidth = 1.5,linetype = 'dotted') +
  #geom_errorbar(aes(ymin = (F_Annual_DegD_Anomaly - F_Annual_DegD_Anomaly_SE)/coeff1, ymax = (F_Annual_DegD_Anomaly + F_Annual_DegD_Anomaly_SE)/coeff1), width = 0.2) +
  #geom_line(aes(x=End_Year, y=predicted_Eros_Anomaly), color = "#bd925a", linewidth = 1.5,linetype = 'dotted') +
  #scale_color_carto_c(name = "Latitude Bin: ",
  scale_y_continuous(
    #name = "Erosion Anomaly \n",
    sec.axis = sec_axis(~., name="NDVI and Discharge Anomalies \n")
  ) +
  theme_pubr +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = 'black'), 
    strip.text.y = element_text(size = 20, face = "bold", color = 'black'),
    legend.position = "none",
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="bold"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=20, face= 'bold', color = 'black')
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  #stat_regline_equation(label.y = .75, size=2) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  # the next line adds the slope and p-values to each facet
  #geom_text(data = summary_table, aes(x=1991, y=0.4, label = paste0(round(estimate, 3), '%/yr ', significance)), hjust = 0) + 
  labs(
    title = "Erosion, NDVI, and Discharge vs. Time",
    subtitle = "Mean annual values across all sites",
    x = '',
    y = 'Percent Area Lost & SSC Anomaly \n' #\n (Land Lost / Total Water at Start) \n'
  ) 

print(NDVI_erosion_allYrs)



### Single Thiessen example Plot for Case Study 

#subset Master1 to just include Thiessen 6331463

R1 <- Master1[Thiessen_ID == 6331463]

R1 <- R1[, Perc_Lost_Anomaly := ((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006)]
R1 <- R1[, Perc_Lost_Anomaly_SE := (sd((Percent_Area_Lost-mean_Area_Lost_1990_2006)/SD_Area_Lost_1990_2006))/sqrt(length(Percent_Area_Lost))]
R1 <- R1[, DegDay_Anomaly := (degreeDays-mean_degreeDays_79_95)/sd_degreeDays_79_95] 
R1 <- R1[, Above0_Anomaly := (above0Days-mean_above0_79_95)/sd_above0_79_95]

#Meander Data: This is from manaully drawing lines and measuring the offset in stream centerline year to year in illustrator/GEE... very annoying, takes a long time 
CS3_m <- c(110, 95, 5, 15, 20, 180);
CS3_m1 <- c(1997, 2002, 2007, 2011, 2017, 2023); 

CS3_new <- data.table(cbind(CS3_m, CS3_m1)); 
CS3_new <- setnames(CS3_new, c('CS3_m', 'CS3_m1'), c('Meander', 'Year'))
CS3_new <- CS3_new[, Perc_Meander := c(0.137, 0.097, 0.003, 0.012, 0.014, 0.135)]
CS3_new <- rbind(CS3_new, list(100, 1992, 0.101))
CS3_new <- CS3_new[, Meander := Meander/5] # change units to m/yr
CS3_new <- CS3_new[, Perc_Meander := Perc_Meander*100] # change units to m/yr

#organized meander data:
setwd(paste0(wd_root))
Meander_Data <- fread("Case_Study/Site3/Manually_Measured_Data.csv")

Meander_Data[, Width_pt := NULL]
Meander_Data[, Width_m := NULL]
Meander_Data[, Meander_pt := NULL]
Meander_Data[, Position := NULL]

Meander_summary <- Meander_Data[, list( mean_Meander = mean(Meander_yr),
                                        mean_Perc_Meander = mean(Perc_Meander)), 
                                        by = .(Year)]
                                

#Case Study Site 3 Plot
CS3 <- ggplot(R1, aes(x = End_Year, y = Percent_Area_Lost)) +
  geom_point(pch = 19, color = '#bd925a', size = 4) +
  geom_line(linewidth = 1, color = '#bd925a') +
  geom_smooth(method = 'lm', lty = 'dashed', linewidth = 2, color = '#0554F2') + #for the Percent Area Lost Plot 
  #geom_line(aes(x=End_Year, y=DegDay_Anomaly), color = "#F21905", linewidth = 2) +
  scale_y_continuous(
    #name = "Erosion Anomaly \n",
    #sec.axis = sec_axis(~./coeff1, name="Degree Day Anomaly \n") #~./coeff1,
  ) +
  theme_pubr +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    strip.text.y = element_text(size = 15, face = "bold"),
    legend.position = "none",
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="bold"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=18, face = "bold", color = 'black')
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  #stat_regline_equation(label.y = .75, size=2) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  # the next line adds the slope and p-values to each facet
  #geom_text(data = summary_table, aes(x=1991, y=0.4, label = paste0(round(estimate, 3), '%/yr ', significance)), hjust = 0) + 
  labs(
    title = "Percent Area Lost vs. Time",
    subtitle = "Thiessen Polygon 6331463",
    x = '',
    y = 'Percent Area Lost \n' #\n (Land Lost / Total Water at Start) \n'
  ) 

print(CS3)


#Case Study Site 3 Plot

CS3_Meander <- ggplot(Meander_summary, aes(x = Year, y = mean_Meander)) +
  geom_point(pch = 21, color = 'black', size = 6, stroke = 3) +
  geom_line(linewidth = 4, color = '#bd925a') +
  geom_smooth(method = 'lm', lty = 'dashed', linewidth = 2, color = '#0554F2') + #for the Percent Area Lost Plot 
  #geom_line(aes(x=End_Year, y=DegDay_Anomaly), color = "#F21905", linewidth = 2) +
  scale_y_continuous(
    #name = "Erosion Anomaly \n",
    #sec.axis = sec_axis(~./coeff1, name="Degree Day Anomaly \n") #~./coeff1,
  ) +
  theme_pubr +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    strip.text.y = element_text(size = 15, face = "bold"),
    legend.position = "none",
    #title=element_text(size=22,face="bold"),
    axis.title = element_text(size=20, face="bold"),
    #subtitle = element_text(size=12, face = "plain"),
    axis.text = element_text(size=18, face = "bold", color = 'black')
  )+
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(strip.text.x.top = element_text(size = 14, face = "bold")) +
  #stat_regline_equation(label.y = .75, size=2) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(face = "plain", size =20)) +
  # the next line adds the slope and p-values to each facet
  #geom_text(data = summary_table, aes(x=1991, y=0.4, label = paste0(round(estimate, 3), '%/yr ', significance)), hjust = 0) + 
  labs(
    title = "Manually Measured Meander vs. Time",
    subtitle = "Thiessen Polygon 6331463",
    x = '',
    y = 'Meander (m/yr)' #\n (Centerline offset / Stream Width) \n'
  ) 

print(CS3_Meander)




#### FIN. ####