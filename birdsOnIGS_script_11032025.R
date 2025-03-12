library(auk)
library(sf)
library(tidyverse)
library(janitor)
library(readxl)
library(arrow)
library(terra)
library(glmmTMB)

################################################################################

# download ebird data #

################################################################################

###### prep objects required to extract data from auk ######

### create object of area data is to be extracted from

## load case study area

sa1_estUrb_clip <- st_read("inputData/sa1_estUrb_clip.shp") %>% 
  # change coordinates to 4326 to be consistent with auk system.
  st_transform(4326)

## create object of bounding box of case study area
# only simple box polygon can be used as input

bbox <- st_bbox(sa1_estUrb_clip)

coords <- matrix(c(bbox[1], bbox[2],
                   bbox[1], bbox[4],
                   bbox[3], bbox[4],
                   bbox[3], bbox[2],
                   bbox[1], bbox[2]), ncol = 2, byrow = T)

bbox_sf <- st_sfc(st_polygon(list(coords)), crs = 4326) %>%   # needs to be 4326 to be consistent with auk system
  st_sf()

### create input and output files as required by auk

# text file downloaded from https://ebird.org/data/download
input_file <- "inputData/ebd_AU-VIC_smp_relMay-2024.txt"

output_file <- paste0("eBird_output/eBird_filtered", Sys.Date(), ".txt")

###### extract eBird data ######

### prep filters for extracting

ebd <- input_file %>%
  auk_ebd() %>% 
  auk_date(date = c("2014-01-01", "2024-05-31")) %>%
  auk_bbox(bbox = bbox) %>% 
  auk_duration(duration = c(5,240)) %>%  # consistent with Callaghan et al (2019). Filter added to control for citizen science nature of eBird.
  # only include complete lists
  auk_complete() %>% 
  # only survey protocols with limited travelling distance to be included to ensure that observations occurred on a given site 
  auk_protocol(c("Stationary", "Traveling"))

### execute filters to extract specified data
ebd <- ebd %>% 
  auk_filter(file = output_file, overwrite = T) %>% 
  read_ebd(unique = T, rollup = T)

### apply further filters that couldn't be applied until in df format

ebd <- ebd %>% 
  filter(
    (protocol_type == "Stationary") |
      # effort distance of 30m ensures that observations occurred on the recorded site
      (protocol_type == "Traveling" & effort_distance_km <= 0.03))

## create nativeness field and remove escapee and exotic species
# coded like this to make it easier if you want to include exotics into consideration

ebd <- ebd %>% 
  mutate(nativity = case_when(is.na(exotic_code) ~ "Native",
                              exotic_code == "N" ~ "Exotic",
                              exotic_code == "X" ~ "Escaped")) %>% 
  filter(!(nativity %in% c("Escaped", "Exotic")))

###### Convert data to an sf format ######

### create function to prevent having to save geometry to the environment
txt_to_sf <- function(obj){
  
  geometry <- lapply(1:nrow(obj),function(x){
    st_point(x = c(obj$longitude[x], obj$latitude[x]))}) %>% 
    st_sfc(crs = 4326) #needs to be in 4326 as this is what the coordinates are written in in the csv
  
  a <- obj %>% 
    mutate("geometry" = geometry) %>%
    st_sf()
  
  return(a)
}

### execute function

ebd <- txt_to_sf(ebd)

###### extract points in case study ######

### create function to avoid needing to save data in environment

clip_pts <- function(pts, clip){
  
  # create is.intersect column in clip to facilitate st_join function 
  clip <- mutate(clip, is.intersect = "intersect")
  
  # run st_join to only extract points within identified clip
  a <- st_join(pts, clip, join = st_intersects)
  
  b <- filter(a, !is.na(is.intersect))
  
  return(b)
}

### execute function

ebd <- clip_pts(ebd, sa1_estUrb_clip) %>% 
  # convert to EPSG:28355 to be consistent with all other data sets used
  st_transform(28355)

###### clean discrepancies in ebd data ######
# this is important as i will be splitting and re-joining vectors and need to
# make sure they align correctly in the first place to avoid bringing in errors.

### test for discrepancies
# make sure each survey only has one corresponding value in relevant vector

## create function

test_vectors <- function(sf, vector1, vector2){
  
  # convert sf to a df to facilitate analysis
  df <- as.data.frame(sf)
  
  # calculate number of unique values in vector 1
  a <- nrow(
    df %>% 
      group_by({{vector1}}) %>% 
      summarise() %>% 
      ungroup()
  )
  
  # calculate number of unique values in vector when grouped with secondary vector
  # if a value in vector 1 appears multiple times with different values in vector 2, it will create a new row in b, repeating the survey ID
  b <- nrow(
    df %>%
      group_by({{vector1}}, {{vector2}}) %>%
      summarise() %>%
      ungroup()
  )
  
  # make results easy to read
  results <- data.frame(names = c("unique vector 1", "unique vector1 + vector2", "match"),
                        results = c(a, b, a==b))
  
  print(results)
}

test_vectors(ebd, checklist_id, duration_minutes)
test_vectors(ebd, checklist_id, geometry)
test_vectors(ebd, common_name, scientific_name)
test_vectors(ebd, common_name, nativity)

### fix discrepancies
# duration_minutes showed the only discrepancy

## identify survey ID of discrepancy

discrep_df <- ebd %>% 
  as.data.frame() %>% 
  group_by(checklist_id, duration_minutes) %>% 
  summarise() %>% 
  ungroup()

identified_discrep <- duplicated(discrep_df$checklist_id)

checklist_for_removal <- as.character(discrep_df[identified_discrep, 'checklist_id'])

## remove survey from eBird data

ebd <- ebd %>% 
  filter(checklist_id != checklist_for_removal)

### remove irrelevant records
# records are of species identified as unrelated to IGS through a manual checking process

## remove pelagic seabird
# list based on list of present species

pelagic <- (c("Great Frigatebird", "Parasitic Jaeger", "Short-tailed Shearwater", "Australasian Gannet", "White-faced Storm-Petrel"))

ebd <- ebd %>% 
  filter(!(common_name %in% pelagic))

################################################################################

# allocate eBird point data to corresponding IGS areas #

################################################################################

###### prep ebd data for allocation  ######

### change ebd object name to facilitate subsequent processes

ebd_raw <- ebd

### convert to df to facilitate processing (it's more complicated as an sf)

ebd_raw <- ebd_raw %>% 
  as.data.frame()

###### create base structure of data frame ######

### create object with each species repeated for each survey
# each species will have a presence/absence value assocaited with it

## create object of all unique species
# every check list (survey) needs to have a row for each unique species
# common_name used as variable. Could be scientific name as both are equal/aligned

unique_sp <- ebd_raw %>%
  group_by(common_name) %>% 
  summarise() %>% 
  ungroup()

## create object of each unique survey

unique_srvy <- ebd_raw %>%
  group_by(checklist_id) %>% 
  summarise() %>% 
  ungroup()

## join unique_sp object to unique_srvy object
# expand.grid creates a data frame of all possible combination between listed vectors

ebd <- expand.grid(checklist_id = unique_srvy$checklist_id, common_name = unique_sp$common_name)

### re-join necessary vectors

## create object of vectors to be rejoined by survey ID
# can only join if there is one value per survey ID. scientific_name + nativity will be joined in subsiquent step.

rejoin_vectors1 <- ebd_raw %>% 
  group_by(checklist_id, duration_minutes, geometry) %>% 
  summarise() %>% 
  ungroup()

## create object of vectors to be rejoined by common name

rejoin_vectors2 <- ebd_raw %>% 
  group_by(common_name, scientific_name, nativity) %>% 
  summarise() %>% 
  ungroup()

## create object of vectors to be rejoined by chkls_ and common name

rejoin_vectors3 <- ebd_raw %>% 
  group_by(checklist_id, common_name, observation_count) %>% 
  summarise() %>% 
  ungroup()

## rejoin vectors

ebd <- ebd %>% 
  left_join(rejoin_vectors1, by = "checklist_id") %>% 
  left_join(rejoin_vectors2, by = "common_name") %>% 
  left_join(rejoin_vectors3, by = c("checklist_id", "common_name"))

## convert ebd back to sf object for spatial joins

ebd <-  ebd %>% 
  st_sf()

###### allocate variables data to each point ######

### allocate land sub-type classification to data

## load and process data

# load igs/fgs_site data 

custom_read <- function(name){
  st_read(paste0("inputData/igs_site/igsVeg18_", name,".shp"))
}

igs_ver <- custom_read("duplicatesRemoved_ver")
igs_brwn <- custom_read("duplicatesRemoved_brwn")
igs_util <- custom_read("duplicatesRemoved_util")
igs_rail <- custom_read("duplicatesRemoved_rail")
igs_otherIgs <- custom_read("otherIgs")
fgs_rpgs <- custom_read("rpgs")
fgs_otherFgs <- custom_read("otherFgs")
nonGs <- custom_read("nonGs")

# give each igs or fgs object a land_st value to facilitate allocation of land_st value via spatial merge

igs_ver <- mutate(igs_ver, land_st = "ver") %>% 
  select(!FID)
igs_brwn <- mutate(igs_brwn, land_st = "brwn") %>% 
  select(!FID)
igs_util <- mutate(igs_util, land_st = "util") %>% 
  select(!FID)
igs_rail <- mutate(igs_rail, land_st = "rail") %>% 
  select(!FID)
igs_otherIgs <- mutate(igs_otherIgs, land_st = "otherIgs") %>% 
  select(!FID)
fgs_rpgs <- mutate(fgs_rpgs, land_st = "rpgs") %>% 
  select(!FID)
fgs_otherFgs <- mutate(fgs_otherFgs, land_st = "otherFgs") %>% 
  select(!FID)
nonGs <- mutate(nonGs, land_st = "nonGs") %>% 
  select(!FID)

# amalgamate igs and fgs objects into single sf object to facilitate spatial join

land_st <- rbind(igs_ver, igs_brwn, igs_util, igs_rail, igs_otherIgs, fgs_rpgs, fgs_otherFgs, nonGs)

## allocate

ebd <- st_join(ebd, land_st, join = st_intersects)

### allocate land type classification to data

ebd <- ebd %>% 
  mutate(land_t = case_when(land_st %in% c("ver", "brwn", "util", "rail", "otherIgs") ~ "IGS",
                            land_st %in% c("rpgs", "otherFgs") ~ "formal",
                            land_st %in% c("nonGs") ~ "nonGs"))

### allocate neighbourhood to data

## load data
# data is SA1 boundaries from the Australian Bureau of Statistics

sa1 <- st_read("inputData/sa1_estUrb.shp") %>% 
  st_transform(28355) %>%
  st_make_valid() %>%
  clean_names("lower_camel") %>% 
  select(sa1Cod)

## allocate

ebd <- st_join(ebd, sa1, join = st_intersects)

### remove geometry vector
# geometry vector no longer needed once spatial joins have been executed

ebd <- ebd %>% 
  as.data.frame() %>% 
  select(!geometry)

### allocate foraging class and foraging domain and conservation status to eBird data

## load data

forage <- read_xlsx("inputData/speciesList.xlsx",
                    sheet = 1) %>% 
  rename(cnsvt_ = "ffgAct_june2024",
         frg_gld =  "foraging_guild",
         frg_dmn = "foraging_domain",
         common_name = "cmmn_nm") %>% 
  select(common_name, frg_gld, frg_dmn, cnsvt_)

## allocate

ebd <- ebd %>% 
  left_join(forage, by = "common_name")

## amend conservation status so that it is a binary between threatened and not-threatened

ebd <- ebd %>% 
  mutate(cnsvt_ = case_when(cnsvt_ %in% c("Critically endangered", "Endangered", "Vulnerable") ~ "threatened",
                            cnsvt_ == "Not threatened" ~ "Not threatened"))

## set threatened species as the reference term to facilitate building glmm
# threatened as the reference term because we want to know the difference for threatened species (kept constant) if we change the land type)

ebd <- ebd %>% 
  mutate(cnsvt_ = factor(cnsvt_, levels = c("threatened", "Not threatened")))

###### create and remove remaining vectors ######

### create vector of truncated presence/absence

ebd <- ebd %>% 
  mutate(presAbs = case_when(!is.na(observation_count) ~ 1,
                             is.na(observation_count) ~ 0)) %>% 
  select(!observation_count)

################################################################################

# calculate covarates

################################################################################

### convert sa1_estUrb_clip to 28355 to be consistant with all other data sources
# to be used throughout covariate calculations

sa1_estUrb_clip <- sa1_estUrb_clip %>% 
  st_transform(28355)

###### calculate tree per unit of analysis ######

### prepare vegetation raster data 

## load veg raster
# input data aggregated to 0.8m with height divisions as 0-50cm, 50cm - 3m, 3m+

veg1 <- rast("inputData/UM_vegData_melb2018_vht_raw/raw/veg18.tif")
veg2 <- rast("inputData/UM_vegData_melb2018_vht_raw/raw/veg18_periUrban_1.tif")
veg3 <- rast("inputData/UM_vegData_melb2018_vht_raw/raw/veg18_periUrban_2.tif")
veg4 <- rast("inputData/UM_vegData_melb2018_vht_raw/raw/veg18_periUrban_3.tif")

## combine all relevant raster layers (takes several hours and requires significant RAM)
veg <- mosaic(veg1, veg2, veg3, veg4)

## clip to established urban area boundary

# load case study boundary

sa1_estUrb_clip_sv <- vect(sa1_estUrb_clip)

# mask veg raster to only include veg WITHIN case study

veg <- mask(veg, sa1_estUrb_clip_sv, inverse = F)

# crop veg raster to exclude any areas outside of case study, speeding up calculations

estUrb_extent <- ext(sa1_estUrb_clip_sv)

veg <- crop(veg, estUrb_extent)

## reclassify into tree and ground cover groups (20 min)
# done to address limitation of veg data in observing midstory veg under canopies

tree <- veg
tree[tree == 1] <- NA
tree[tree %in% c(2,3)] <- 2

## aggregate data to facilitate later computations (several hours)

# create aggregating function so NA's don't impact aggregation, meaning aggregating outputs will not incorrectly diminish amount of tree

agg_fun_tree <- function(x) {
  ifelse(any(x == 2, na.rm = TRUE), 1, NA)
}

# aggregate groundcover and treecover data

tree <- aggregate(tree, fact = 2, fun = agg_fun_tree)

### calculate proportion of tree and ground cover per spatial unit

## make each cell's value equal the size of the cell (several hours)

tree <- cellSize(tree, mask = T, unit = "m")

## create object of area of each neighbourhood for use in calculations of densities

sa1_area <- sa1 %>% 
  mutate(sa1_area = st_area(.),
         sa1Cod = as.character(sa1Cod)) %>%  
  as.data.frame() %>% 
  select(sa1Cod, sa1_area)

## create spatialVector version of neighbourhoods for use with Raster function

sa1_sv <- vect(sa1)

## run calculations for tree cover (~ 1hr)

# split sa1 into 3 parts to facilitate extract
# current extract exceeds RAM limits of my computer and needs to be reduced in size

sa1_split_1 <- sa1[1:(nrow(sa1)/2),]
sa1_split_1_sv <- vect(sa1_split_1)

sa1_split_2 <- sa1[((nrow(sa1)/2)+1):nrow(sa1),]
sa1_split_2_sv <- vect(sa1_split_2)

# run calculations for each sa1_split area

total_tree_1 <- terra::extract(tree, sa1_split_1_sv) %>% 
  group_by(ID) %>% #needs to be ID as extract creates an ID value from which to group
  summarise(total_tree = sum(area, na.rm = T))

total_tree_2 <- terra::extract(tree, sa1_split_2_sv) %>% 
  group_by(ID) %>% #needs to be ID as extract creates an ID value from which to group
  summarise(total_tree = sum(area, na.rm = T))

# amend ID numbers of total_tree_2 as these re-start at 1 instead of 4541
total_tree_2 <- total_tree_2 %>% 
  mutate(ID = ID + nrow(sa1)/2)

# combine the outputs

total_tree <- rbind(total_tree_1, total_tree_2)

# calculate proportion of tree cover for each sa1 area

prop_tree <- total_tree %>%  
  mutate(prop_tree = total_tree/as.numeric(sa1_area$sa1_area)) %>% 
  select(ID, prop_tree)

## allocate sa1 identified to each row of and prop_tree object

# update sa1 object to include ID numbers
# these need to be the row numbers as this is how the extract() function for and tree calculate ID

sa1 <- sa1 %>% 
  mutate(ID = row_number())

# join sa1 object to and prop_tree objects

prop_tree <- prop_tree %>% 
  left_join(sa1, by = "ID") %>% 
  select(!geometry)

###### Calculate distance to water from survey points ######

### Create unique spatial unit for water calculations
# spatial unit should be the individual survey point

##convert edb_raw back to sf object ##

ebd_raw <- st_sf(ebd_raw)

## allocate neighbourhood data to ebd data

ebd1 <- st_join(ebd_raw, sa1, join = st_intersects)

## group by survey ID, SA1 code and geometry to ensure only 1 point per survey

ebd1 <- ebd1 %>%  
  group_by(checklist_id, sa1Cod, geometry) %>% 
  summarise() %>% 
  ungroup()

# make sure sa1Cod is in character form for later joining

ebd1 <- ebd1 %>% 
  mutate(sa1Cod = as.character(sa1Cod))

### prepare water location data

lake <- st_read("inputData/gis_osm_water_a_free_1.shp") %>% 
  st_transform(crs = 28355)

bay <-  st_read("inputData/portPhilipBay.shp") %>% 
  st_transform(crs = 28355) %>% 
  select(geometry)

creek <- st_read("inputData/gis_osm_waterways_free_1.shp") %>%
  st_transform(crs = 28355)

## clip to case study area

# add buffer to include water bodies that neighbour case study area

sa1_estUrb_buffer <- sa1_estUrb_clip %>% 
  st_buffer(200)

# run clipping function

lake <- st_intersection(lake, sa1_estUrb_buffer)

creek <- st_intersection(creek, sa1_estUrb_buffer)

## filter to only include appropriate water sources

lake <- lake %>% 
  filter(fclass %in% c("water", "wetland")) %>% 
  st_union() %>% 
  st_make_valid() %>% 
  st_sf()

creek <- creek %>%
  filter(fclass %in% c("stream", "river")) %>%
  st_union() %>%
  st_make_valid() %>%
  st_sf()

## combine water objects

water <- rbind(lake, creek, bay) %>% 
  st_union() %>% 
  st_make_valid() %>% 
  st_sf()

### measure distance between water and survey points

dist_water <- ebd1 %>% 
  mutate(distWater = as.vector(st_distance(., water)))

# convert into a data frame to remove geometry and sa1 vector
# object will be joined to master df via checklist_id. inclusion of Sa1Cod
# messes with left_join function later in code.

dist_water <- dist_water %>% 
  as.data.frame() %>% 
  select(!c(geometry, sa1Cod))

###### calculate dwelling density per spatial unit ######

### Create object of the number of dwellings per neighbourhood

## load abs data

dwel <- read.csv("inputData/2021Census_G36_VIC_SA1.csv") %>% 
  rename(sa1Cod = "SA1_CODE_2021") %>% 
  # change to character format to be consistant with format of sa1Cod in other objects 
  mutate(sa1Cod = as.character(sa1Cod))

## remove irrelevant dwelling types
# ABS data includes various types of dwellings including caravans and tents. These don't relate to perminant built structures and therefore need to be removed

dwel <- dwel %>% 
  mutate(total_dwel = Total_PDs_Dwellings
         # dwelling structure unknown
         - OPDs_Dwlling_structur_NS_Dwgs
         # "other dwellings"
         - OPDs_Other_dwelling_Tot_Dwgs
         # improvised home, tent, sleepers out etc...
         - OPDs_Ot_dwg_Im_hm_tnt_SO_Ds
         # cabin or houseboat
         - OPDs_Oth_dwg_cab_hboat_Ds
         # caravan
         - OPDs_Oth_dwg_Cvn_Ds) %>% 
  select(sa1Cod, total_dwel)

### calculate density of dwellings

## combine abs data with neighbourhood area data frames

dens_dwel <- dwel %>% 
  right_join(sa1_area, by = "sa1Cod")

## calculate density

dens_dwel <- dens_dwel %>% 
  mutate(den_dwel = total_dwel/as.numeric(sa1_area/10000)) %>% 
  select(sa1Cod, den_dwel)


################################################################################

# prepare master df of all variables

################################################################################

### combine objects into a single df

df <- left_join(ebd, prop_tree, by = "sa1Cod") %>%
  left_join(dens_dwel, by = "sa1Cod") %>% 
  left_join(dist_water, by = "checklist_id") %>% 
  mutate(land_st = factor(.$land_st, levels = c("ver", "brwn", "util", "rail", "otherIgs", "rpgs", "otherFgs", "nonGs")),
         common_name = factor(.$common_name, levels = c(unique(ebd$common_name))))

### scale the covariate data to make it easier to compare across your different covariates
# scaling function divides by mean and standard deviation --> this is just what needs to be done to scale

df$prop_tree <- scale(df$prop_tree) 
df$den_dwel <- scale(df$den_dwel)
df$distWater <- scale(df$distWater)

### re-level land_t to make it is easier for the model to pick a level to act as the intercept term (baseline from which everything gets compared)

df$land_t <- factor(df$land_t, levels = c("IGS","formal", "nonGs"))

################################################################################

# build GLMM models

################################################################################

m1 <- glmmTMB(presAbs ~ land_t * frg_gld +
                (1 | land_st) +
                (1 | common_name) +
                prop_tree +
                den_dwel +
                distWater,
              data = df,
              offset = log(drtn_mn),
              family = binomial(link = "cloglog"))

m2 <- glmmTMB(presAbs ~ land_t * cnsvt_ +
                (1 | land_st) +
                (1 | common_name)+
                prop_tree +
                den_dwel +
                distWater,
              data = df,
              offset = log(drtn_mn),
              family = binomial(link = "cloglog"))

m3 <- glmmTMB(presAbs ~ land_t * frg_dmn +
                (1 | land_st) +
                (1 | common_name)+
                prop_tree +
                den_dwel +
                distWater,
              data = df,
              offset = log(drtn_mn),
              family = binomial(link = "cloglog"))

################################################################################

# create graphs from GLMM models
# example for foraging guild provided. Other graphs replicate code, chaning
# underlying model in response to avenue of enquiry. 

################################################################################

###### graph foraging guild by land type ######

### generate predictions to graph using Monte Caro simulation ###

# Process simulates new log(odds) values based on the predicted log(odds) from
# the model. The process simulations many (1000s) log(odds) values and then converts
# these into probabilities. This is more accurate than trying to convert the predicted
# log(odds) directly back into probability scale, which has biasing issues (for maths
# reasons).

## define the number of Monte Carlo simulations to run to get the estimate and CIs

n_sims <- 10000

## create new data frame with predictions in it

new_df <- expand.grid(
  # values need to be organised with the same levels as the model so the output
  # of the prediction() function can be aligned to the correct land_st.
  land_t = levels(df$land_t),
  common_name = levels(df$common_name),
  # include land_st as na. land_st needs to be in the model but this script is
  # focusing on land_t specifically.
  land_st = NA,
  # covariate values held at either mean or at 0. This says "if we hold the
  # covariates constant, what is the impact of the remaining, non-constant variable"
  prop_tree = mean(df$prop_tree, na.rm = TRUE),
  den_dwel = mean(df$den_dwel, na.rm = TRUE),
  distWater = mean(df$distWater, na.rm = TRUE),
  # keep offset term constant.
  drtn_mn = mean(df$drtn_mn, na.rm = TRUE)
)

## add frg_gld but ensure it perfectly aligns with common_name

# create vector of names

fg_names <- df %>% 
  group_by(common_name, frg_gld) %>% 
  summarise(.groups = "drop")

# run join function to attribute frg_gld names based on common_name

new_df <- new_df %>% 
  left_join(fg_names, by = "common_name")

## predict values from model
# do this prior to simulating many dummy rows of data as this is much faster computationally

pred <- predict(m1,
                newdata = new_df,
                type = "link",
                se.fit = TRUE)

## attribute predicted log(odds) for each species x land use combination

new_df <- new_df %>%
  # for each land_st and common_name, get the mean and SE on the log(odds) scale
  mutate(
    mu = pred$fit,
    se = pred$se.fit,
  )

## prepare dummy rows to give simulated log(odds) values

new_df <- new_df %>%
  # duplicates all values (n() gives the size of new_df)
  slice(rep(1:n(), times = n_sims)) %>%
  # gives a number reflecting the which duplication each row is part of
  mutate(sim = rep(1:n_sims, each = nrow(new_df)))

## simulate log(odds) values for each dummy row (log(odds) for observing each
## species on each land use sub-type)

# the simulated values are based on a normal distribution around the mean and se
# from the model predicted log(odds). simulated values are very similar to the
# modeled values. However, the simulated values take into consideration the error
# that is reflected in the model predicted se values (effectively, it combines
# the mu value and se value into a new, single number). The error is needed as
# uncertainty is a core requirement of the Monte Carlo simulation; to use the mu
# value directly risks biasing the process. We get a better approximation of the
# probability if we incorporate the error.

new_df <- new_df %>% 
  mutate(
    # create a random number for each row based on a normal distribution around
    # a mean of "mu" and a standard error of "se".
    eta = rnorm(n(), mu, se))

## convert to probability of observing each species on each land use sub-type
## for each dummy row

# convert through cloglog to a probability of observing that species on that
# site, in random observation.

new_df <- new_df %>%
  mutate(
    probability = 1 - exp(-exp(eta))
  )

## aggregate to calculate probability of observing any species in each frg_gld
# on each land use type for each simulation

# the probability of observing 1 or more species = one minus the probability of
# observing none of the species

new_df <- new_df %>%
  group_by(
    frg_gld,
    land_t,
    sim
  ) %>%
  summarise(
    probability_any = 1 - prod(1 - probability)
  )

## find the mean and 95%CI by amalgamating all the simulations

new_df <- new_df %>%
  # now we have n_sims simulations of the probability of observing any species,
  # for each land cover type, we can summarise them for plotting with a
  # prediction and 95% confidence interval
  group_by(
    frg_gld,
    land_t
  ) %>%
  summarise(
    probability_mean = mean(probability_any),
    probability_lower = quantile(probability_any, 0.025),
    probability_upper = quantile(probability_any, 0.975)
  )

###### prep graphing aesthetics: foraging guild by land type ######

## prep colours palette

# create an object for each colour
colour1 <- "#4335A7"
colour2 <- "#80C4E9"
colour3 <- "#ffd1adff"

## prep labels

# change land type names

new_df <- new_df %>% 
  # convert to factor to ensure that leveling stays as desired
  mutate(land_t = factor(case_when(land_t == "IGS" ~ "IGS",
                                   land_t == "formal" ~ "Formal GS",
                                   land_t == "nonGs" ~ "Not GS"),
                         levels = c("IGS", "Formal GS", "Not GS")))

###### graph data: Foraging guild by land type ######

ggplot(new_df, aes(x = land_t, y = probability_mean, color = land_t)) +
  geom_point(size = 2) +
  geom_pointrange(aes(ymin = probability_lower, ymax = probability_upper), size = 0, linewidth = 0.8)+
  facet_wrap(vars(frg_gld), ncol = 5)+
  scale_color_manual(labels = c(levels(new_df$land_t)), values = c(colour1, colour2, colour3))+
  labs(x = "Land Type", y = "Probability of presence (%)", color = "Land Type") +
  theme_minimal()+
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "mm")),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "mm")),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))




###### graph foraging guild by land sub-type ######

### generate predictions to graph using Monte Caro simulation ###

# Process simulates new log(odds) values based on the predicted log(odds) from
# the model. The process simulations many (1000s) log(odds) values and then converts
# these into probabilities. This is more accurate than trying to convert the predicted
# log(odds) directly back into probability scale, which has biasing issues (for maths
# reasons).

## define the number of Monte Carlo simulations to run to get the estimate and CIs

n_sims <- 10000

## create new data frame with predictions in it

new_df <- expand.grid(
  # values need to be organised with the same levels as the model so the output
  # of the prediction() function can be aligned to the correct land_st.
  land_st = levels(df$land_st),
  common_name = levels(df$common_name),
  # covariate values held at either mean or at 0. This says "if we hold the
  # covariates constant, what is the impact of the remaining, non-constant variable"
  prop_tree = mean(df$prop_tree, na.rm = TRUE),
  den_dwel = mean(df$den_dwel, na.rm = TRUE),
  distWater = mean(df$distWater, na.rm = TRUE),
  # keep offset term constant.
  drtn_mn = mean(df$drtn_mn, na.rm = TRUE)
)

## add frg_gld but ensure it perfectly aligns with common_name

# create vector of names

fg_names <- df %>% 
  group_by(common_name, frg_gld) %>% 
  summarise(.groups = "drop")

# run join function to attribute frg_gld names based on common_name

new_df <- new_df %>% 
  left_join(fg_names, by = "common_name")

## add land_t but ensure it perfectly aligns with land_st

lt_names <- df %>% 
  group_by(land_st, land_t) %>% 
  summarise(.groups = "drop")

# run join function to attribute land_t names based on land_st

new_df <- new_df %>% 
  left_join(lt_names, by = "land_st")

## predict values from model
# do this prior to simulating many dummy rows of data as this is much faster computationally

pred <- predict(m1,
                newdata = new_df,
                type = "link",
                se.fit = TRUE)

## attribute predicted log(odds) for each species x land use combination

new_df <- new_df %>%
  # for each land_st and common_name, get the mean and SE on the log(odds) scale
  mutate(
    mu = pred$fit,
    se = pred$se.fit,
  )

## prepare dummy rows to give simulated log(odds) values

new_df <- new_df %>%
  # duplicates all values (n() gives the size of new_df)
  slice(rep(1:n(), times = n_sims)) %>%
  # gives a number reflecting the which duplication each row is part of
  mutate(sim = rep(1:n_sims, each = nrow(new_df)))

## simulate log(odds) values for each dummy row (log(odds) for observing each
## species on each land use sub-type)

# the simulated values are based on a normal distribution around the mean and se
# from the model predicted log(odds). simulated values are very similar to the
# modeled values. However, the simulated values take into consideration the error
# that is reflected in the model predicted se values (effectively, it combines
# the mu value and se value into a new, single number). The error is needed as
# uncertainty is a core requirement of the Monte Carlo simulation; to use the mu
# value directly risks biasing the process. We get a better approximation of the
# probability if we incorporate the error.

new_df <- new_df %>% 
  mutate(
    # create a random number for each row based on a normal distribution around
    # a mean of "mu" and a standard error of "se".
    eta = rnorm(n(), mu, se))

## convert to probability of observing each species on each land use sub-type
## for each dummy row

# convert through cloglog to a probability of observing that species on that
# site, in random observation.

new_df <- new_df %>%
  mutate(
    probability = 1 - exp(-exp(eta))
  )

## aggregate to calculate probability of observing any species in each frg_gld
# on each land use type for each simulation

# the probability of observing 1 or more species = one minus the probability of
# observing none of the species

new_df <- new_df %>%
  group_by(
    frg_gld,
    land_st,
    sim
  ) %>%
  summarise(
    probability_any = 1 - prod(1 - probability)
  )

## find the mean and 95%CI by amalgamating all the simulations

new_df <- new_df %>%
  # now we have n_sims simulations of the probability of observing any species,
  # for each land cover type, we can summarise them for plotting with a
  # prediction and 95% confidence interval
  group_by(
    frg_gld,
    land_st
  ) %>%
  summarise(
    probability_mean = mean(probability_any),
    probability_lower = quantile(probability_any, 0.025),
    probability_upper = quantile(probability_any, 0.975)
  )

###### prep graphing aesthetics: foraging guild by land sub-type ######

## prep labels

# run join function to attribute land_t names based on land_st

new_df <- new_df %>% 
  left_join(lt_names, by = "land_st")

# change land type names

new_df <- new_df %>% 
  # convert to factor to ensure that leveling stays as desired
  mutate(land_t = factor(case_when(land_t == "IGS" ~ "IGS",
                                   land_t == "formal" ~ "Formal GS",
                                   land_t == "nonGs" ~ "Not GS"),
                         levels = c("IGS", "Formal GS", "Not GS")))

## change land subtype names

new_df <- new_df %>% 
  mutate(land_st = factor(case_when(land_st == "ver" ~ "Verge",
                                    land_st == "brwn" ~ "Brownfield",
                                    land_st == "util" ~ "Utility",
                                    land_st == "rail" ~ "Railway",
                                    land_st == "otherIgs" ~ "Other IGS",
                                    land_st == "otherFgs" ~ "Garden areas",
                                    land_st == "rpgs" ~ "RPGS",
                                    land_st == "nonGs" ~ "Not GS"),
                          levels = c("Verge",
                                     "Brownfield",
                                     "Utility",
                                     "Railway",
                                     "Other IGS",
                                     "Garden areas",
                                     "RPGS",
                                     "Not GS")))

## filter for selected graphs to be highlighted"

frug <- new_df %>% 
  filter(frg_gld == "Frugivore")

gran <- new_df %>% 
  filter(frg_gld == "Granivore")

nect <- new_df %>% 
  filter(frg_gld == "Nectarivore")

omni <- new_df %>% 
  filter(frg_gld == "Omnivore")

###### graph data: foraging guild by land_sub-type ######

### create a graph for each targetted foraging guild

## frugivore graph

ggplot(frug, aes(x = land_st, y = probability_mean, color = land_t)) +
  geom_point(size = 2) +
  geom_pointrange(aes(ymin = probability_lower, ymax = probability_upper), size = 0.5, linewidth = 0.8)+
  facet_wrap(vars(frg_gld), ncol = 5)+
  scale_color_manual(labels = c(levels(new_df$land_t)), values = c(colour1, colour2, colour3))+
  labs(x = "Land Sub-type", y = "Probability of presence (%)", color = "Land Type") +
  ylim(0,0.4)+
  theme_minimal()+
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "mm")),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "mm")),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

## granivore graph

ggplot(gran, aes(x = land_st, y = probability_mean, color = land_t)) +
  geom_point(size = 2) +
  geom_pointrange(aes(ymin = probability_lower, ymax = probability_upper), size = 0.5, linewidth = 0.8)+
  facet_wrap(vars(frg_gld), ncol = 5)+
  scale_color_manual(labels = c(levels(new_df$land_t)), values = c(colour1, colour2, colour3))+
  labs(x = "Land Sub-type", y = "Probability of presence (%)", color = "Land Type") +
  ylim(0.5,.9)+
  theme_minimal()+
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "mm")),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "mm")),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

## nectarivore graph

ggplot(nect, aes(x = land_st, y = probability_mean, color = land_t)) +
  geom_point(size = 2) +
  geom_pointrange(aes(ymin = probability_lower, ymax = probability_upper), size = 0.5, linewidth = 0.8)+
  facet_wrap(vars(frg_gld), ncol = 5)+
  scale_color_manual(labels = c(levels(new_df$land_t)), values = c(colour1, colour2, colour3))+
  labs(x = "Land Sub-type", y = "Probability of presence (%)", color = "Land Type") +
  ylim(0.6,1)+
  theme_minimal()+
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "mm")),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "mm")),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

## omnivore graph

ggplot(omni, aes(x = land_st, y = probability_mean, color = land_t)) +
  geom_point(size = 2) +
  geom_pointrange(aes(ymin = probability_lower, ymax = probability_upper), size = 0.5, linewidth = 0.8)+
  facet_wrap(vars(frg_gld), ncol = 5)+
  scale_color_manual(labels = c(levels(new_df$land_t)), values = c(colour1, colour2, colour3))+
  labs(x = "Land Sub-type", y = "Probability of presence (%)", color = "Land Type") +
  ylim(0.6,1)+
  theme_minimal()+
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "mm")),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "mm")),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))