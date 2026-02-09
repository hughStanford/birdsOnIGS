library(pacman)

p_load(tidyverse, glmmTMB, arrow, DHARMa, sf, mgcv, janitor)


# Prepare data sample -----------------------------------------------------

infile <- "spatial_autocorrelation_testing/"

df <- st_read(
  file.path(infile,
            "master_df.gpkg")) 

# Sample 250 surveys ensuring spatial diversity while preserving eBird's 
# natural location duplication patterns. Iteratively resamples to accomodate 
# duplicate geometries until 250 unique locations are achieved (without removing surveys with duplicated geometries).

set.seed(123)
unique_surveys <- unique(df$chckls_)
n_surveys <- 250

sampled_surveys <- sample(unique_surveys, n_surveys)
df_subset <- df %>% 
  filter(chckls_ %in% sampled_surveys)

unique_geom_count <- df_subset %>%
  group_by(geom) %>%
  slice(1) %>%
  ungroup() %>%
  nrow()

max_iterations <- 100
iteration <- 0

while(unique_geom_count < n_surveys && iteration < max_iterations) {
  iteration <- iteration + 1
  
  duplicate_count <- n_surveys - unique_geom_count
  
  available_surveys <- setdiff(unique_surveys, sampled_surveys)
  
  if(length(available_surveys) < duplicate_count) {
    warning("Not enough unique surveys remaining to reach target number of unique geometries")
    break
  }
  
  new_sampled_surveys <- sample(available_surveys, duplicate_count)
  
  new_df_subset <- df %>% 
    filter(chckls_ %in% new_sampled_surveys)
  
  df_subset <- rbind(df_subset, new_df_subset)
  
  unique_geom_count <- df_subset %>%
    group_by(geom) %>%
    slice(1) %>%
    ungroup() %>%
    nrow()
  
  cat("Iteration:", iteration, "| Unique geometries:", unique_geom_count, "\n")
}

if(iteration >= max_iterations) {
  warning("Maximum iterations reached without achieving target number of unique geometries")
}

coords <- st_coordinates(df_subset)
df_subset$x_coord <- coords[, "X"]
df_subset$y_coord <- coords[, "Y"]

# 2. build model ---------------------------------------------------------

# m <- glmmTMB(presAbs ~ land_t * frg_gld +
#                (1 | land_st) +
#                (1 | cmmn_nm) +
#                prop_tree +
#                den_dwel +
#                distWater,
#              data = df_subset,
#              offset = log(drtn_mn),
#              family = binomial(link = "cloglog"))

m <- glmmTMB(presAbs ~ land_t * cnsvt_ +
                (1 | land_st) +
                (1 | cmmn_nm) +
                prop_tree +
                den_dwel +
                distWater,
              data = df_subset,
              offset = log(drtn_mn),
              family = binomial(link = "cloglog"))

# m <- glmmTMB(presAbs ~ land_t * frg_dmn +
#                 (1 | land_st) +
#                 (1 | cmmn_nm)+
#                 prop_tree +
#                 den_dwel +
#                 distWater,
#               data = df_subset,
#               offset = log(drtn_mn),
#               family = binomial(link = "cloglog"))

# 3. Run spatial autocorrelation tests ------------------------------------

# Create simulated residuals and group by unique survey ID (chckls_)
sim_residuals <- simulateResiduals(fittedModel = m, n = 500)

sim_residuals_agg <- recalculateResiduals(sim_residuals,
                                          group = df_subset$chckls_)

# Ensure unique coordinates per survey for testSpatialAutocorrelation() 
# which requires distinct geometries.
coords_unique <- df_subset %>% 
  group_by(chckls_) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(
    x_jitter = x_coord + runif(n(), -1, 1),
    y_jitter = y_coord + runif(n(), -1, 1)
  ) %>% 
  arrange(chckls_)

# Test for spatial autocorrelation
spatial_test <- testSpatialAutocorrelation(sim_residuals_agg,
                                           x = coords_unique$x_jitter,
                                           y = coords_unique$y_jitter)

spatial_test
