## ----load-packages, message=FALSE, warning=FALSE, results="hold"-------------------------
library(dplyr)
library(reticulate)
library(sf)
use_condaenv("Space")

## ----set-seed----------------------------------------------------------------------------
set.seed(24601)


## ----load-data---------------------------------------------------------------------------
dat_S <- readRDS("data/dat_S.rds")
dat_A <- readRDS("data/dat_A.rds")


## ----source-python-----------------------------------------------------------------------
source_python("dist_between_vectors.py")


## ----dist-AS-----------------------------------------------------------------------------
dat_D <- dist_between(as.matrix(dat_S %>% select(s_id,latitude,longitude)),
                      as.matrix(dat_A %>% select(a_id,latitude,longitude)))
colnames(dat_D) <- c("s_id","a_id","dist_m")
dat_D <- as_tibble(dat_D) %>%
  mutate(s_id = as.integer(s_id),
         a_id = as.integer(a_id),
         dist_m = as.integer(dist_m),
         dist_km = dist_m / 1000,
         dist_mi = dist_m*0.62137119223733/1000)
saveRDS(dat_D,"data/dat_D.rds")


## ----dist-SS-----------------------------------------------------------------------------
dat_DS <- dist_between(as.matrix(dat_S %>% select(s_id,latitude,longitude)),
                       as.matrix(dat_S %>% select(s_id,latitude,longitude)))
colnames(dat_DS) <- c("s1_id","s2_id","dist_m")
dat_DS <- as_tibble(dat_DS) %>%
  mutate(s1_id = as.integer(s1_id),
         s2_id = as.integer(s2_id),
         dist_m = as.integer(dist_m),
         dist_km = dist_m / 1000,
         dist_mi = dist_m*0.62137119223733/1000) %>%
  filter(s1_id != s2_id)
saveRDS(dat_DS,"data/dat_DS.rds")


## ----pick-isolated-----------------------------------------------------------------------
dat_S_isolated <- dat_D %>%
  # find distance to nearest grocery store
  group_by(a_id) %>% 
  summarize(dist_mi = min(dist_mi)) %>%
  # keep in the relevant range
  filter(between(dist_mi,0.2,2))


## ----move-random-------------------------------------------------------------------------
dat_S_isolated_random <- dat_S_isolated %>% 
  inner_join(dat_A, by="a_id") %>% 
  # randomly shift by approximately 0.025 - 0.05 miles to move away from business
  # flip needed due to mean shift!
  mutate(random_shock_lat = rnorm(n=n(),mean=0.0004,sd=0.0001),
         random_shock_lon = rnorm(n=n(),mean=0.0004,sd=0.0001),
         sign_flip_lat = (runif(n=n()) > 0.5)*2-1,
         sign_flip_lon = (runif(n=n()) > 0.5)*2-1
  ) %>% 
  mutate(latitude=latitude+random_shock_lat*sign_flip_lat,
         longitude=longitude+random_shock_lon*sign_flip_lon) %>% 
  mutate(s_id = as.integer(1000 + row_number())) %>% 
  select(s_id, latitude, longitude)


## ----skip-close--------------------------------------------------------------------------
tmp <- dist_between(as.matrix(dat_S_isolated_random),as.matrix(dat_S_isolated_random))
colnames(tmp) <- c("s1_id","s2_id","dist_m")
tmp <- as_tibble(tmp) %>% 
  filter(s2_id < s1_id) %>% 
  filter(dist_m < 100) %>% 
  arrange(s2_id) %>% 
  .$s2_id
dat_S_candidate <- dat_S_isolated_random %>% filter(!(s_id %in% tmp))


## ----save-candidates---------------------------------------------------------------------
saveRDS(dat_S_candidate,"data/dat_S_candidate_random.rds")


## ----calc-2d-proj------------------------------------------------------------------------
# use NAD83(2011) projection, EPSG:6419 for California 3 zone (~ Bay Area)
mat_S_xy <- cbind(dat_S$s_id,
                  sf_project(from="WGS84", to="EPSG:6419",
                             pts=cbind(dat_S$longitude,dat_S$latitude)))
colnames(mat_S_xy) <- c("s_id","x","y")
mat_S_candidate_xy <- cbind(dat_S_candidate$s_id,
                            sf_project(from="WGS84", to="EPSG:6419",
                                       pts=cbind(dat_S_candidate$longitude,
                                                 dat_S_candidate$latitude)))
colnames(mat_S_candidate_xy) <- c("s_id","x","y")
mat_A_xy <- cbind(dat_A$a_id,
                  sf_project(from="WGS84", to="EPSG:6419",
                             pts=cbind(dat_A$longitude,dat_A$latitude)))
colnames(mat_A_xy) <- c("a_id","x","y")
# convert to tables instead of matrices for easier handling
dat_S_xy <- as_tibble(mat_S_xy)
dat_S_candidate_xy <- as_tibble(mat_S_candidate_xy)
dat_A_xy <- as_tibble(mat_A_xy)


## ----define-nearby-----------------------------------------------------------------------
dist_keep_mi <-  2 * sqrt(2)  # * sqrt(2) fills a square with base 4 miles instead of a circle with diameter 4


## ----all-near-real-----------------------------------------------------------------------
dat_S_A <- rbind(
  # other businesses
  dat_D %>% 
    filter(dist_mi < dist_keep_mi) %>% 
    select(s_id,a_id) %>%
    # grid position of the grocery store
    inner_join(dat_S_xy, by="s_id") %>% 
    rename(x_s = x, y_s = y) %>% 
    # grid position of other business
    inner_join(dat_A_xy, by="a_id") %>%
    # relative positions
    mutate(x = x-x_s, y = y-y_s) %>% 
    # industry of other business
    inner_join(dat_A %>% select(a_id,naics_code), by="a_id") %>% 
    select(s_id, a_id, x, y, naics_code),
  # other treatment locations
  dat_DS %>% 
    filter(dist_mi < dist_keep_mi) %>% 
    select(s1_id,s2_id) %>%
    # grid position of the grocery store
    inner_join(dat_S_xy, by=c("s1_id"="s_id")) %>% 
    rename(x_s = x, y_s = y) %>% 
    # grid position of other business
    inner_join(dat_S_xy, by=c("s2_id"="s_id")) %>%
    # relative positions
    mutate(x = x-x_s, y = y-y_s) %>% 
    # industry of other business
    inner_join(dat_S %>% select(s_id,naics_code), by=c("s2_id"="s_id")) %>% 
    mutate(s_id = s1_id, a_id = -s2_id) %>% 
    select(s_id, a_id, x, y, naics_code) %>%
    arrange(s_id)
) %>% 
  arrange(s_id) %>% 
  mutate(x = round(x,2),
         y = round(y,2))


## ----distances-to-random-----------------------------------------------------------------
# distance to other businesses
dat_D_candidate <- dist_between(as.matrix(dat_S_candidate),
                                as.matrix(dat_A %>% select(a_id,latitude,longitude)))
colnames(dat_D_candidate) <- c("s_id","a_id","dist_m")
dat_D_candidate <- as_tibble(dat_D_candidate) %>% 
  mutate(s_id = as.integer(s_id),
         a_id = as.integer(a_id),
         dist_mi = dist_m*0.62137119223733/1000)
# distance to grocery stores
dat_DS_candidate <- dist_between(as.matrix(dat_S_candidate),
                                 as.matrix(dat_S %>% select(s_id,latitude,longitude)))
colnames(dat_DS_candidate) <- c("s1_id","s2_id","dist_m")
dat_DS_candidate <- as_tibble(dat_DS_candidate) %>% 
  mutate(s_id = as.integer(s1_id),
         a_id = as.integer(s2_id),
         dist_mi = dist_m*0.62137119223733/1000)


## ----all-near-random---------------------------------------------------------------------
dat_S_candidate_A <- rbind(
  # other businesses
  dat_D_candidate %>% 
    filter(dist_mi < dist_keep_mi) %>% 
    select(s_id,a_id) %>%
    # grid position of the grocery store
    inner_join(dat_S_candidate_xy, by="s_id") %>% 
    rename(x_s = x, y_s = y) %>% 
    # grid position of other business
    inner_join(dat_A_xy, by="a_id") %>%
    # relative positions
    mutate(x = x-x_s, y = y-y_s) %>% 
    # industry of other business
    inner_join(dat_A %>% select(a_id,naics_code), by="a_id") %>% 
    select(s_id, a_id, x, y, naics_code),
  # real treatment locations
  dat_DS_candidate %>% 
    filter(dist_mi < dist_keep_mi) %>% 
    select(s1_id,s2_id) %>%
    # grid position of the grocery store
    inner_join(dat_S_candidate_xy, by=c("s1_id"="s_id")) %>% 
    rename(x_s = x, y_s = y) %>% 
    # grid position of other business
    inner_join(dat_S_xy, by=c("s2_id"="s_id")) %>%
    # relative positions
    mutate(x = x-x_s, y = y-y_s) %>% 
    # industry of other business
    inner_join(dat_S %>% select(s_id,naics_code), by=c("s2_id"="s_id")) %>% 
    mutate(s_id = s1_id, a_id = -s2_id) %>% 
    select(s_id, a_id, x, y, naics_code) %>%
    arrange(s_id)
) %>% 
  arrange(s_id) %>% 
  mutate(x = round(x,2),
         y = round(y,2))


## ----save-csv----------------------------------------------------------------------------
write.csv(x = dat_S_A %>% 
            filter(floor(naics_code/100) != 4451) %>% 
            select(s_id,a_id,x,y,naics_code),
          file = "neural-net/grid_S_I.csv",
          row.names=FALSE, quote=FALSE)
write.csv(x = dat_S_A %>% 
            filter(floor(naics_code/100) == 4451) %>% 
            select(s_id,a_id,x,y),
          file = "neural-net/grid_S_S.csv",
          row.names=FALSE, quote=FALSE)
write.csv(x = dat_S_candidate_A %>% 
            filter(floor(naics_code/100) != 4451) %>% 
            select(s_id,a_id,x,y,naics_code),
          file = "neural-net/grid_S_random_I.csv",
          row.names=FALSE, quote=FALSE)
write.csv(x = dat_S_candidate_A %>% 
            filter(floor(naics_code/100) == 4451) %>% 
            select(s_id,a_id,x,y),
          file = "neural-net/grid_S_random_S.csv",
          row.names=FALSE, quote=FALSE)

