################################################################################
#### Identificar los estados de los países ocn los que limita Colombia
################################################################################
library(dplyr)
library(tidyverse)
library(mapsf)
library(purrr)
library(lubridate)
library(stringr)
library(parallel)
library(doSNOW)
library(raster)
library(sf)
library(stringr)

# Leer el archivo shapefile
shp_depto <- st_read("Departamentos202208_shp/Depto.shp")
shp_mundo <- st_read("admin00/admin00.shp") # quiet=TRUE
COLMB_ = shp_mundo$FIPS_ADMIN[str_detect(shp_mundo$FIPS_ADMIN, "CO")]
COLMB_2 <- shp_mundo %>% 
  filter(CNTRY_NAME %in% c("Colombia")) %>% 
  filter(FIPS_ADMIN %in% c(COLMB_))

mundocol <- shp_mundo %>% 
  filter(CNTRY_NAME %in% c("Colombia","Peru","Brazil","Venezuela","Ecuador","Panama")) %>% 
  filter(FIPS_ADMIN %in% c(COLMB_, "EC22", "PE16", "BR04", "VE01", "VE06", "VE03", "VE20", "VE23"))

# Exportar el objeto sf a un archivo shapefile
#st_write(mundocol, "admin00/admin00_vecinos.shp", append = FALSE)
st_write(COLMB_2, "admin00/admin00_colombia.shp", append = FALSE)

mundocol <- st_read("admin00/admin00_vecinos.shp") # quiet=TRUE

ggplot() + 
  geom_sf(data = COLMB_2) #+
  #geom_sf_text(aes(label = FIPS_ADMIN), data = mundocol, size = 2) 

  
# Functions to use --------------------------------------------------------
download_chirps <- function(year){
  
  print(year)
  dates <- seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")), by="months")
  #dates <- gsub('-', '.', as.character(dates))
  dates <- paste(year(dates), ifelse(month(dates)<10, paste0(0,month(dates)), month(dates)), sep ='.')
  paths <- paste0(pth, '/chirps-v2.0.', dates, '.tif.gz')
  
  lapply(1:length(paths), function(k){
    download.file(url = paths[k],
                  destfile = paste0('./chirps/', basename(paths[k])),
                  mode = 'wb')
  })
}    

# Load data ---------------------------------------------------------------
#pth <- 'ftp://ftp.chg.ucsb.edu/pub/org/chg/products/CHIRPS-2.0/global_daily/tifs/p05/'
#pth <- 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/tifs/p05/'
pth <- 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/tifs/'
dir.create('./chirps/', recursive = TRUE)

#yrs <- 1981:2022
yrs <- 1999:2022
map(.x = yrs, .f = download_chirps)

# https://stackoverflow.com/questions/65404992/decompress-tif-gz-file-using-r
chirps <- list.files('./chirps', full.names = TRUE, pattern = '.tif.gz$') 
dir.create('./chirps/stack_v2/', recursive = TRUE)
lapply(1:length(chirps), function(k){
  R.utils::gunzip(filename = chirps[k],
                  destname = paste0('./chirps/stack_v2/', basename(gsub('.gz', '',chirps[k]))),
                  remove = FALSE)
})

# Functions to use --------------------------------------------------------
extract_mask <- function(year){
  
  print(year)
  fle <- grep(year, fles, value = TRUE)
  
  print('To extract by mask')
  cl <- makeCluster(2)
  registerDoSNOW(cl)
  
  foreach(i = 1:length(fle), .verbose = TRUE) %dopar% {
    print(paste("mes", i))  
    #shpf <- shapefile('//dapadfs/workspace_cluster_9/Coffee_Cocoa2/_coffeeColombia/shp/base/MGN_DPTO_POLITICO.shp')
    #shpf <- shapefile('/content/drive/MyDrive/Curso/MGN_DPTO_POLITICO.shp')
    #shpf <- sf::st_read("admin00/admin00_vecinos.shp")
    shpf <- sf::st_read("admin00/admin00_colombia.shp")
    rst <- raster::raster(fle[i])
    rst <- raster::crop(rst, shpf)
    rst <- raster::mask(rst, shpf)
    rst <- rst * 1
    out <- 'chirps/colombia/'
    #out <- 'chirps/colombia_vecinos/'
    raster::writeRaster(rst, filename = paste0(out, basename(fle[i])), overwrite = TRUE)
    
  }
  
  stopCluster(cl)
  print(paste0('Done ', year))  
  
}

# Load data ---------------------------------------------------------------
root <- './chirps/stack_v2'
fles <- list.files(root, full.names = TRUE, pattern = '.tif$')

# Get the years available--------------------------------------------------
years <- unique(str_sub(basename(fles), start = 13, end = 16))

#dir.create('chirps/colombia_vecinos/', recursive = TRUE)
dir.create('chirps/colombia/', recursive = TRUE)

# To extract by mask for only Colombia ------------------------------------

for(j in 1:length(years)){
  print(years[j])
  extract_mask(year = years[j])
  print(' Done')
}

#chirps <- list.files('chirps/colombia_vecinos/', full.names = TRUE, pattern = '.tif$') 
chirps <- list.files('chirps/colombia/', full.names = TRUE, pattern = '.tif$') 
tail(chirps)

chirps <- terra::rast(chirps)
chirps

chirps_df <- as.data.frame(chirps, xy = TRUE) # this works with Raster* objects as well

#arrow::write_parquet(chirps_df, "chirps_df.parquet")
arrow::write_parquet(chirps_df, "chirps_colombia_df.parquet")
