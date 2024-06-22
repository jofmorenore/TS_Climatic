library(arrow)
library(dplyr)
library(sf)
library(terra)
library(openxlsx)
library(geodata)
library(ggplot2)
library(ggtext) 
library(tidyverse)

setwd('TS_climatic/')
list.files()

# Carga los datos ----
#sttns = arrow::read_parquet("sttns_pcp_col.parquet")
#dim(sttns)

pcp_col = openxlsx::read.xlsx("pcp_col.xlsx")
dim(pcp_col)
#class(pcp_col)

chirps_df = arrow::read_parquet("chirps_df.parquet")

# Carga datos de elevación de datos abiertos ----
SRTM_30 <- rast("Servicio-159/SRTM30/SRTM_30_Col1.tif")
print(SRTM_30)

SRTM_30 <- project(SRTM_30, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
crs(SRTM_30) 

# Convierte a objeto sf ----
chirps.point <- st_as_sf(x = chirps_df, 
                         coords = c("x", "y"),
                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

sttns = dplyr::distinct(pcp_col, `X` , `Y`, .keep_all = TRUE)

sttns.points <- st_as_sf(x = sttns, 
                         coords = c("X", "Y"),
                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 

#terra::crs(sttns.points)
print(sttns.points)
# Extrae la altitud de los puntos de las estaciones ----
SRTM_30.sttns <- extract(SRTM_30, sttns.points)

SRTM_30.sttns$SRTM_30_Col1[1] = extract(geodata::elevation_3s(lon = pcp_col$X[1], lat = pcp_col$Y[1], path = './tmpr'), 
                                        cbind(pcp_col$X[1], pcp_col$Y[1]))[[1]]

SRTM_30.sttns$SRTM_30_Col1[2] = extract(geodata::elevation_3s(lon = pcp_col$X[2], lat = pcp_col$Y[2], path = './tmpr'), 
                                        cbind(pcp_col$X[2]-0.006, pcp_col$Y[2]))[[1]]

# Define los cuartiles para determinar n_chirps ----
Q1 = as.vector(quantile(SRTM_30.sttns$SRTM_30_Col1)[2])
Q2 = as.vector(quantile(SRTM_30.sttns$SRTM_30_Col1)[3])
Q3 = as.vector(quantile(SRTM_30.sttns$SRTM_30_Col1)[4])

SRTM_30.sttns$n_chirps = ifelse(SRTM_30.sttns$SRTM_30_Col1 <= Q1,4,
                                ifelse(SRTM_30.sttns$SRTM_30_Col1 > Q1 & SRTM_30.sttns$SRTM_30_Col1 <= Q2, 3,
                                       ifelse(SRTM_30.sttns$SRTM_30_Col1 > Q2 & SRTM_30.sttns$SRTM_30_Col1 <= Q3,2,
                                              ifelse(SRTM_30.sttns$SRTM_30_Col1 > Q3, 1, NA))))


table(SRTM_30.sttns$n_chirps)
sttns.points$n_chirps = SRTM_30.sttns$n_chirps

# Función para extraer los pixeles chirps más cercanos ----
#chirps.point#[1]

pol_ = function(x, r){
  #set.seed(25052024)
  #x = sample(1:length(sttns.points$ID), size = 1)
  #print(x)
  pt = sttns.points[x,]
  pol = data.frame(matrix(c(pt$geometry[[1]] +c(-r,-r),
                            pt$geometry[[1]] +c(-r,r),
                            pt$geometry[[1]] +c(r,r),
                            pt$geometry[[1]] +c(r,-r)
                            ), nc = 2, byrow = F
                          )
                   )
  names(pol) = c("X", "Y")
  pol = pol %>%
    st_as_sf(coords = c("X", "Y"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON")                            
  join_chirps_pol = st_join(pol, chirps.point)
  chirps_ = st_join(chirps.point, join_chirps_pol)
  chirps_ = na.omit(chirps_) %>% distinct(geometry, .keep_all = T)
  distances <- st_distance(pt, chirps_)
  s_distances <- as.vector(distances)
  closest <- order(s_distances)[1:pt$n_chirps]
  chirps_closest = chirps_[chirps_$geometry %in% chirps_$geometry[closest],1:505] %>% bind_cols(select(as.data.frame(pt), n_chirps))
  colnames(chirps_closest) = gsub(".x","",colnames(chirps_closest))
  
  #ggplot() + 
  #  geom_sf(data = chirps_, aes(fill = `chirps-v2.0.1981.01.x`, col = `chirps-v2.0.1981.01.x`)) +
  #  geom_sf(data = pt, fill = NA, col = 'red', lwd = 0.2) +
  #  geom_sf(data = chirps_closest, fill = NA, col = 'green', lwd = 0.2 )
  
  return(chirps_closest)
}

sttns.points$r = 0.05
table(sttns.points$r)

pol_list = list()
for(i in 1:nrow(sttns.points)){
  pol_list[[i]]  = pol_(sttns.points[i,], sttns.points$r[i])  
}

df_ = data.frame(nchirps_ = sttns.points$n_chirps, nrow_ = NA)
df_$ID = sttns.points$ID
  
for(i in 1:nrow(sttns.points)){
  df_$nrow_[i] =  nrow(pol_list[[i]])
}

df_ %>% count(nchirps_, nrow_)
sum(df_$nrow_) # 2062 pixeles para validar 846 estaciones (~2.5 pix, x estación)

# Check point 1 ----
save.image('df_pol_chirps.RData')
load('df_pol_chirps.RData')

# Corrección estaciones con CHIRPS recortados cn shp depto ----
# 15 estaciones donde son diferentes el número de pixeles a seleccionar que 
ID_pix = df_[df_$nchirps != df_$nrow_, 'ID']

chirps_df = arrow::read_parquet("chirps_df_v1.parquet")

chirps.point <- st_as_sf(x = chirps_df, 
                         coords = c("x", "y"),
                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

sttns.points[1,'r'] = 0.1
sttns.points[2,'r'] = 0.85
ID_pix = which(sttns.points$ID %in% ID_pix)
#ID_pix = ID_pix[3:length(ID_pix)]

for(i in ID_pix){
  print(i)
  pol_list[[i]]  = pol_(sttns.points[i,], sttns.points$r[i])  
}

for(i in ID_pix){
  df_$nrow_[i] =  nrow(pol_list[[i]])
}

for(i in 1:length(pol_list)){
  pol_list[[i]]$ID =  df_$ID[i]
}

# 12 filas donde son diferentes el número de pixeles a seleccionar que 
ID_pix = df_[df_$nchirps != df_$nrow_, 'ID']

shp_depto <- st_read("Departamentos202208_shp/Depto.shp")
mundocol <- st_read("admin00/admin00_vecinos.shp")

ggplot() + 
  geom_sf(data = shp_depto) +
  geom_sf(data = sttns.points[sttns.points$ID%in% ID_pix,], col = 'red', lwd = 0.2) +
  geom_sf_text(aes(label = ID), data = sttns.points[sttns.points$ID %in% ID_pix,], size = 3) 

sum(df_$nrow_) # 2093 pixeles para validar 846 estaciones (~2.5 pix, x estación)

# join con sttns.points por ID ----
pol2_list = list()
for(i in 1:length(pol_list)){
  pol2_list[[i]] = left_join(as.data.frame(pol_list[[i]]), 
                             as.data.frame(select(sttns.points,-c('n_chirps', 'r'))), by = 'ID')
  names(pol2_list[[i]]) = gsub("-", ".",names(pol2_list[[i]]))
}

# Pivotar de columnas a filas ----

for(i in 1:length(pol2_list)){
  pol2_list[[i]] = pol2_list[[i]] %>%
    select(-matches("^(19|20|geometry.y)")) %>%
    pivot_longer(
      cols = starts_with("chirps.v2.0."),
      names_to = "Date",
      values_to = "chirps.v2.0"
    ) %>% mutate(Date = gsub('chirps.v2.0.','', Date)) %>% 
    left_join(
      pol2_list[[i]] %>% 
        select(-matches("^(chirps.v2.0.|geometry.x)")) %>% 
        pivot_longer(
          cols = matches("^(19|20)"),
          names_to = "Date",
          values_to = "sttns"
        )%>% distinct(), by = 'Date'
    ) %>% select(-c('ID.y', 'n_chirps.y')) %>% 
    rename('ID'='ID.x', 'n_chirps'='n_chirps.x')
}

# Calcula las medidas de desempeño ----
# r pearson ----
R_pearson = c()
for(i in 1:length(pol2_list)){
  R_pearson[i] = cor(pol2_list[[i]]$chirps.v2.0, pol2_list[[i]]$sttns, use = "complete.obs")
}

quantile(R_pearson)
mean(R_pearson) # 0.619
# indica que en promedio es una correlación aceptable 0.6 < r < 0.75
hist(R_pearson)

# MAE (pred vs obs) ----
f_MAE = function(df){
  MAE = df %>%
    mutate(abs_error = abs(chirps.v2.0 - sttns)) %>%
    summarise(MAE = mean(abs_error, na.rm = TRUE)) %>%
    pull(MAE)
  return(MAE)
}

MAE_ = c()
for(i in 1:length(pol2_list)){
  MAE_[i] = f_MAE(pol2_list[[i]])
}

quantile(MAE_)
mean(MAE_) # 69.07
# indica que, en promedio, la diferencia absoluta entre Chirps y las observaciones
## es de aproximadamente 69mm.
hist(MAE_)

# PBIAS ----
pol2_list

# R^2 (pretende explicar la variabilidad del modelo) ----
f_R2 = function(df){
  r2 <- df %>%
    summarise(
      ss_total = sum((sttns - mean(sttns, na.rm = TRUE))^2, na.rm = TRUE),
      ss_res = sum((sttns - chirps.v2.0)^2, na.rm = TRUE)
    ) %>%
    mutate(R2 = 1 - (ss_res / ss_total)) %>%
    pull(R2)
  return(r2)
}

R2_ = c()
for(i in 1:length(pol2_list)){
  R2_[i] = f_R2(pol2_list[[i]])
}

quantile(R2_)
mean(R2_) # 0.33
# indica que, en promedio, la diferencia absoluta entre Chirps y las observaciones
## es de aproximadamente 69mm.
hist(R2_)
# NSE (Nash-Sutcliffe Efficiency) ----

f_NSE = function(df){
  NSE = df %>%
    summarise(
      mean_obs = mean(sttns, na.rm = TRUE),
      num = sum((sttns - chirps.v2.0)^2, na.rm = TRUE),
      den = sum((sttns - mean_obs)^2, na.rm = TRUE)
    ) %>%
    mutate(NSE_ = 1 - (num/den)) %>%
    pull(NSE_)
  return(NSE)
}

NSE_ = c()
for(i in 1:length(pol2_list)){
  NSE_[i] = f_NSE(pol2_list[[i]])
}

quantile(NSE_)
mean(NSE_) # 0.33
# indica que, Chirps tiene una buena predicción, mejor que usar 
## la media de los valores observados como predicción, sin embargo
### el promedio está lejos de 1, aún hay un margen significativo para mejorar
hist(NSE_)

# Check point 2 ----
save.image('df_eval_metrics.RData')
load('df_eval_metrics.RData')

########################################################
#set.seed(123)
#n = sample(1:nrow(sttns.points), 1)
x = 95 #745
#which(sttns.points$ID %in% 15075501)
r = 0.1
pt = sttns.points[x,]
print(sttns.points$n_chirps[x])

pol = data.frame(matrix(c(pt$geometry[[1]] +c(-r,-r),
                          pt$geometry[[1]] +c(-r,r),
                          pt$geometry[[1]] +c(r,r),
                          pt$geometry[[1]] +c(r,-r)), nc = 2, byrow = F
                        )
                 )
names(pol) = c("X", "Y")
pol = pol %>%
  st_as_sf(coords = c("X", "Y"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")                            
join_chirps_pol = st_join(pol, chirps.point)
chirps_ = st_join(chirps.point, join_chirps_pol)
chirps_ = na.omit(chirps_) %>% distinct(geometry, .keep_all = T)
distances <- st_distance(pt, chirps_)
s_distances <- as.vector(distances)
closest <- order(s_distances)[1:pt$n_chirps]
chirps_closest = chirps_[chirps_$geometry %in% chirps_$geometry[closest],1:505] %>% bind_cols(select(as.data.frame(pt), n_chirps))
colnames(chirps_closest) = gsub(".x","",colnames(chirps_closest))

chirps_1981_01 = filter(chirps.point[1],`chirps-v2.0.1981.01`>0)

join_chirps_pol = st_join(pol, chirps_1981_01)
chirps_1981_01 = st_join(chirps_1981_01, join_chirps_pol)
chirps_1981_01 = filter(chirps_1981_01, !is.na(`chirps-v2.0.1981.01.y`)) %>% 
  select(-`chirps-v2.0.1981.01.y`)

ggplot() + 
  geom_sf(data = shp_depto %>% filter(DeNombre %in% "La Guajira")) +
  #geom_sf(data = chirps.point, aes(fill = `chirps-v2.0.1981.01`, col = `chirps-v2.0.1981.01`)) +
  geom_sf(data = pol, fill = NA, col = 'blue', lwd = 0.2) +
  geom_sf(data = pt, fill = NA, col = 'red', lwd = 0.2) +
  geom_sf(data = chirps_closest, fill = NA, col = 'green', lwd = 0.2 )
  


# Descarga `Elevación` del SRTM 90m ----

SRTM <- numeric(length(pcp_col$X))
for(i in 1:length(pcp_col$X)){
  SRTM[[i]] <- extract(geodata::elevation_3s(lon = pcp_col$X[i], lat = pcp_col$Y[i], path = './tmpr'), 
                       cbind(pcp_col$X[i], pcp_col$Y[i]))[[1]]
}

SRTM[2] = extract(geodata::elevation_3s(lon = pcp_col$X[2], lat = pcp_col$Y[2], path = './tmpr'), 
                  cbind(pcp_col$X[2]-0.006, pcp_col$Y[2]))[[1]]

SRTM[38] = extract(geodata::elevation_3s(lon = pcp_col$X[38], lat = pcp_col$Y[38], path = './tmpr'), 
                   cbind(pcp_col$X[38], pcp_col$Y[38]-0.006))[[1]]

SRTM[91] = extract(geodata::elevation_3s(lon = pcp_col$X[91], lat = pcp_col$Y[91], path = './tmpr'), 
                   cbind(pcp_col$X[91], pcp_col$Y[91]-0.03))[[1]]

#Q1 = as.vector(quantile(SRTM)[2])
#Q2 = as.vector(quantile(SRTM)[3])
#Q3 = as.vector(quantile(SRTM)[4])

SRTM_30.sttns$geodata = SRTM
SRTM_30.sttns %>% filter(ID %in% sample(SRTM_30.sttns$ID, size = 6, replace = FALSE))

#pcp_col$SRTM = SRTM
#pcp_col$n_chirps = ifelse(pcp_col$SRTM <= Q1,4,
#                          ifelse(pcp_col$SRTM > Q1 & pcp_col$SRTM <= Q2, 3,
#                                 ifelse(pcp_col$SRTM > Q2 & pcp_col$SRTM <= Q3,2,
#                                        ifelse(pcp_col$SRTM > Q3, 1, NA))))
#table(pcp_col$n_chirps)

