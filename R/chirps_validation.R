library(arrow)
library(dplyr)
library(sf)
library(terra)
library(fpp3)
library(openxlsx)
library(geodata)
library(ggplot2)
library(ggtext) 
library(tidyverse)
library(tsibble)
library(tseries)
library(viridis)
library(nortest)
library(car)
library(lmtest)

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

sttns.points$n_chirps = 1 # (oct/2024)
sttns.points$r = 0.05
table(sttns.points$r)

pol_list = list()
for(i in 1:nrow(sttns.points)){
  pol_list[[i]]  = pol_(sttns.points[i,], sttns.points$r[i])  
}

df_ = data.frame(nchirps_ =  sttns.points$n_chirps, nrow_ = NA)
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

# Evalua la estacionalidad de las series ----
# La estacionariedad de una serie temporal implica que sus propiedades estadísticas
## no cambian con el tiempo, en otras palabras no hay tendencia.

idx_nchrips1 = which(df_$nchirps_==1)
pol3_list = pol2_list[idx_nchrips1]
# ADF test: serie chirps original ----
# crear el ciclo for para calcular los ordenes de cada serie CHIRPS
# Tsay recomienda diferenciar ordinalmente la serie para ajustar el modelo AR
m_order = c()
k = 1
for(i in 1:length(pol3_list)){
  m_order[k] <- tryCatch(
    {
      # Intenta ajustar el modelo AR
      ar_obj <- ar(diff(pol3_list[[i]]$chirps.v2.0), method = 'mle')
      ar_obj$order  # Devuelve el orden del modelo ajustado
    },
    error = function(e) {
      # Si hay un error, devuelve -1
      -1
    }
  )
  # Incrementa el índice k
  k <- k + 1
}

table(m_order)

# calcula el orden sin diferenciar ordinalmente
m_order2 = c()
k = 1
for(i in 1:length(pol3_list)){
  m_order2[k] <- tryCatch(
    {
      # Intenta ajustar el modelo AR
      ar_obj <- ar(pol3_list[[i]]$chirps.v2.0, method = 'mle')
      ar_obj$order  # Devuelve el orden del modelo ajustado
    },
    error = function(e) {
      # Si hay un error, devuelve -1
      -1
    }
  )
  # Incrementa el índice k
  k <- k + 1
}

table(m_order2)

## gráficos de las estaciones que arrojan error al ajustar un modelo AR ----
no_order = which(m_order==-1)
ID_no_order = c()
k = 1
for(i in no_order){
  ID_no_order[k] = unique(pol3_list[[i]]$ID)
  k <- k+1
}

plot(pol3_list[[7]]$chirps.v2.0, type = 'l', main = unique(pol3_list[[7]]$ID))
ggplot() + 
  geom_sf(data = shp_depto) +
  geom_sf(data = sttns.points[sttns.points$ID%in% ID_no_order,], col = 'red', lwd = 0.2) +
  geom_sf_text(aes(label = ID), data = sttns.points[sttns.points$ID %in% ID_no_order,], size = 3) 

## prueba ADF librería aTSA ----
m_order = ifelse(m_order==-1, 12, m_order)

len_ = length(pol3_list)
df_ADF = data.frame(lag_t1 = rep(0, len_), p.value_t1 = rep(0, len_),
                    lag_t2 = rep(0, len_), p.value_t2 = rep(0, len_),
                    lag_t3 = rep(0, len_), p.value_t3 = rep(0, len_))

k = 1
for(i in 1:length(pol3_list)){
  # cambiar m_order por m_order2
  l_ADF = aTSA::adf.test(pol3_list[[i]]$chirps.v2.0, nlag = m_order2[i], output = FALSE)
  df_ADF$lag_t1[k] = which(l_ADF$type1[,'p.value']==max(l_ADF$type1[,'p.value']))
  df_ADF$p.value_t1[k] = max(l_ADF$type1[,'p.value'])
  df_ADF$lag_t2[k] = which(l_ADF$type2[,'p.value']==max(l_ADF$type2[,'p.value']))
  df_ADF$p.value_t2[k] = max(l_ADF$type2[,'p.value'])
  df_ADF$lag_t3[k] = which(l_ADF$type3[,'p.value']==max(l_ADF$type3[,'p.value']))
  df_ADF$p.value_t3[k] = max(l_ADF$type3[,'p.value'])
  k = k+1
}

df_ADF$RD_t1 = ifelse(df_ADF$p.value_t1<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF$RD_t2 = ifelse(df_ADF$p.value_t2<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF$RD_t3 = ifelse(df_ADF$p.value_t3<0.05, 'R_H0 - estacionaria', 'NO R_H0')

table(df_ADF$RD_t1)
#df_ADF[(df_ADF$RD_t1==df_ADF$RD_t2) & (df_ADF$RD_t2==df_ADF$RD_t3),]
#plot(pol3_list[[91]]$chirps.v2.0, type  ='l')
#table(df_ADF$RD_t1,df_ADF$RD_t2, df_ADF$RD_t3)

# prueba ADF librería urca
interp_urdf_FM <- function(urdf, level="5pct") {
  if(class(urdf) != "ur.df") stop('parameter is not of class ur.df from urca package')
  if(!(level %in% c("1pct", "5pct", "10pct") ) ) stop('parameter level is not one of 1pct, 5pct, or 10pct')
  
  #cat( paste("At the", level, "level:\n") )
  if(urdf@model == "none") {
    #cat("The model is of type none\n")
    tau1_crit = urdf@cval["tau1",level]
    tau1_teststat = urdf@teststat["statistic","tau1"]
    tau1_teststat_wi_crit = tau1_teststat > tau1_crit
    if(tau1_teststat_wi_crit) {
      rslt = "NO R_H0"
      #cat("tau1: The null hypothesis is not rejected, unit root is present\n")
    } else {
      #cat("tau1: The null hypothesis is rejected, unit root is not present\n")
      rslt ="R_H0 - estacionaria"
    }}
  return(rslt)}

df_ur = data.frame(RD_t1 = rep('', 846))
for(i in 1:length(pol3_list)){
  ADF_ur = urca::ur.df(pol3_list[[i]]$chirps.v2.0, type  = "none", lags = m_order2[i], selectlags = "AIC")
  df_ur$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur$RD_t1)

# diff estacional ----
# Periodograma
len_ = length(pol3_list)
df_periodograma = data.frame(freq = rep(0, len_),
                             period = rep(0, len_))

# encontrar el periodo del ciclo estacional
for(i in 1:length(pol3_list)){
  Periodgrama_ = spectrum(pol3_list[[i]]$chirps.v2.0,log='no')
  loc_ = which.max(Periodgrama_$spec)
  df_periodograma$freq[i] = Periodgrama_$freq[loc_]
  df_periodograma$period[i] = 1/Periodgrama_$freq[loc_]
}

#spectrum(pol3_list[[91]]$chirps.v2.0,log='no')
#abline(v = 1/12, col = 'red')
## plot muestra pixeles con periodos 4 y 6 ----
p4 = which(round(df_periodograma$period)==4)
p6 = which(round(df_periodograma$period)==6)
set.seed(22112024)
s_p6  = sample(p6, size = 1, replace = FALSE)

spectrum(pol3_list[[s_p6]]$chirps.v2.0,log='no')
abline(v = 1/12, col = 'red')
abline(v = 1/6, col = 'red')
abline(v = 1/4, col = 'red')

y <- tsibble::tsibble(pol3_list[[s_p6]] %>% mutate(Date = as.Date(str_c(Date, '.01'), format = '%Y.%m.%d')) 
                      %>%  select(c("Date" ,  "chirps.v2.0")), index = Date)

y = y %>% mutate(Date = yearmonth(Date)) %>% 
  as_tsibble(index = Date)

gg_season(y, `chirps.v2.0`)
gg_subseries(y, `chirps.v2.0`)

table(round(df_periodograma$period))

## diferencias estacionales ----
d_chirps.v2.0 = list()
for(i in 1:length(pol3_list)){
  d_chirps.v2.0[[i]] = diff(pol3_list[[i]]$chirps.v2.0, 
                            lag = round(df_periodograma$period[i]))
}

for(i in 1:length(pol3_list)){
  pol3_list[[i]]$d_chirps.v2.0 = c(rep(NA, round(df_periodograma$period[i])),
                                   diff(pol3_list[[i]]$chirps.v2.0, 
                                        lag = round(df_periodograma$period[i])))
}

# crear el ciclo for para calcular los ordenes de cada serie CHIRPS diferenciadas estacionalmente
# !!! dos escenarios sin diferenciar 
m_order4 = c()
k = 1
for(i in 1:length(pol3_list)){
  m_order4[k] <- tryCatch(
    {
      # Intenta ajustar el modelo AR
      #ar_obj <- ar(diff(d_chirps.v2.0[[i]]), method = 'mle')
      ar_obj <- ar(d_chirps.v2.0[[i]], method = 'mle')
      ar_obj$order  # Devuelve el orden del modelo ajustado
    },
    error = function(e) {
      # Si hay un error, devuelve -1
      -1
    }
  )
  # Incrementa el índice k
  k <- k + 1
}

table(m_order4)
#m_order4 = ifelse(m_order3==-1,12, m_order4)
## prueba ADF librería aTSA ----
df_ADF2 = data.frame(lag_t1 = rep(0, len_), p.value_t1 = rep(0, len_),
                     lag_t2 = rep(0, len_), p.value_t2 = rep(0, len_),
                     lag_t3 = rep(0, len_), p.value_t3 = rep(0, len_))

k = 1
for(i in 1:length(pol3_list)){
  l_ADF = aTSA::adf.test(d_chirps.v2.0[[i]], nlag = m_order4[i], output = FALSE)
  df_ADF2$lag_t1[k] = which(l_ADF$type1[,'p.value']==max(l_ADF$type1[,'p.value']))
  df_ADF2$p.value_t1[k] = max(l_ADF$type1[,'p.value'])
  df_ADF2$lag_t2[k] = which(l_ADF$type2[,'p.value']==max(l_ADF$type2[,'p.value']))
  df_ADF2$p.value_t2[k] = max(l_ADF$type2[,'p.value'])
  df_ADF2$lag_t3[k] = which(l_ADF$type3[,'p.value']==max(l_ADF$type3[,'p.value']))
  df_ADF2$p.value_t3[k] = max(l_ADF$type3[,'p.value'])
  k = k+1
}

df_ADF2$RD_t1 = ifelse(df_ADF2$p.value_t1<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF2$RD_t2 = ifelse(df_ADF2$p.value_t2<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF2$RD_t3 = ifelse(df_ADF2$p.value_t3<0.05, 'R_H0 - estacionaria', 'NO R_H0')

table(df_ADF2$RD_t1)

## prueba ADF librería urca ----
df_ur = data.frame(RD_t1 = rep('', 846))
for(i in 1:length(pol3_list)){
  ADF_ur = urca::ur.df(d_chirps.v2.0[[i]], type  = "none", lags = m_order4[i], selectlags = "AIC")
  df_ur$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur$RD_t1)

# Removiendo componente estacional usando LOESS ----
set.seed(20112024)
s = sample(1:length(pol3_list), size = 1)

y <- tsibble::tsibble(pol3_list[[s]] %>% mutate(Date = as.Date(str_c(Date, '.01'), format = '%Y.%m.%d')) 
                      %>%  select(c("Date" ,  "chirps.v2.0")), index = Date)

y = y %>% mutate(Date = yearmonth(Date)) %>% 
  as_tsibble(index = Date)

# descomposición STL usando librería base ----
# t.window = 21 (fpp3) auto
#nextodd <- function(x) {
#  x <- round(x)
#  if (x%%2 == 0) 
#    x <- x + 1
#  as.integer(x)}

#nextodd(ceiling(1.5 * frequency(y)/(1 - 1.5/(10 * as.integer(length(y)) + 1))))

stl_chirps.v2.0 = list()
for(i in 1:length(pol3_list)){
  y <- tsibble::tsibble(pol3_list[[i]] %>% mutate(Date = as.Date(str_c(Date, '.01'), format = '%Y.%m.%d')) 
                        %>%  select(c("Date" ,  "chirps.v2.0")), index = Date)
  
  y = y %>% mutate(Date = yearmonth(Date)) %>% 
    as_tsibble(index = Date)
  # t.window = 21 (fpp3) auto
  #nextodd(ceiling(1.5 * frequency(y)/(1 - 1.5/(10 * as.integer(length(y)) + 1))))
  fit1 <- stl(y, s.window = "periodic")
  # t.window = 13 (fpp2)
  #fit11<- stl(y, t.window = 13, s.window = "periodic")
  # Selecting a shorter trend window, avoid produces a trend-cycle component that is too rigid
  #fit12<- stl(y, t.window = 7, s.window = "periodic")
  stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit1$time.series[, "seasonal"]
  #stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit11$time.series[, "seasonal"]
  #stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit12$time.series[, "seasonal"]
}

m_order5 = c()
k = 1
for(i in 1:length(pol3_list)){
  m_order5[k] <- tryCatch(
    {
      # Intenta ajustar el modelo AR
      ar_obj <- ar(diff(stl_chirps.v2.0[[i]]), method = 'mle')
      #ar_obj <- ar(stl_chirps.v2.0[[i]], method = 'mle')
      ar_obj$order  # Devuelve el orden del modelo ajustado
    },
    error = function(e) {
      # Si hay un error, devuelve -1
      -1
    }
  )
  # Incrementa el índice k
  k <- k + 1
}

table(m_order5)
m_order5 = ifelse(m_order5==-1,12, m_order5)

## prueba ADF librería aTSA ----
len_ = length(pol3_list)
df_ADF3 = data.frame(lag_t1 = rep(0, len_), p.value_t1 = rep(0, len_),
                     lag_t2 = rep(0, len_), p.value_t2 = rep(0, len_),
                     lag_t3 = rep(0, len_), p.value_t3 = rep(0, len_))

k = 1
for(i in 1:length(pol3_list)){
  l_ADF = aTSA::adf.test(stl_chirps.v2.0[[i]], nlag = m_order5[i], output = FALSE)
  df_ADF3$lag_t1[k] = which(l_ADF$type1[,'p.value']==max(l_ADF$type1[,'p.value']))
  df_ADF3$p.value_t1[k] = max(l_ADF$type1[,'p.value'])
  df_ADF3$lag_t2[k] = which(l_ADF$type2[,'p.value']==max(l_ADF$type2[,'p.value']))
  df_ADF3$p.value_t2[k] = max(l_ADF$type2[,'p.value'])
  df_ADF3$lag_t3[k] = which(l_ADF$type3[,'p.value']==max(l_ADF$type3[,'p.value']))
  df_ADF3$p.value_t3[k] = max(l_ADF$type3[,'p.value'])
  k = k+1
}

df_ADF3$RD_t1 = ifelse(df_ADF3$p.value_t1<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF3$RD_t2 = ifelse(df_ADF3$p.value_t2<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF3$RD_t3 = ifelse(df_ADF3$p.value_t3<0.05, 'R_H0 - estacionaria', 'NO R_H0')

table(df_ADF3$RD_t1)

## prueba ADF librería urca ----
df_ur = data.frame(RD_t1 = rep('', 846))
for(i in 1:length(pol3_list)){
  ADF_ur = urca::ur.df(stl_chirps.v2.0[[i]], type  = "none", lags = m_order5[i], selectlags = "AIC")
  df_ur$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur$RD_t1)

# descomposición STL usando librería feast ----
f_stl_chirps.v2.0 = list()
for(i in 1:length(pol3_list)){
  y <- tsibble::tsibble(pol3_list[[i]] %>% mutate(Date = as.Date(str_c(Date, '.01'), format = '%Y.%m.%d')) 
                        %>%  select(c("Date" ,  "chirps.v2.0")), index = Date)
  
  y = y %>% mutate(Date = yearmonth(Date)) %>% 
    as_tsibble(index = Date)
  
  # The default setting for monthly data is trend(window=21)
  # t.window = 13 (fpp2)
  # Selecting a shorter trend window, avoid produces a trend-cycle component that is too rigid
  fit2 = y |>
    model(
      #STL(`chirps.v2.0` ~ trend(window = 13) +
      STL(`chirps.v2.0` ~ 
            season(window = "periodic"),
          robust = TRUE)) |>
    components() 
  f_stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit2$season_year
}

set.seed(24112024)
s = sample(1:length(pol3_list), size = 1)

min_length <- min(length(stl_chirps.v2.0[[s]]), length(d_chirps.v2.0[[s]]))
s1_adj <- window(stl_chirps.v2.0[[s]], end = time(stl_chirps.v2.0[[s]])[min_length])
y_adj = window(pol3_list[[s]]$chirps.v2.0, end = time(pol3_list[[s]]$chirps.v2.0)[min_length])
s2_adj <- d_chirps.v2.0[[s]]
s3_adj <- window(f_stl_chirps.v2.0[[s]], end = time(f_stl_chirps.v2.0[[s]])[min_length])

plot(s2_adj, type = "l", col = "red", lwd = 2, lty = 2, xlab = "Index", ylim  = c(-200,700))
lines(1:length(s2_adj), s1_adj, type = "l", col = "blue", lwd = 2, lty =2)
lines(1:length(s2_adj), y_adj, type = "l", col = "black", lwd = 2)
lines(1:length(s2_adj), s3_adj, type = "l", col = "gray", lwd = 2, lty =2)

legend("topleft", legend = c("serie original CHIRPS", "sin comp. estacional - stats::stl", 
                             "sin comp. estacional - feast::STL", "dif estacional"), 
       col = c("black", "blue", "gray", "red"), lwd = 2, lty = c(1,2,2,2), bty = "n")


m_order6 = c()
k = 1
for(i in 1:length(pol3_list)){
  m_order6[k] <- tryCatch(
    {
      # Intenta ajustar el modelo AR
      ar_obj <- ar(diff(f_stl_chirps.v2.0[[i]]), method = 'mle')
      #ar_obj <- ar(f_stl_chirps.v2.0[[i]], method = 'mle')
      ar_obj$order  # Devuelve el orden del modelo ajustado
    },
    error = function(e) {
      # Si hay un error, devuelve -1
      -1
    }
  )
  # Incrementa el índice k
  k <- k + 1
}

table(m_order6)
m_order6 = ifelse(m_order6==-1,12, m_order6)

## prueba ADF librería aTSA ----
len_ = length(pol3_list)
df_ADF4 = data.frame(lag_t1 = rep(0, len_), p.value_t1 = rep(0, len_),
                     lag_t2 = rep(0, len_), p.value_t2 = rep(0, len_),
                     lag_t3 = rep(0, len_), p.value_t3 = rep(0, len_))

k = 1
for(i in 1:length(pol3_list)){
  l_ADF = aTSA::adf.test(f_stl_chirps.v2.0[[i]], nlag = m_order6[i], output = FALSE)
  df_ADF4$lag_t1[k] = which(l_ADF$type1[,'p.value']==max(l_ADF$type1[,'p.value']))
  df_ADF4$p.value_t1[k] = max(l_ADF$type1[,'p.value'])
  df_ADF4$lag_t2[k] = which(l_ADF$type2[,'p.value']==max(l_ADF$type2[,'p.value']))
  df_ADF4$p.value_t2[k] = max(l_ADF$type2[,'p.value'])
  df_ADF4$lag_t3[k] = which(l_ADF$type3[,'p.value']==max(l_ADF$type3[,'p.value']))
  df_ADF4$p.value_t3[k] = max(l_ADF$type3[,'p.value'])
  k = k+1
}

df_ADF4$RD_t1 = ifelse(df_ADF4$p.value_t1<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF4$RD_t2 = ifelse(df_ADF4$p.value_t2<0.05, 'R_H0 - estacionaria', 'NO R_H0')
df_ADF4$RD_t3 = ifelse(df_ADF4$p.value_t3<0.05, 'R_H0 - estacionaria', 'NO R_H0')

table(df_ADF4$RD_t1)

## prueba ADF librería urca ----
df_ur = data.frame(RD_t1 = rep('', 846))
for(i in 1:length(pol3_list)){
  ADF_ur = urca::ur.df(stl_chirps.v2.0[[i]], type  = "none", lags = m_order6[i], selectlags = "AIC")
  df_ur$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur$RD_t1)

#fit2 %>% autoplot()

deseasoned2 <- y$chirps.v2.0 - fit2$season_year
deseasoned21 <- y$chirps.v2.0 - fit21$season_year
deseasoned22 <- y$chirps.v2.0 - fit22$season_year

cbind(deseasoned, deseasoned2)

plot.ts(deseasoned,
        cex.main = 0.85)
lines(deseasoned, col = 'red')

#set.seed(11102024)
#ts_sample = sample(1:length(pol2_list), 10, replace = FALSE)
#pol2_sample = pol2_list[ts_sample]
#saveRDS(pol2_sample, 'Muestra_CHIRPS.RDS')


# Calculas medidas de desempeño ----
# r pearson ----
for(i in 1:length(pol3_list)){
  pol3_list[[i]]$mes = format(as.Date(paste0(pol3_list[[i]]$Date, ".01"), format = "%Y.%m.%d"), "%B")
}

# crear variable categorica (DJF/MAM/JJA/SON) a partir de la variable mes 
DJF_ = c("diciembre", "enero", "febrero" )
MAM_ = c("marzo", "abril", "mayo")
JJA_ = c("junio", "julio", "agosto")
SON_ = c("septiembre", "octubre", "noviembre")
for(i in 1:length(pol3_list)){
  pol3_list[[i]]$cat_mes = ifelse(pol3_list[[i]]$mes %in% DJF_, 
                                  'DJF', ifelse(pol3_list[[i]]$mes %in% MAM_,
                                                'MAM', ifelse(pol3_list[[i]]$mes %in% JJA_,
                                                              'JJA', ifelse(pol3_list[[i]]$mes %in% SON_, 'SON', NA))))
}

R_pearson = c()
for(i in 1:length(pol2_list)){
  R_pearson[i] = cor(pol2_list[[i]]$chirps.v2.0, pol2_list[[i]]$sttns, use = "complete.obs")
}

R_Spearman = c()
for(i in 1:length(pol3_list)){
  R_Spearman[i] = cor(pol3_list[[i]]$chirps.v2.0, pol3_list[[i]]$sttns, 
                      use = "complete.obs", method = "spearman")
}

d_sttns = list()
for(i in 1:length(pol3_list)){
  d_sttns[[i]] = diff(pol3_list[[i]]$sttns, 
                      lag = round(df_periodograma$period[i]))
}

for(i in 1:length(pol3_list)){
  pol3_list[[i]]$d_sttns = c(rep(NA, round(df_periodograma$period[i])),
                             diff(pol3_list[[i]]$sttns, 
                                  lag = round(df_periodograma$period[i])))
}

# Boxplot por mes ----
R_Spearman_mes = list()
for(i in 1:length(pol3_list)){
  meses <- unique(pol3_list[[i]]$mes)
  # Calculamos el coeficiente de correlación de Spearman para cada mes
  corrs <- sapply(meses, function(m){
    # Filtramos los datos para el mes 'm'
    data_mes <- subset(pol3_list[[i]], mes == m)
    # Calculamos el coeficiente de correlación de Spearman
    cor_value <- cor(data_mes$d_chirps.v2.0, data_mes$d_sttns, use = "complete.obs", method = "spearman")
    return(cor_value)
  })
  # Almacenamos los resultados en el listado de R_Spearman
  R_Spearman_mes[[i]] <- corrs
}

R_Spearman_df <- as.data.frame(R_Spearman_mes)
R_Spearman_df <- t(R_Spearman_df)
R_Spearman_df <- as.data.frame(R_Spearman_df)
colnames(R_Spearman_df) <- unique(pol3_list[[1]]$mes)
rownames(R_Spearman_df) <- 1:nrow(R_Spearman_df)

R_Spearman_long <- R_Spearman_df %>%
  pivot_longer(cols = everything(), 
               names_to = "Mes", 
               values_to = "Spearman_R")

R_Spearman_long$Mes <- factor(R_Spearman_long$Mes, 
                              levels = c("enero", "febrero", "marzo", "abril", "mayo", 
                                         "junio", "julio", "agosto", "septiembre", "octubre", 
                                         "noviembre", "diciembre"))

ggplot(R_Spearman_long, aes(x = Mes, y = Spearman_R)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = expression("Boxplot de r"["s"] ~ "por mes"),
       x = "Mes", 
       y = expression(r["s"])) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot por cat_mes (DJF, MAM, JJA, SON) ----
R_Spearman_mes = list()
for(i in 1:length(pol3_list)){
  meses <- unique(pol3_list[[i]]$cat_mes)
  # Calculamos el coeficiente de correlación de Spearman para cada mes
  corrs <- sapply(meses, function(m){
    # Filtramos los datos para el mes 'm'
    data_mes <- subset(pol3_list[[i]], cat_mes == m)
    # Calculamos el coeficiente de correlación de Spearman
    cor_value <- cor(data_mes$d_chirps.v2.0, data_mes$d_sttns, use = "complete.obs", method = "spearman")
    return(cor_value)
  })
  # Almacenamos los resultados en el listado de R_Spearman
  R_Spearman_mes[[i]] <- corrs
}

R_Spearman_df <- as.data.frame(R_Spearman_mes)
R_Spearman_df <- t(R_Spearman_df)
R_Spearman_df <- as.data.frame(R_Spearman_df)
colnames(R_Spearman_df) <- unique(pol3_list[[1]]$cat_mes)
rownames(R_Spearman_df) <- 1:nrow(R_Spearman_df)

R_Spearman_long <- R_Spearman_df %>%
  pivot_longer(cols = everything(), 
               names_to = "Cat_mes", 
               values_to = "Spearman_R")

R_Spearman_long$Cat_mes <- factor(R_Spearman_long$Cat_mes, 
                                  levels = c("DJF", "MAM", "JJA", "SON"))

# ANOVA por grupo de meses ----
modelo <- aov(Spearman_R ~ Cat_mes, data = R_Spearman_long)
modelo_ = summary(modelo)
p_value = modelo_[[1]]$`Pr(>F)`[1]
print(ifelse(p_value <= 0.05, "R. H0, concluye dif significativas entre al menos dos grupos",
             "NO R. H0, concluye no hay dif significativas"))

tukey <- TukeyHSD(modelo)
#par(mar = c(5, 10, 4, 2)) 
plot(tukey, las =1)

# validación de supuestos
## Normalidad
plot(density(modelo$residuals))
qqnorm(modelo$residuals)
qqline(modelo$residuals, col = "red")

lillie.test(modelo$residuals)
shapiro.test(modelo$residuals)

## Homocedasticidad
leveneTest(Spearman_R ~ Cat_mes, data = R_Spearman_long)
bptest(modelo)

## Autocorrelación
dwtest(modelo)

# Calcular los valores atípicos por mes
outliers_df <- R_Spearman_long %>%
  group_by(Cat_mes) %>%
  summarise(
    Q1 = quantile(Spearman_R, 0.25, na.rm = TRUE),
    Q3 = quantile(Spearman_R, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    outlier_count = sum(Spearman_R < lower_bound | Spearman_R > upper_bound, na.rm = TRUE),
    total_count = n(),
    outlier_pct = round((outlier_count / total_count) * 100, 1) # Porcentaje de atípicos
  )

# Unir el dataframe con los datos originales para graficar
R_Spearman_long <- left_join(R_Spearman_long, outliers_df, by = "Cat_mes")

ggplot(R_Spearman_long, aes(x = Cat_mes, y = Spearman_R)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  #geom_hline(yintercept = 0.68, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(label = paste0(outlier_pct, "% \n outliers"), y = min(Spearman_R, na.rm = TRUE) + 0.05), 
            size = 3, color = "black", fontface = "bold") +
  theme_minimal() +
  labs(title = expression("Boxplot de r"["s"] ~ "por grupo de meses"),
       x = "Meses", 
       y = expression(r["s"])) 

d_R_pearson = c()
for(i in 1:length(pol2_list)){
  d_R_pearson[i] = cor(d_chirps.v2.0[[i]], d_sttns[[i]], use = "complete.obs")
}

d_R_Spearman = c()
for(i in 1:length(pol3_list)){
  d_R_Spearman[i] = cor(d_chirps.v2.0[[i]], d_sttns[[i]], 
                      use = "complete.obs", method = "spearman")
}

# test de significancia coef corr spearman ----

test_d_R_Spearman = c()
for(i in 1:length(pol3_list)){
  obj_ = cor.test(d_chirps.v2.0[[i]], d_sttns[[i]], 
                             method = "spearman", exact = FALSE)
  test_d_R_Spearman[i] = obj_$p.value
}

# todos significativos a un nivel de confianza del 95%
summary(test_d_R_Spearman)


length(test_d_R_Spearman[test_d_R_Spearman>=0.05])
which(test_d_R_Spearman>=0.01)
# dos valores de rho (0.10 y 0.11) con p-value 0.02 y 0.01
d_R_Spearman[which(test_d_R_Spearman>=0.01)] # valor de rho
test_d_R_Spearman[which(test_d_R_Spearman>=0.01)] # valor de los p-value

df_test = data.frame(test = test_d_R_Spearman, dif_est =  d_R_Spearman)
table(df_test$test==df_test$dif_est)

# Coef de corr spearman cruzado ----
ccf_spearman <- function(x, y, max_lag = 10) {
  #x_rank <- rank(x)
  #y_rank <- rank(y)
  x_rank <- x
  y_rank <- y
  
  lags <- seq(-max_lag, max_lag)
  cor_values <- sapply(lags, function(lag) {
    if (lag < 0) {
      cor(x_rank[1:(length(x) + lag)],y_rank[(1 - lag):length(y)], 
          method = "spearman", use = "complete.obs")
    } else if (lag > 0) {
      cor(x_rank[(1 + lag):length(x)], y_rank[1:(length(y) - lag)], 
          method = "spearman", use = "complete.obs")
    } else {
      cor(x_rank, y_rank, method = "spearman", use = "complete.obs")
    }
  })
  df_cor <- as.data.frame(t(cor_values))
  colnames(df_cor) <- paste0("Lag_", lags)
  
  return(df_cor)
}

d_ccf_R_Spearman = list()
for(i in 1:length(pol3_list)){
  d_ccf_R_Spearman[[i]] = ccf_spearman(d_chirps.v2.0[[i]], d_sttns[[i]], max_lag = 12)
}

df_d_ccf_R_Spearman = d_ccf_R_Spearman %>% bind_rows()
#df_d_ccf_R_Spearman = df_d_ccf_R_Spearman %>% bind_cols(R_s = d_R_Spearman)
#table(df_d_ccf_R_Spearman$Lag_0 == df_d_ccf_R_Spearman$R_s)

# Pivotear el data frame a formato largo
df_d_ccf_R_Spearman <- df_d_ccf_R_Spearman %>%
  pivot_longer(cols = everything(), 
               names_to = "Lag", 
               values_to = "Spearman_Correlation")

df_d_ccf_R_Spearman =df_d_ccf_R_Spearman %>% 
  mutate(Lag = factor(Lag, levels = paste0("Lag_", -12:12)))

# Crear el gráfico de boxplot ccf R Spearman ----
ggplot(df_d_ccf_R_Spearman, aes(x = Lag, y = Spearman_Correlation)) +
  geom_boxplot(fill = "lightblue", color = "black") + 
  theme_minimal() +
  labs(title = "Distribución de la Correlación de Spearman por Lag",
       x = "Lag",
       y = "Coeficiente de Spearman") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotar etiquetas en eje X

# Scatterplot r_s vs altitud ----
df_pearson = data.frame(r = c(R_pearson, d_R_pearson),
                        serie = rep(c("original", "dif estacional"), each = length(R_pearson)))

df_Spearman = data.frame(r_s = c(R_Spearman, d_R_Spearman),
                         serie = rep(c("original", "dif estacional"), each = length(R_Spearman)))

geometry_ <- lapply(pol3_list, function(df_) unique(df_$geometry.y))
geometry_ <- do.call(c, geometry_)  # Combina las geometrías en un solo objeto sf
geometry_ <- st_sfc(geometry_) 
geometry_ <- geometry_[!st_is_empty(geometry_)]
geometry_df <- data.frame(geometry = geometry_)
class(geometry_df$geometry)

# se requiere traer las columnas región, mes y cat_mes
# cómo se traen las regiones más arriba?
df_Spearman_alt = data.frame(r_s = d_R_Spearman,
                             altitud = SRTM_30.sttns$SRTM_30_Col1,
                             geometry = geometry_df$geometry)

# scatterplot v1
ggplot(df_Spearman_alt, aes(x = r_s, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión r"["s"] ~ "vs altitud "),
       x = expression(r["s"]),
       y = "Altitud") +
  theme_minimal()

dim(pol3_list[[1]])

df_stats_r <- df_pearson %>%
  group_by(serie) %>%
  summarise(
    Q1 = quantile(r, 0.25),
    Q2 = median(r),
    Q3 = quantile(r, 0.75)
  )

df_stats <- df_Spearman %>%
  group_by(serie) %>%
  summarise(
    Q1 = quantile(r_s, 0.25),
    Q2 = median(r_s),
    Q3 = quantile(r_s, 0.75)
  )

# Boxplot Corr Pearson ----
ggplot(df_pearson, aes(x = serie, y = r, fill = serie)) +
  geom_boxplot() +
  geom_text(data = df_stats_r, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = df_stats_r, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = df_stats_r, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(title = "Boxplots del coef de Pearson", y = expression(r)) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14) ) 

# Boxplot Corr Spearman ----
ggplot(df_Spearman, aes(x = serie, y = r_s, fill = serie)) +
  geom_boxplot() +
  geom_text(data = df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(title = "Boxplots del coef de Spearman", y = expression(r[s])) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14) ) 

# Distribución espacial de R_s ----
shp_depto <- st_read("Departamentos202208_shp/Depto.shp")
mundocol <- st_read("admin00/admin00_vecinos.shp")

# necesito un df con las columnas: geometry, dpto, Región natural, r_s
# grafico auxiliar de boxplot de r_s por región natural
# grafico dispersión R_s vs Altitud (tendencia a mayor altitud)

df_Spearman_g <- st_as_sf(df_Spearman_alt$geometry, wkt = "geometry")
df_Spearman_g <- st_set_crs(df_Spearman_g, 4326)
df_Spearman_g <- st_transform(df_Spearman_g, st_crs(shp_depto))

df_Spearman_g <- st_join(df_Spearman_g, shp_depto["DeNombre"])
df_Spearman_g$ID = 1:nrow(df_Spearman_g)
df_Spearman_g$DeNombre[df_Spearman_g$ID==2] = "San Andrés Providencia y Santa Catalina"
df_Spearman_g$DeNombre[df_Spearman_g$ID==38] = "Atlántico"
df_Spearman_g$DeNombre[df_Spearman_g$ID==91] = "La Guajira"
df_Spearman_g$DeNombre[df_Spearman_g$ID==191] = "Arauca"
df_Spearman_g$DeNombre[df_Spearman_g$ID==592] = "Vichada"
df_Spearman_g$DeNombre[df_Spearman_g$ID==626] = "Putumayo"
df_Spearman_g$DeNombre[df_Spearman_g$ID==636] = "Putumayo"
df_Spearman_g$DeNombre[df_Spearman_g$ID==274] = "Huila"

caribe_ = c("Atlántico", "Bolívar", "Cesar", "Córdoba" , "La Guajira", "Magdalena", 
         "San Andrés Providencia y Santa Catalina", "Sucre")
pacifica_ = c("Chocó", "Cauca", "Nariño", "Valle del Cauca")
orinoquia_ = c("Arauca", "Casanare", "Meta", "Vichada")
amazonia_ = c("Amazonas", "Caquetá", "Guainía", "Guaviare", "Putumayo", "Vaupés")
andina_ = c("Antioquia", "Boyacá", "Caldas", "Cundinamarca", "Huila", 
            "Norte de Santander", "Quindío", "Risaralda", "Santander", "Tolima")

df_Spearman_g$Region = ifelse(df_Spearman_g$DeNombre %in% caribe_, "Caribe",
                              ifelse(df_Spearman_g$DeNombre %in% pacifica_, "Pacífica",
                                     ifelse(df_Spearman_g$DeNombre %in% orinoquia_, "Orinoquía",
                                            ifelse(df_Spearman_g$DeNombre %in% amazonia_, "Amazonía",
                                                   ifelse(df_Spearman_g$DeNombre %in% andina_, "Andina", "")))))

df_Spearman_g$Region2 = ifelse(df_Spearman_g$DeNombre %in% caribe_, "CAR",
                              ifelse(df_Spearman_g$DeNombre %in% pacifica_, "PAC",
                                     ifelse(df_Spearman_g$DeNombre %in% orinoquia_, "ORN",
                                            ifelse(df_Spearman_g$DeNombre %in% amazonia_, "AMZ",
                                                   ifelse(df_Spearman_g$DeNombre %in% andina_, "AND", "")))))

table(df_Spearman_g$Region)

df_Spearman_g$r_s = df_Spearman_alt$r_s
df_Spearman_g$altitud = df_Spearman_alt$altitud

# panel de scatterplot para cada categoria de región
ggplot(df_Spearman_g, aes(x = r_s, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión r"["s"] ~ "vs altitud por región"),
       x = expression(r["s"]),
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ Region, scales = "free")

# panel de scatterplot para cada categoria de mes
df_Spearman_g2 = df_Spearman_g %>% bind_cols(R_Spearman_df)

R_Spearman_long <- df_Spearman_g2 %>%
  pivot_longer(#cols = enero:diciembre,
               cols = c('DJF', 'MAM', 'JJA', 'SON'),
               #names_to = "Mes", 
               names_to = "Cat_mes", 
               values_to = "Spearman_R") %>%
  #select(ID, altitud, Mes, Spearman_R)
  select(ID, altitud, Cat_mes, Spearman_R) 

#R_Spearman_long$Mes <- factor(R_Spearman_long$Mes, 
#                              levels = c("enero", "febrero", "marzo", "abril", "mayo", 
#                                         "junio", "julio", "agosto", "septiembre", "octubre", 
#                                         "noviembre", "diciembre"))

R_Spearman_long$Cat_mes <- factor(R_Spearman_long$Cat_mes, 
                                  levels = c("DJF", "MAM", "JJA", "SON"))

ggplot(R_Spearman_long, aes(x = Spearman_R, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  #labs(title = expression("Diagrama de dispersión r"["s"] ~ "vs altitud por mes"),
  labs(title = expression("Diagrama de dispersión r"["s"] ~ "vs altitud por grupo de meses"),
       x = expression(r["s"]),
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ Cat_mes, scales = "free")
  #facet_wrap(~ Mes, scales = "free")

# ANOVA no paramétrico ----
kruskal_test <- kruskal.test(r_s ~ Region, data = df_Spearman_g)
print(kruskal_test)

# ANOVA por Region ----
modelo <- aov(r_s ~ Region, data = df_Spearman_g)
modelo_ = summary(modelo)
p_value = modelo_[[1]]$`Pr(>F)`[1]
print(ifelse(p_value <= 0.05, "R. H0, concluye dif significativas entre al menos dos grupos",
             "NO R. H0, concluye no hay dif significativas"))

tukey <- TukeyHSD(modelo)
#par(mar = c(5, 10, 4, 2)) 
plot(tukey, las =1)

# validación de supuestos
## Normalidad
plot(density(modelo$residuals))
qqnorm(modelo$residuals)
qqline(modelo$residuals, col = "red")

lillie.test(modelo$residuals)
shapiro.test(modelo$residuals)

## Homocedasticidad
leveneTest(r_s ~ Region, data = df_Spearman_g)
bptest(modelo)

## Autocorrelación
dwtest(modelo)

# Boxplot por región ----
ggplot(df_Spearman_g) + 
  geom_boxplot(aes(x = Region, y = r_s)) +  # Usamos Región como eje X y r_s como eje Y
  theme_minimal() +  # Estilo de tema minimalista
  labs(x = "Región", y = expression(r["s"]), title = expression("Boxplot de r"["s"] ~ "por Región natural")) #+  # Etiquetas de los ejes y título
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar las etiquetas del eje X si es necesario

ggplot() + 
  geom_sf(data = shp_depto) +  # Muestra el mapa de los departamentos
  geom_sf(data = df_Spearman_g, aes(fill = r_s, color = r_s), size = 1) +
  scale_fill_gradient(low = "red", high = "green", name = "r_s") +# Muestra los puntos con valores de r_s
  scale_color_gradient(low = "red", high = "green") +
  theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de r_s", fill = "Valor de r_s") +  # Título y leyenda
  theme(legend.position = "right") 

quantile(R_Spearman)
mean(R_Spearman) # 0.785275 # 0.619

quantile(R_pearson)
mean(R_pearson) # 0.785275 # 0.619

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
for(i in 1:length(pol3_list)){
  MAE_[i] = f_MAE(pol3_list[[i]])
}

df_list = list()
for(i in 1:length(d_chirps.v2.0)){
  df_list[[i]] = data.frame('chirps.v2.0' = d_chirps.v2.0[[i]],
                            'sttns' =  d_sttns[[i]]) 
}
  
d_MAE_ = c()
for(i in 1:length(pol3_list)){
  d_MAE_[i] = f_MAE(df_list[[i]])
}

df_MAE_ = data.frame(MAE = c(MAE_, d_MAE_),
                     serie = rep(c("original", "dif estacional"), each = length(MAE_)))

df_stats <- df_MAE_ %>%
  group_by(serie) %>%
  summarise(
    Q1 = quantile(MAE, 0.25),
    Q2 = median(MAE),
    Q3 = quantile(MAE, 0.75)
  )

ggplot(df_MAE_, aes(x = serie, y = MAE, fill = serie)) +
  geom_boxplot() +
  geom_text(data = df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(title = "Boxplots del MAE", y = expression(r[s])) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14) ) 

# Scatterplot MAE vs altitud ----
df_MAE_alt = data.frame(MAE = d_MAE_,
                     altitud = SRTM_30.sttns$SRTM_30_Col1,
                     geometry = geometry_df$geometry)

ggplot(df_MAE_alt, aes(x = MAE, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "Diagrama de dispersión  MAE vs altitud",
       x = "MAE",
       y = "Altitud") +
  theme_minimal()

# Boxplot por región natural ----
df_MAE_alt_g <- st_as_sf(df_MAE_alt$geometry, wkt = "geometry")
df_MAE_alt_g <- st_set_crs(df_MAE_alt_g, 4326)
df_MAE_alt_g <- st_transform(df_MAE_alt_g, st_crs(shp_depto))

df_MAE_alt_g <- st_join(df_MAE_alt_g, shp_depto["DeNombre"])
df_MAE_alt_g$ID = 1:nrow(df_MAE_alt_g)
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==2] = "San Andrés Providencia y Santa Catalina"
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==38] = "Atlántico"
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==91] = "La Guajira"
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==191] = "Arauca"
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==592] = "Vichada"
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==626] = "Putumayo"
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==636] = "Putumayo"
df_MAE_alt_g$DeNombre[df_MAE_alt_g$ID==274] = "Huila"

caribe_ = c("Atlántico", "Bolívar", "Cesar", "Córdoba" , "La Guajira", "Magdalena", 
            "San Andrés Providencia y Santa Catalina", "Sucre")
pacifica_ = c("Chocó", "Cauca", "Nariño", "Valle del Cauca")
orinoquia_ = c("Arauca", "Casanare", "Meta", "Vichada")
amazonia_ = c("Amazonas", "Caquetá", "Guainía", "Guaviare", "Putumayo", "Vaupés")
andina_ = c("Antioquia", "Boyacá", "Caldas", "Cundinamarca", "Huila", 
            "Norte de Santander", "Quindío", "Risaralda", "Santander", "Tolima")

df_MAE_alt_g$Region = ifelse(df_MAE_alt_g$DeNombre %in% caribe_, "Caribe",
                              ifelse(df_MAE_alt_g$DeNombre %in% pacifica_, "Pacífica",
                                     ifelse(df_MAE_alt_g$DeNombre %in% orinoquia_, "Orinoquía",
                                            ifelse(df_MAE_alt_g$DeNombre %in% amazonia_, "Amazonía",
                                                   ifelse(df_MAE_alt_g$DeNombre %in% andina_, "Andina", "")))))

table(df_MAE_alt_g$Region)

df_MAE_alt_g$MAE = df_MAE_alt$MAE
df_MAE_alt_g$altitud = df_MAE_alt$altitud

# ANOVA por Region ----
modelo <- aov(MAE ~ Region, data = df_MAE_alt_g)
modelo_ = summary(modelo)
p_value = modelo_[[1]]$`Pr(>F)`[1]
print(ifelse(p_value <= 0.05, "R. H0, concluye dif significativas entre al menos dos grupos",
             "NO R. H0, concluye no hay dif significativas"))

tukey <- TukeyHSD(modelo)
#par(mar = c(5, 10, 4, 2)) 
plot(tukey, las =1)

# validación de supuestos
## Normalidad
plot(density(modelo$residuals))
qqnorm(modelo$residuals)
qqline(modelo$residuals, col = "red")

lillie.test(modelo$residuals)
shapiro.test(modelo$residuals)

## Homocedasticidad
leveneTest(MAE ~ Region, data = df_MAE_alt_g)
bptest(modelo)

## Autocorrelación
dwtest(modelo)

# panel de scatterplot para cada categoria de región
ggplot(df_MAE_alt_g, aes(x = MAE, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión MAE vs altitud por región"),
       x = "MAE",
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ Region, scales = "free")

ggplot(df_MAE_alt_g) + 
  geom_boxplot(aes(x = Region, y = MAE)) +  # Usamos Región como eje X y r_s como eje Y
  theme_minimal() +  # Estilo de tema minimalista
  labs(x = "Región", y = "MAE", title = "Boxplot de MAE por Región natural") #+  # Etiquetas de los ejes y título
#theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar las etiquetas del eje X si es necesario


ggplot() + 
  geom_sf(data = shp_depto) +  # Muestra el mapa de los departamentos
  geom_sf(data = df_MAE_alt_g, aes(fill = MAE, color = MAE), size = 1) +
  scale_fill_gradient(low = "green", high = "red", name = "MAE") +# Muestra los puntos con valores de r_s
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de MAE", fill = "Valor de MAE") +  # Título y leyenda
  theme(legend.position = "right") 



# Boxplot por mes ----
d_MAE_ = list()
for(i in 1:length(pol3_list)){
  meses <- unique(pol3_list[[i]]$mes)
  #meses <- unique(pol3_list[[i]]$cat_mes)
  # Calculamos el coeficiente de correlación de Spearman para cada mes
  MAE <- sapply(meses, function(m){
    # Filtramos los datos para el mes 'm'
    data_mes <- subset(pol3_list[[i]], mes == m)
    #data_mes <- subset(pol3_list[[i]], cat_mes == m)
    # Calculamos el coeficiente de correlación de Spearman
    MAE_value <- f_MAE(data_mes)
    #MAE_value <- f_MAE(pol3_list[[i]])
    return(MAE_value)
  })
  # Almacenamos los resultados en el listado de R_Spearman
  d_MAE_[[i]] <- MAE
}

d_MAE_df <- as.data.frame(d_MAE_)
d_MAE_df <- t(d_MAE_df)
d_MAE_df <- as.data.frame(d_MAE_df)
colnames(d_MAE_df) <- unique(pol3_list[[1]]$mes)
#colnames(d_MAE_df) <- unique(pol3_list[[1]]$cat_mes)
rownames(d_MAE_df) <- 1:nrow(d_MAE_df)

d_MAE_long <- d_MAE_df %>%
  pivot_longer(cols = everything(), 
               #cols = c("DJF", "MAM", "JJA", "SON"), 
               names_to = "Mes", 
               #names_to = "Cat_mes", 
               values_to = "MAE")

d_MAE_long$Mes <- factor(d_MAE_long$Mes, 
                              levels = c("enero", "febrero", "marzo", "abril", "mayo", 
                                         "junio", "julio", "agosto", "septiembre", "octubre", 
                                         "noviembre", "diciembre"))

d_MAE_long$Cat_mes <- factor(d_MAE_long$Cat_mes, 
                             levels = c("DJF", "MAM", "JJA", "SON"))

ggplot(d_MAE_long, aes(x = Mes, y = MAE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = expression("Boxplot de  MAE por mes"),
       x = "Mes", 
       y = "MAE") #+

ggplot(d_MAE_long, aes(x = Cat_mes, y = MAE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = expression("Boxplot de MAE por grupo de meses"),
       x = "Meses", 
       y = "MAE") #+


#theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calcular los valores atípicos por mes
outliers_df <- d_MAE_long %>%
  group_by(Mes) %>%
  #group_by(Cat_mes) %>%
  summarise(
    Q1 = quantile(MAE, 0.25, na.rm = TRUE),
    Q3 = quantile(MAE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    outlier_count = sum(MAE < lower_bound | MAE > upper_bound, na.rm = TRUE),
    total_count = n(),
    outlier_pct = round((outlier_count / total_count) * 100, 1) # Porcentaje de atípicos
  )

# Unir el dataframe con los datos originales para graficar
d_MAE_long2 <- left_join(d_MAE_long, outliers_df, by = "Mes")
#d_MAE_long2 <- left_join(d_MAE_long, outliers_df, by = "Cat_mes")

ggplot(d_MAE_long2, aes(x = Mes, y = MAE)) +
#ggplot(d_MAE_long2, aes(x = Cat_mes, y = MAE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  #geom_hline(yintercept = 0.68, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(label = paste0(outlier_pct, "% \n outliers"), y = max(MAE, na.rm = TRUE) + 0.05), 
            size = 3, color = "black", fontface = "bold") +
  theme_minimal() +
  labs(title = expression("Boxplot de MAE por mes"),
  #labs(title = expression("Boxplot de MAE por grupo de meses"),
       x = "Meses", 
       y = "MAE") 


# panel de scatterplot para cada categoria de mes
df_MAE_alt_g2 = df_MAE_alt_g %>% bind_cols(d_MAE_df)

MAE_long <- df_MAE_alt_g2 %>%
  pivot_longer(cols = enero:diciembre,
               #cols = c('DJF', 'MAM', 'JJA', 'SON'),
               names_to = "Mes", 
               #names_to = "Cat_mes", 
               values_to = "MAE_") %>%
  select(ID, altitud, Mes, MAE_)
  #select(ID, altitud, Cat_mes, MAE_) 

MAE_long$Mes <- factor(MAE_long$Mes, 
                       levels = c("enero", "febrero", "marzo", "abril", "mayo", 
                                  "junio", "julio", "agosto", "septiembre", "octubre", 
                                  "noviembre", "diciembre"))

MAE_long$Cat_mes <- factor(MAE_long$Cat_mes, 
                           levels = c("DJF", "MAM", "JJA", "SON"))

ggplot(MAE_long, aes(x = MAE_, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión MAE vs altitud por mes"),
  #labs(title = expression("Diagrama de dispersión MAE vs altitud por grupo de meses"),
       x = "MAE",
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ Mes, scales = "free")
  #facet_wrap(~ Cat_mes, scales = "free")

quantile(MAE_)
mean(MAE_) # 69.07
# indica que, en promedio, la diferencia absoluta entre Chirps y las observaciones
## es de aproximadamente 69mm.
hist(MAE_)

# PBIAS ----
f_PBIAS = function(df){
  PBIAS = df %>%
    mutate(error_ = chirps.v2.0 - sttns) %>%
    summarise(PBIAS = (sum(error_, na.rm = TRUE)/sum(sttns, na.rm = TRUE))*100) %>%
    pull(PBIAS)
  return(PBIAS)
}

PBIAS_ = c()
for(i in 1:length(pol3_list)){
  PBIAS_[i] = f_PBIAS(pol3_list[[i]])
}

df_list = list()
for(i in 1:length(d_chirps.v2.0)){
  df_list[[i]] = data.frame('chirps.v2.0' = d_chirps.v2.0[[i]],
                            'sttns' =  d_sttns[[i]]) 
}

d_PBIAS_ = c()
for(i in 1:length(pol3_list)){
  d_PBIAS_[i] = f_PBIAS(df_list[[i]])
}

df_PBIAS_ = data.frame(PBIAS = c(PBIAS_, d_PBIAS_),
                      serie = rep(c("original", "dif estacional"), each = length(PBIAS_)))

df_stats <- df_PBIAS_ %>%
  group_by(serie) %>%
  summarise(
    Q1 = quantile(PBIAS, 0.25),
    Q2 = median(PBIAS),
    Q3 = quantile(PBIAS, 0.75)
  ) %>% mutate(IQR = Q3 - Q1,
               lim_inf = Q1 - 1.5 * IQR,
               lim_sup = Q3 + 1.5 * IQR)

ggplot(df_PBIAS_, aes(x = serie, y = PBIAS, fill = serie)) +
  geom_boxplot() +
  geom_text(data = df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(title = "Boxplots del PBIAS", y = 'PBIAS') +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14) ) +
  coord_cartesian(ylim = c(min(df_stats$lim_inf), max(df_stats$lim_sup)))

# Scatterplot MAE vs altitud ----
df_PBIAS_alt = data.frame(PBIAS = d_PBIAS_,
                          altitud = SRTM_30.sttns$SRTM_30_Col1,
                          geometry = geometry_df$geometry)

ggplot(df_PBIAS_alt, aes(x = PBIAS, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "Diagrama de dispersión  PBIAS vs altitud",
       x = "PBIAS",
       y = "Altitud") +
  theme_minimal()+
  coord_cartesian(xlim = c(min(df_stats$lim_inf), max(df_stats$lim_sup)))

# Boxplot por región natural ----
df_PBIAS_alt_g <- st_as_sf(df_PBIAS_alt$geometry, wkt = "geometry")
df_PBIAS_alt_g <- st_set_crs(df_PBIAS_alt_g, 4326)
df_PBIAS_alt_g <- st_transform(df_PBIAS_alt_g, st_crs(shp_depto))

df_PBIAS_alt_g <- st_join(df_PBIAS_alt_g, shp_depto["DeNombre"])
df_PBIAS_alt_g$ID = 1:nrow(df_PBIAS_alt_g)
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==2] = "San Andrés Providencia y Santa Catalina"
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==38] = "Atlántico"
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==91] = "La Guajira"
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==191] = "Arauca"
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==592] = "Vichada"
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==626] = "Putumayo"
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==636] = "Putumayo"
df_PBIAS_alt_g$DeNombre[df_PBIAS_alt_g$ID==274] = "Huila"

caribe_ = c("Atlántico", "Bolívar", "Cesar", "Córdoba" , "La Guajira", "Magdalena", 
            "San Andrés Providencia y Santa Catalina", "Sucre")
pacifica_ = c("Chocó", "Cauca", "Nariño", "Valle del Cauca")
orinoquia_ = c("Arauca", "Casanare", "Meta", "Vichada")
amazonia_ = c("Amazonas", "Caquetá", "Guainía", "Guaviare", "Putumayo", "Vaupés")
andina_ = c("Antioquia", "Boyacá", "Caldas", "Cundinamarca", "Huila", 
            "Norte de Santander", "Quindío", "Risaralda", "Santander", "Tolima")

df_PBIAS_alt_g$Region = ifelse(df_PBIAS_alt_g$DeNombre %in% caribe_, "Caribe",
                             ifelse(df_PBIAS_alt_g$DeNombre %in% pacifica_, "Pacífica",
                                    ifelse(df_PBIAS_alt_g$DeNombre %in% orinoquia_, "Orinoquía",
                                           ifelse(df_PBIAS_alt_g$DeNombre %in% amazonia_, "Amazonía",
                                                  ifelse(df_PBIAS_alt_g$DeNombre %in% andina_, "Andina", "")))))

table(df_PBIAS_alt_g$Region)

df_PBIAS_alt_g$PBIAS = df_PBIAS_alt$PBIAS
df_PBIAS_alt_g$altitud = df_PBIAS_alt$altitud

# panel de scatterplot para cada categoria de región
ggplot(df_PBIAS_alt_g, aes(x = PBIAS, y = altitud)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión PBIAS vs altitud por región"),
       x = "PBIAS",
       y = "Altitud") +
  theme_minimal()+
  coord_cartesian(xlim = c(min(df_stats$lim_inf), max(df_stats$lim_sup)))+
  facet_wrap(~ Region, scales = "free")

ggplot(df_PBIAS_alt_g) + 
  geom_boxplot(aes(x = Region, y = PBIAS)) +  # Usamos Región como eje X y r_s como eje Y
  theme_minimal() +  # Estilo de tema minimalista
  labs(x = "Región", y = "PBIAS", title = "Boxplot de PBIAS por Región natural")  +
  coord_cartesian(ylim = c(min(df_stats$lim_inf), max(df_stats$lim_sup)))#+  # Etiquetas de los ejes y título
#theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar las etiquetas del eje X si es necesario


ggplot() + 
  geom_sf(data = shp_depto) +  # Muestra el mapa de los departamentos
  geom_sf(data = df_PBIAS_alt_g %>% filter(PBIAS > min(df_stats$lim_inf) & 
                                             PBIAS < max(df_stats$lim_sup)), aes(fill = PBIAS, color = PBIAS), size = 1) +
  scale_fill_gradient(low = "green", high = "red", name = "PBIAS") +# Muestra los puntos con valores de r_s
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de PBIAS", fill = "Valor de PBIAS") +  # Título y leyenda
  theme(legend.position = "right") 

# Boxplot por mes ----
d_PBIAS_ = list()
for(i in 1:length(pol3_list)){
  meses <- unique(pol3_list[[i]]$mes)
  #meses <- unique(pol3_list[[i]]$cat_mes)
  # Calculamos el coeficiente de correlación de Spearman para cada mes
  PBIAS <- sapply(meses, function(m){
    # Filtramos los datos para el mes 'm'
    data_mes <- subset(pol3_list[[i]], mes == m)
    #data_mes <- subset(pol3_list[[i]], cat_mes == m)
    # Calculamos el coeficiente de correlación de Spearman
    PBIAS_value <- f_PBIAS(data_mes)
    #MAE_value <- f_MAE(pol3_list[[i]])
    return(PBIAS_value)
  })
  # Almacenamos los resultados en el listado de R_Spearman
  d_PBIAS_[[i]] <- PBIAS
}

d_PBIAS_df <- as.data.frame(d_PBIAS_)
d_PBIAS_df <- t(d_PBIAS_df)
d_PBIAS_df <- as.data.frame(d_PBIAS_df)
colnames(d_PBIAS_df) <- unique(pol3_list[[1]]$mes)
#colnames(d_MAE_df) <- unique(pol3_list[[1]]$cat_mes)
rownames(d_PBIAS_df) <- 1:nrow(d_PBIAS_df)

d_PBIAS_long <- d_PBIAS_df %>%
  pivot_longer(cols = everything(), 
               #cols = c("DJF", "MAM", "JJA", "SON"), 
               names_to = "Mes", 
               #names_to = "Cat_mes", 
               values_to = "PBIAS")

d_PBIAS_long$Mes <- factor(d_PBIAS_long$Mes, 
                         levels = c("enero", "febrero", "marzo", "abril", "mayo", 
                                    "junio", "julio", "agosto", "septiembre", "octubre", 
                                    "noviembre", "diciembre"))

d_MAE_long$Cat_mes <- factor(d_MAE_long$Cat_mes, 
                             levels = c("DJF", "MAM", "JJA", "SON"))

ggplot(d_MAE_long, aes(x = Mes, y = PBIAS)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = expression("Boxplot de  PBIAS por mes"),
       x = "Mes", 
       y = "PBIAS") +
  coord_cartesian(ylim = c(min(df_stats$lim_inf), max(df_stats$lim_sup))) #+

ggplot(d_MAE_long, aes(x = Cat_mes, y = MAE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = expression("Boxplot de MAE por grupo de meses"),
       x = "Meses", 
       y = "MAE") #+


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
for(i in 1:length(pol3_list)){
  NSE_[i] = f_NSE(pol3_list[[i]])
}

d_NSE_ = c()
for(i in 1:length(pol3_list)){
  d_NSE_[i] = f_NSE(df_list[[i]])
}

df_NSE_ = data.frame(NSE = c(NSE_, d_NSE_),
                     serie = rep(c("original", "dif estacional"), each = length(NSE_)))

df_stats <- df_NSE_ %>%
  group_by(serie) %>%
  summarise(
    Q1 = quantile(NSE, 0.25),
    Q2 = median(NSE),
    Q3 = quantile(NSE, 0.75)
  )

ggplot(df_NSE_, aes(x = serie, y = NSE, fill = serie)) +
  geom_boxplot() +
  geom_text(data = df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(title = "Boxplots del NSE", y = expression(r[s])) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14) ) 


quantile(NSE_)
mean(NSE_) # 0.33
# indica que, Chirps tiene una buena predicción, mejor que usar 
## la media de los valores observados como predicción, sin embargo
### el promedio está lejos de 1, aún hay un margen significativo para mejorar
hist(NSE_)

# Check point 2 ----
#save.image('df_eval_metrics_ADF.RData')
load('df_eval_metrics_ADF.RData')
#load('df_eval_metrics.RData')

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

