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
m_order4 = ifelse(m_order3==-1,12, m_order4)
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
nextodd <- function(x) {
  x <- round(x)
  if (x%%2 == 0) 
    x <- x + 1
  as.integer(x)}

nextodd(ceiling(1.5 * frequency(y)/(1 - 1.5/(10 * as.integer(length(y)) + 1))))
fit1 <- stl(y, s.window = "periodic")
# t.window = 13 (fpp2)
fit11<- stl(y, t.window = 13, s.window = "periodic")
# Selecting a shorter trend window, avoid produces a trend-cycle component that is too rigid
fit12<- stl(y, t.window = 7, s.window = "periodic")

stl_chirps.v2.0 = list()
for(i in 1:length(pol3_list)){
  y <- tsibble::tsibble(pol3_list[[i]] %>% mutate(Date = as.Date(str_c(Date, '.01'), format = '%Y.%m.%d')) 
                        %>%  select(c("Date" ,  "chirps.v2.0")), index = Date)
  
  y = y %>% mutate(Date = yearmonth(Date)) %>% 
    as_tsibble(index = Date)
  fit1 <- stl(y, s.window = "periodic")
  #fit11<- stl(y, t.window = 13, s.window = "periodic")
  #fit12<- stl(y, t.window = 7, s.window = "periodic")
  stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit1$time.series[, "seasonal"]
  #stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit11$time.series[, "seasonal"]
  #stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit12$time.series[, "seasonal"]
}

plot(stl_chirps.v2.0[[s]])

m_order5 = c()
k = 1
for(i in 1:length(pol3_list)){
  m_order5[k] <- tryCatch(
    {
      # Intenta ajustar el modelo AR
      #ar_obj <- ar(diff(stl_chirps.v2.0[[i]]), method = 'mle')
      ar_obj <- ar(stl_chirps.v2.0[[i]], method = 'mle')
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
m_order5 = ifelse(m_order5==0,12, m_order5)

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
# The default setting for monthly data is trend(window=21)
fit2 = y |>
  model(
    STL(`chirps.v2.0` ~ 
          season(window = "periodic"),
        robust = TRUE)) |>
  components() 

# t.window = 13 (fpp2)
fit21 = y |>
  model(
    STL(`chirps.v2.0` ~ trend(window = 13) + 
          season(window = "periodic"),
        robust = TRUE)) |>
  components() 

# Selecting a shorter trend window, avoid produces a trend-cycle component that is too rigid
fit22 = y |>
  model(
    STL(`chirps.v2.0` ~ trend(window = 7) + 
          season(window = "periodic"),
        robust = TRUE)) |>
  components() 

stl_chirps.v2.0 = list()
for(i in 1:length(pol3_list)){
  y <- tsibble::tsibble(pol3_list[[i]] %>% mutate(Date = as.Date(str_c(Date, '.01'), format = '%Y.%m.%d')) 
                        %>%  select(c("Date" ,  "chirps.v2.0")), index = Date)
  
  y = y %>% mutate(Date = yearmonth(Date)) %>% 
    as_tsibble(index = Date)
  
  fit2 = y |>
    model(
      #STL(`chirps.v2.0` ~ trend(window = 13) +
      STL(`chirps.v2.0` ~ 
            season(window = "periodic"),
          robust = TRUE)) |>
    components() 
  stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit2$season_year
}

head(y$chirps.v2.0 - fit2$season_year)
head(fit2$season_adjust)

m_order6 = c()
k = 1
for(i in 1:length(pol3_list)){
  m_order6[k] <- tryCatch(
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

table(m_order6)
m_order6 = ifelse(m_order6==-1,12, m_order6)

## prueba ADF librería aTSA ----
len_ = length(pol3_list)
df_ADF4 = data.frame(lag_t1 = rep(0, len_), p.value_t1 = rep(0, len_),
                     lag_t2 = rep(0, len_), p.value_t2 = rep(0, len_),
                     lag_t3 = rep(0, len_), p.value_t3 = rep(0, len_))

k = 1
for(i in 1:length(pol3_list)){
  l_ADF = aTSA::adf.test(stl_chirps.v2.0[[i]], nlag = m_order6[i], output = FALSE)
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
R_pearson = c()
for(i in 1:length(pol2_list)){
  R_pearson[i] = cor(pol2_list[[i]]$chirps.v2.0, pol2_list[[i]]$sttns, use = "complete.obs")
}

R_Spearman = c()
for(i in 1:length(pol3_list)){
  R_Spearman[i] = cor(pol3_list[[i]]$chirps.v2.0, pol3_list[[i]]$sttns, use = "complete.obs", method = "spearman")
}

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
#save.image('df_eval_metrics_ADF.RData')
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

