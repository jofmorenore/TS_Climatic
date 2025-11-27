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
library(viridis) # colores(?)
library(nortest) # pruebas de hipótesis ANOVA
library(car) # pruebas de hipótesis ANOVA
library(lmtest) # pruebas de hipótesis ANOVA
library(qmap) # correción de sesgo
library(hydroGOF) # medidas de desempeño
library(rlang) # strings como symbols para función sym
library(scales) # agrega simbolo % al eje y
library(missForest)
library(imputeTS) # na_kalman, na_seadec 
library(units) # convertir a km
library(future) # procesamiento en paralelo
library(future.apply)


setwd('~/TS_climatic/')
list.files()

# 1. Carga los datos ----
#sttns = arrow::read_parquet("sttns_pcp_col.parquet")
#dim(sttns)

#pcp_col = openxlsx::read.xlsx("pcp_col.xlsx")
pcp_col = arrow::read_parquet("./IDEAM/sttns_pcp_col_train.parquet")
dim(pcp_col)
#class(pcp_col)

chirps_df = arrow::read_parquet("chirps_df.parquet")

# 2. Convierte a objeto sf ----
chirps.point <- st_as_sf(x = chirps_df, 
                         coords = c("x", "y"),
                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

## remueve duplicados por LON y LAT (2 sttns; Valle y Meta)
sttns = dplyr::distinct(pcp_col, `LONGITUD` , `LATITUD`, .keep_all = TRUE)

sttns.points <- st_as_sf(x = sttns, 
                         coords = c("LONGITUD", "LATITUD"),
                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 

#terra::crs(sttns.points)
print(sttns.points)

# 3. Función para extraer los pixeles chirps más cercanos ----
#chirps.point#[1]

sttns.points$n_chirps = 1 # (oct/2024)
sttns.points$r = 0.05
table(sttns.points$r)

pol_ = function(x, r, CHIRPS = chirps.point){
  # Argumentos: ----------------------------------------------------------------
  # 
  # x df 1 x ncol datos de entrada de c/sttn
  # r num col radio (default 0.05)
  # CHIRPS sf con los pxls de CHIRPS
  #----------------------------------------------------------------------------
  
  #En caso de error con la función correr manualmente --------------------------
  # set.seed(25052024)
  # i = sample(1:nrow(sttns.points), size = 1)
  # print(x)
  #x = sttns.points[i,]
  #r = sttns.points$r[i]
  #----------------------------------------------------------------------------
  
  # 1. crea un df 2(lon, lat) x 4(Norte, Sur,...) con las coordenadas ---------
  pol = data.frame(matrix(c(x$geometry[[1]] +c(-r,-r),
                            x$geometry[[1]] +c(-r,r),
                            x$geometry[[1]] +c(r,r),
                            x$geometry[[1]] +c(r,-r)
                            ), nc = 2, byrow = F
                          )
                   )
  names(pol) = c("X", "Y")
  
  # 2. convierte a objeto sf y conforma un polígono ---------------------------
  pol = pol %>%
    st_as_sf(coords = c("X", "Y"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON")                            
  
  # 3. join espacial ----------------------------------------------------------
  ## join entre el polígono y los pixeles de CHIRPS
  join_chirps_pol = st_join(pol, CHIRPS)
  ## join para quedarse con los pxls CHIRPS contenidos en el polígono
  chirps_ = st_join(CHIRPS, join_chirps_pol) %>% 
    na.omit() %>% distinct(geometry, .keep_all = T)
  
  # 4. Calcula el pxl CHIRP más cercano a la sttn (new) -----------------------
  ## !!! st_nearest_feature: toma mucho tiempo en un loop
  closest_idx <- st_nearest_feature(x, chirps_)
  chirps_closest <- chirps_[closest_idx[1], 1:ncol(CHIRPS)]
  colnames(chirps_closest) = gsub(".x","",colnames(chirps_closest))
  
  # Calcula distancias (old) --------------------------------------------------
  #distances <- st_distance(x, chirps_)
  #s_distances <- as.vector(distances)
  #closest <- order(s_distances)[1] 
  #geom_closest = chirps_$geometry[closest]
  #chirps_closest = chirps_[chirps_$geometry %in% geom_closest, 1:ncol(CHIRPS)]
  #
  #chirps_closest = chirps_[chirps_$geometry %in% chirps_$geometry[closest],1:505] %>% 
  #  bind_cols(select(as.data.frame(pt), n_chirps))
  #colnames(chirps_closest) = gsub(".x","",colnames(chirps_closest))
  
  # gráfico para ilustrar lo que devuelve la función --------------------------
  #ggplot() + 
  #  geom_sf(data = chirps_, aes(fill = `chirps-v2.0.1981.01.x`, col = `chirps-v2.0.1981.01.x`)) +
  #  geom_sf(data = x, fill = NA, col = 'red', lwd = 0.2) +
  #  geom_sf(data = chirps_closest, fill = NA, col = 'green', lwd = 0.2 )
  
  return(chirps_closest)
}

pol_list = list()
for(i in 1:nrow(sttns.points)){
  pol_list[[i]]  = pol_(sttns.points[i,], sttns.points$r[i])
  pol_list[[i]]$ID =  sttns.points$CodigoEstacion[i]
}

## número de pxls para c/sttn
n_pxls = c()
for(i in 1:length(pol_list)){
  n_pxls[i] = nrow(pol_list[[i]])
}
table(n_pxls)

## conteo de nulos (CHIRPPS no tiene nulos)
n_na = c()
for(i in 1:length(pol_list)){
  n_na[i] = sum(is.na(pol_list[[i]]))
}
table(n_na) 
v_na = which(n_na!=0) 

## elimina elementos de la lista que tienen nulos 
pol_list = pol_list[-v_na]


# Check point 1 ----
save.image('df_pol_chirps_25.RData')
load('df_pol_chirps_25.RData')

# 4. Carga datos elevación de datos abiertos ----
SRTM_30 <- rast("Servicio-159/SRTM30/SRTM_30_Col1.tif")
print(SRTM_30)

SRTM_30 <- project(SRTM_30, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
crs(SRTM_30) 

# 5. Extrae la altitud de los puntos de las estaciones ----
SRTM_30.sttns <- terra::extract(SRTM_30, sttns.points, ID = FALSE)

SRTM_30.sttns =  SRTM_30.sttns %>% 
  mutate(CodigoEstacion = sttns.points$CodigoEstacion) %>%  
  filter(!is.na(SRTM_30_Col1)) %>% 
  rename(SRTM_30 = SRTM_30_Col1)

# 6. join con sttns.points por ID ----
## join por ID para traer las cols de sttns (Obs.)
pol2_list = list()

# crea vector con nombres de cols a eliminar (fuera del hist de 40 años )
test_cols = paste0("chirps-v2.0.",
                   format(seq(from = as.Date("2021-01-01"), 
                              to = as.Date("2022-12-01"), by = "month"), 
                          '%Y.%m'))

for(i in 1:length(pol_list)){
  # (1/oct) se incluye filter para homogeneizar el rango de tiempo
  # se requerirá CHIRPS después del 2020-12  para pronóstico ?
  pol2_list[[i]] = left_join(as.data.frame(pol_list[[i]]) %>% dplyr::select(-all_of(test_cols)), 
                             sttns.points %>% 
                               rename('ID' = CodigoEstacion) %>% 
                               dplyr::select(-c('n_chirps', 'r')), 
                             by = 'ID')
  names(pol2_list[[i]]) = gsub('\\.','-',names(pol2_list[[i]]))
}
# OBS: las sttns van hasta 2020-12, mientras que CHIRPS va hasta 2022

# 7. Pivotar de columnas a filas ----

for(i in 1:length(pol2_list)){
  pol2_list[[i]] = pol2_list[[i]] %>%
    # selecciona sólo las cols de CHIRPS
    # pivota cols a filas
    dplyr::select(-matches("^(19|20|geometry-y)")) %>%
    pivot_longer(
      cols = starts_with("chirps-v2-0-"),     
      names_to = "Date",
      values_to = "chirps-v2-0"
    ) %>% mutate(Date = paste0(gsub('chirps-v2-0-','', Date), '-01')) %>% 
    left_join(
      pol2_list[[i]] %>% 
        dplyr::select(-matches("^(chirps-v2-0-|geometry-x)")) %>% 
        pivot_longer(
          cols = matches("^(19|20)"),
          names_to = "Date",
          values_to = "sttns"
        )%>% distinct(), by = c('Date', 'ID', 'NombreEstacion',
                                'DEPARTAMENTO', 'n', 'ALTITUD')
    ) 
}

# 8. Evalua la estacionariedad de las series ----
# La estacionariedad de una serie temporal implica que sus propiedades estadísticas
## no cambian con el tiempo, en otras palabras no hay tendencia.

## 8.1 get_ar_order: función obtener orden modelo AR (arg lag de prueba ADF) ----
# crear el ciclo for para calcular los ordenes de cada serie CHIRPS
# Tsay recomienda diferenciar ordinalmente la serie para ajustar el modelo AR

get_ar_order <- function(x, diff_ = TRUE){
  # Argumentos: ----------------------------------------------------------------
  # 
  # x lista con dfs con col num 'sttns' 1 x n series de tiempo univariadas
  # diff boolean deafult TRUE aplica diferenciación ordinaria a x (convierte a estacionaria)
  #----------------------------------------------------------------------------
  
  # len_ num iteraciones (elementos de la lista)
  len_ = length(x)
  k = 1
  m_order = c()
  # cond if/else especifica bucle for por separado para no evaluar en c/iter el arg diff_
  
  if(diff_ == TRUE){
    for(i in 1:len_){
      m_order[k] <- tryCatch(
        {
          # Intenta ajustar el modelo AR
          # para la ts diff estacionalmente crear nueva lista con col `chirps-v2-0`
          ar_obj <- ar(diff(x[[i]]$`chirps-v2-0`), method = 'mle') # , na.action = na.approx
          # sin na.action arroja error porque no puede ajustar AR en TS con NAs
          # na.aprox replace NA by interpolation
          # ar_obj <- ar(diff(x[[i]]$sttns), method = 'mle', na.action = na.approx) 
          ar_obj$order  # devuelve el orden del modelo ajustado
        },
        error = function(e){
          # Si hay un error, devuelve -1
          -1
          }
      )
      # Incrementa el índice k
      k <- k + 1
    }
  }else{
    # si diff_ = FALSE no aplica diferencia ordinaria
    for(i in 1:len_){
      m_order[k] <- tryCatch(
        {
          # Intenta ajustar el modelo AR sin diferenciar ordinaria 
          ar_obj <- ar(x[[i]]$`chirps-v2-0`, method = 'mle') #, na.action = na.approx
          # ar_obj <- ar(x[[i]]$sttns, method = 'mle', na.action = na.approx) 
          ar_obj$order  # devuelve el orden del modelo ajustado
        },
        error = function(e) {
          # Si hay un error, devuelve -1
          -1
        }
      )
      # Incrementa el índice k
      k <- k + 1
    }
  }
  return(m_order)
}

## 8.2 ADF_test: función prueba aTSA::adf.test (Raiz Unitaria) ----

ADF_test <- function(x, m_order, type = 'type1'){
  # Argumentos: ----------------------------------------------------------------
  # 
  # x lista con dfs con col num 'chirps-v2-0' 1 x n ts CHIRPS univariadas
  # m_order vector num ordenes modelo AR (arg lag para prueba ADF)
  # type str col de la salida aTSA::adf.test default type1 (no drift no trend)
  #----------------------------------------------------------------------------
  
  if(!(type %in% c("type1", "type2", "type3") ) ) stop('parameter type is not one of type1, type2, or type3')
  # si no obtiene orden del modelo AR reemplazar por orden 12 (default ts mensuales)
  m_order = ifelse(m_order==-1, 12, m_order)
  # num elementos de la lista para iter bucle for
  len_ = length(x)
  
  # pt1 type1 ----
  if(type == 'type1'){
    
    df_ADF = data.frame(p.value_t1 = rep(0, len_), lag_t1 = rep(0, len_))  
    for(i in 1:len_){
      # para la ts diff estacionalmente crear nueva lista con col `chirps-v2-0`
      l_ADF = aTSA::adf.test(x[[i]]$`chirps-v2-0`, 
                             nlag = m_order[i], output = FALSE)
      # valor máx: rezago donde posiblemente No R. H0
      # por qué toma p.value como el valor max.? no debería ser el del rezago más alto?
      df_ADF$p.value_t1[i] = max(l_ADF$type1[,'p.value'])
      df_ADF$lag_t1[i] = max(which(l_ADF$type1[,'p.value'] == df_ADF$p.value_t1[i]))
      df_ADF$RD_t1[i] = ifelse(df_ADF$p.value_t1[i]< 0.05, 'R H0: estacionaria', 
                               'No R. H0: Raíz unitaria (no-estacionaria)')
    }
  } else if(type == 'type2'){
    df_ADF = data.frame(lag_t2 = rep(0, len_), p.value_t2 = rep(0, len_))  
    for(i in 1:len_){
      # para la ts diff estacionalmente crear nueva lista con col `chirps-v2-0`
      l_ADF = aTSA::adf.test(x[[i]]$`chirps-v2-0`, 
                             nlag = m_order[i], output = FALSE)
      df_ADF$p.value_t2[i] = max(l_ADF$type2[,'p.value'])
      df_ADF$lag_t2[i] = which(l_ADF$type2[,'p.value'] == df_ADF$p.value_t2[i])
      df_ADF$RD_t2[i] = ifelse(df_ADF$p.value_t2[i]< 0.05, 'R H0: estacionaria', 
                               'No R. H0: Raíz unitaria (no-estacionaria)')
    }
  } else{
    df_ADF = data.frame(lag_t3 = rep(0, len_), p.value_t3 = rep(0, len_))  
    for(i in 1:len_){
      # para la ts diff estacionalmente crear nueva lista con col `chirps-v2-0`
      l_ADF = aTSA::adf.test(x[[i]]$`chirps-v2-0`, 
                             nlag = m_order[i], output = FALSE)
      df_ADF$p.value_t3[i] = max(l_ADF$type3[,'p.value'])
      df_ADF$lag_t3[i] = which(l_ADF$type3[,'p.value'] == df_ADF$p.value_t3[i])
      df_ADF$RD_t3[i] = ifelse(df_ADF$p.value_t3[i]< 0.05, 'R H0: estacionaria', 
                               'No R. H0: Raíz unitaria (no-estacionaria)')
    }
  }
  return(df_ADF)
}

# Evalua estacionalidad a ts CHIRPS diff estacional (R. H0 - concluye estacionaria) ----
## motivación: para la ts CHIRPS concluye no estacionariedad, tal vez debido a ciclos (ocultos) deterministas
## luego una forma de chequear si la ts deja de ser no-estacionaria es eliminar la comp. estacional
## aplicando la prueba ADF a una ts transformada (desestacionalizada) de la ts CHIRPS

## 8.3 df_pgram: periodograma ----
## herramienta para detección de ciclos (ocultos) deterministas
len_ = length(pol2_list)
df_pgram = data.frame(freq = rep(0, len_), period = rep(0, len_))

# encontrar el periodo del ciclo estacional
for(i in 1:len_){
  pgram = spectrum(pol2_list[[i]]$`chirps-v2-0`,log='no')
  f_max = which.max(pgram$spec)
  df_pgram$freq[i] = pgram$freq[f_max]
  df_pgram$period[i] = round(1/pgram$freq[f_max])
}

#spectrum(pol3_list[[91]]$chirps.v2.0,log='no')
#abline(v = 1/12, col = 'red')
## plot muestra pixeles con periodos 4 y 6 ----
s_p4 = which(round(df_pgram$period)==4)
p6 = which(round(df_pgram$period)==6)
s_p480 = which(round(df_pgram$period)==480)
set.seed(22112024)
s_p6  = sample(p6, size = 1, replace = FALSE)


spectrum(pol2_list[[s_p4]]$`chirps-v2-0`,log='no')
abline(v = 1/12, col = 'red')
abline(v = 1/6, col = 'blue')
abline(v = 1/4, col = 'green')

y <- tsibble::tsibble(pol2_list[[s_p4]] %>% mutate(Date = as.Date(Date)) 
                      %>% dplyr::select(c("Date" ,  "chirps-v2-0")), index = Date)

y = y %>% mutate(Date = yearmonth(Date)) %>% 
  as_tsibble(index = Date)

gg_season(y, `chirps-v2-0`)
gg_subseries(y, `chirps-v2-0`)

# por qué hay 1 de periodo 480 (valor spectrum más alto cercano a cero)
table(round(df_pgram$period))

## 8.4 d_chirps.v2.0: lista ts CHIRPS diff estacionalmente ----
# cambiar s_p480 a 1/12
df_pgram$period[df_pgram$period==480] = 12

d_chirps.v2.0 = list()
for(i in 1:len_){
  d_chirps.v2.0[[i]] = data.frame(
    chirps_col = diff(pol2_list[[i]]$`chirps-v2-0`, 
                      lag = df_pgram$period[i])) %>% 
    rename(`chirps-v2-0` = chirps_col)
}

for(i in 1:len_){
  pol2_list[[i]]$d_chirps.v2.0 = c(rep(NA, df_pgram$period[i]),
                                   d_chirps.v2.0[[i]]$`chirps-v2-0`)
}

## 8.5 obtiene orden AR usando ts CHIRPS diff estacionalmente ----
## con diferencia ordinaria (diff_ = TRUE)
m_order3 = get_ar_order(x = d_chirps.v2.0, diff_ = TRUE)
table(m_order3)

## sin diferencia ordinaria (diff_ = FALSE)
m_order4 = get_ar_order(x = d_chirps.v2.0, diff_ = FALSE)
table(m_order4)

## 8.6 prueba aTSA::adf.test a ts CHIRPS diff estacionalmente ----
## con m_order3 (arg lag prueba ADF) orden AR obtenido diff ordinaria a la ts CHIRPS diff estacionalmente
## type1: no drift and no trend
df_ADF3 = ADF_test(d_chirps.v2.0, m_order3, type = 'type1')

# si se quiere comparar con los otros args type de la función aTSA::adf.test
#df_ADF3 = ADF_test(d_chirps.v2.0, m_order3, type = 'type1') %>% 
#  # type2: drift but no trend
#  bind_cols(ADF_test(d_chirps.v2.0, m_order3, type = 'type2')) %>% 
#  # type3: drift and linear trend 
#  bind_cols(ADF_test(d_chirps.v2.0, m_order3, type = 'type3'))

table(df_ADF3$RD_t1) # R. H0: concluye no hay Raíz unitaria (estacionaria?)

## con m_order4 (arg lag prueba ADF) orden AR obtenido sin diff ordinalmente a la ts CHIRPS diff estacionalmente
df_ADF4 = ADF_test(d_chirps.v2.0, m_order4, type = 'type1')
table(df_ADF4$RD_t1) # R. H0: concluye no hay Raíz unitaria (estacionaria?)

## 8.7 prueba urca::ur.df a ts CHIRPS diff estacionalmente ----
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

df_ur3 = data.frame(RD_t1 = rep('', len_))
for(i in 1:len_){
  ADF_ur = urca::ur.df(d_chirps.v2.0[[i]]$`chirps-v2-0`, type  = "none", 
                       lags = m_order3[i], selectlags = "AIC")
                       #lags = m_order4[i], selectlags = "AIC")
  
  df_ur3$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur3$RD_t1) # R. H0: concluye no hay Raíz unitaria (estacionaria?)
#table(df_ur3$RD_t1) # R. H0: concluye no hay Raíz unitaria (estacionaria)

# 9 .Calculas medidas de desempeño ----

## 9.1 crea cols mes y cat_mes ----
for(i in 1:len_){
  pol2_list[[i]]$mes = format(as.Date(pol2_list[[i]]$Date), "%B")
}

meses_ = unique(pol2_list[[i]]$mes)
DJF_ = meses_[c(12,1,2)] # c("diciembre", "enero", "febrero" )
MAM_ = meses_[3:5] # c("marzo", "abril", "mayo")
JJA_ = meses_[6:8] # c("junio", "julio", "agosto")
SON_ = meses_[9:11] # c("septiembre", "octubre", "noviembre")

for(i in 1:len_){
  ## crear col cat (DJF/MAM/JJA/SON) a partir de la variable mes 
  pol2_list[[i]]$cat_mes <-
    ifelse(pol2_list[[i]]$mes %in% DJF_, 'DJF', 
           ifelse(pol2_list[[i]]$mes %in% MAM_, 'MAM', 
                  ifelse(pol2_list[[i]]$mes %in% JJA_, 'JJA', 
                         ifelse(pol2_list[[i]]$mes %in% SON_, 'SON', NA))))
}

## 9.2 d_sttns: lista con dfs de ts sttns diff estacionalmente ----
## misma transformación que en 8.4 para ts CHIRPS diff estacionalmente
## calcula R_s con ts diff estacionalmente (juzgadas como estacionarias por prueba ADF)
d_sttns = list()
for(i in 1:len_){
  # aplica la misma diff estacional a partir del periodo obtenido con ts CHIRPS
  # por criterio de misma transformación
  d_sttns[[i]] = data.frame(sttns = diff(pol2_list[[i]]$sttns, 
                                         lag = df_pgram$period[i]))
}

## crea la col para poder calcular corr por mes, cat_mes y región_natural
for(i in 1:len_){
  ## crea col d_sttns en lista pol2_list
  pol2_list[[i]]$d_sttns = c(rep(NA, round(df_pgram$period[i])),
                             d_sttns[[i]]$sttns)
}

## 9.3 df_corr: data frame para graficar coef corr ----

ADZ = "Archipielago De San Andres, Providencia Y Santa Catalina"
df_corr = data.frame(CodigoEstacion = 
                       sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ])

for(i in 1:len_){
  ## R_pearson: vector num con coef corr Pearson 
  df_corr$R_pearson[i] = cor(pol2_list[[i]]$`chirps-v2-0`, pol2_list[[i]]$sttns, 
                             use = "complete.obs")
  ## R_Spearman: vector num con coef corr Spearman a ts original 
  df_corr$R_Spearman[i] = cor(pol2_list[[i]]$`chirps-v2-0`, pol2_list[[i]]$sttns, 
                              use = "complete.obs", method = "spearman")
  ## d_R_pearson: vector num con coef corr Pearson a ts diff estacionalmente 
  df_corr$d_R_pearson[i] = cor(d_chirps.v2.0[[i]]$`chirps-v2-0`, d_sttns[[i]]$sttns, 
                               use = "complete.obs")
  ## d_R_Spearman: vector num con coef corr Spearman a ts diff estacionalmente 
  df_corr$d_R_Spearman[i] = cor(d_chirps.v2.0[[i]]$`chirps-v2-0`, d_sttns[[i]]$sttns, 
                                use = "complete.obs", method = "spearman")
}

## Boxplot comparación Corr Pearson ts diff estacionalmente vs ts original ----
df_pearson = df_corr %>% dplyr::select(c(R_pearson, d_R_pearson)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "r") %>% 
  mutate(serie = factor(serie, labels = c("diff estacionalmente", 
                                          "original")))

bx_text = function(df, serie, r){
  # Argumentos: ----------------------------------------------------------------
  # 
  # df data frame con cols serie y r
  # serie str col para agrupar 
  # r num col para resumir 
  #----------------------------------------------------------------------------
  output = list()
  
  output$df_stats = df %>% 
    group_by({{serie}}) %>% 
    summarise(
      Q1 = quantile({{r}}, 0.25), Q2 = median({{r}}), Q3 = quantile({{r}}, 0.75)
    )
  
  output$df_outliers = df %>%
    group_by({{serie}}) %>%
    summarise(
      min_r = min({{r}}) - 0.1,
      max_r = max({{r}}) - 0.1,
      Q1 = quantile({{r}}, 0.25, na.rm = TRUE),
      Q3 = quantile({{r}}, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower_bound = Q1 - 1.5 * IQR,
      upper_bound = Q3 + 1.5 * IQR,
      outlier_count = sum({{r}} < lower_bound | {{r}} > upper_bound, na.rm = TRUE),
      total_count = n(),
      outlier_pct = round((outlier_count / total_count) * 100, 1) # Porcentaje de atípicos
    )
  return(output)
}

bx_ = bx_text(df_pearson, serie, r)

## Al comparar graficamente el coef corr Pearson de la ts original vs ts diff estacionalmente
## +obs1: se observa que r_pearson para la ts original está inflado y sobreestima el valor de corr
## +obs2: el % de outliers disminuye al transformar diff estacionalmente

ggplot(df_pearson, aes(x = fct_rev(serie), y = r, fill = serie)) +
  geom_boxplot() + 
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del coef corr Pearson", x = "serie", y = expression(r)) +
  theme_minimal() +
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14,face="bold"))

## Boxplot coef corr Spearman ----
# df para graficar boxpot coef corr Spearman ts original vs ts transformada
df_Spearman = df_corr %>% dplyr::select(c(R_Spearman, d_R_Spearman)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "r_s") %>% 
  mutate(serie = factor(serie, labels = c("diff estacionalmente", "original")))

quantile(df_corr$d_R_Spearman)

bx_ = bx_text(df_Spearman, serie, r_s)

## +obs1: se observa que r_Spearman para la ts original está inflado y sobreestima el valor de r_s
## +obs2: el % de outliers disminuye al transformar diff estacionalmente
## +obs3: El RIQ es mayor en la ts diff estacionalmente, la variabilidad de r_s aumenta (amplitud de la caja)
## +obs4: menos % outliers (1.9% a 2.2%) en ts diff estacionalmente que en corr_Pearson

ggplot(df_Spearman, aes(x = fct_rev(serie), y = r_s, fill = serie)) +
  geom_boxplot() +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del coef de Spearman", x = "serie", y = expression(r[s])) +
  theme_minimal() +
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## 9.4 test de significancia coef corr spearman a ts diff estacionalmente ----
test_d_R_Spearman = c()
for(i in 1:len_){
  obj_ = cor.test(d_chirps.v2.0[[i]]$`chirps-v2-0`, d_sttns[[i]]$sttns, 
                  method = "spearman", exact = FALSE)
  test_d_R_Spearman[i] = obj_$p.value
}

# todos significativos a un nivel de confianza del 95%
max(test_d_R_Spearman) # max es 0.03 < 0.05 (con lo cual R H0 y conlcuye R_s != 0 por tanto R_s significativos)
# 1 valor de rho 0.10 con p-value 0.03
df_corr$d_R_Spearman[which(test_d_R_Spearman>=0.01)] # valor de rho
test_d_R_Spearman[which(test_d_R_Spearman>=0.01)] # valor de los p-value

## Boxplot R_Spearman por mes -------------------------------------------------
## transformaciones para poder graficas boxplot
# R_Spearman_df: df c/col es un mes
meses <- unique(pol2_list[[i]]$mes)
df_R_Spearman_mes = data.frame(matrix(numeric(len_), nc = 12, nr = len_)) %>%
  `colnames<-`(meses)

for(i in 1:len_){
  # Calculamos el coef corr Spearman para c/mes
  corrs <- sapply(meses, function(m){
    # Filtra datos para el mes 'm'
    df_mes <- subset(pol2_list[[i]], mes == m)
    # Calcula coef corr Spearman
    cor_value <- cor(df_mes$d_chirps.v2.0, df_mes$d_sttns, 
                     use = "complete.obs", method = "spearman")
    return(cor_value)
  })
  # Almacenamos los resultados en el listado de R_Spearman
  df_R_Spearman_mes[i,]  <- corrs
}

df_R_Spearman_mes$CodigoEstacion = sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ]

df_R_Spearman_mes <- df_R_Spearman_mes %>%
  pivot_longer(cols = all_of(meses), 
               names_to = "mes", 
               values_to = "r_s")

# factor ordenado especifica levels
df_R_Spearman_mes = df_R_Spearman_mes %>% 
  mutate(mes = factor(mes, levels = meses))

bx_ = bx_text(df_R_Spearman_mes, serie = mes, r = r_s)

ggplot(df_R_Spearman_mes, aes(x = mes, y = r_s)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  theme_minimal() +
  labs(title = expression("Boxplot del coef corr Spearman r"["s"] ~ "por mes"),
       x = "Mes", 
       y = expression(r["s"])) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Boxplot R_Spearman por cat_mes (DJF, MAM, JJA, SON) ----
cat_meses <- unique(pol2_list[[i]]$cat_mes)

df_R_Spearman_cat_mes = data.frame(matrix(numeric(len_), nc = 4, nr = len_)) %>%
  `colnames<-`(cat_meses)

for(i in 1:len_){
  # Calculamos el coeficiente de correlación de Spearman para cada mes
  corrs <- sapply(cat_meses, function(m){
    # Filtramos los datos para el mes 'm'
    data_mes <- subset(pol2_list[[i]], cat_mes == m)
    # Calculamos el coeficiente de correlación de Spearman
    cor_value <- cor(data_mes$d_chirps.v2.0, data_mes$d_sttns, use = "complete.obs", method = "spearman")
    return(cor_value)
  })
  # Almacenamos los resultados en el listado de R_Spearman
  df_R_Spearman_cat_mes[i,] <- corrs
  #R_Spearman_mes[[i]] <- corrs
}

df_R_Spearman_cat_mes$CodigoEstacion = sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ]

df_R_Spearman_cat_mes <- df_R_Spearman_cat_mes %>%
  pivot_longer(cols = all_of(cat_meses), 
               names_to = "cat_mes", 
               values_to = "r_s")

df_R_Spearman_cat_mes = mutate(df_R_Spearman_cat_mes, 
                               cat_mes = factor(cat_mes, levels = cat_meses))

bx_ = bx_text(df_R_Spearman_cat_mes, cat_mes, r_s)

ggplot(df_R_Spearman_cat_mes, aes(x = cat_mes, y = r_s)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  theme_minimal() +
  labs(title = expression("Boxplot del coef corr Spearman r"["s"] ~ "por grupo de meses"),
       x = "meses", 
       y = expression(r["s"])) 

## Scatterplot r_s vs altitud ----

# sttns.points: trae geometry
# SRTM_30.sttns: trae SRTM_30 (altitud)
# df_corr: trae los coef corr r_s
sttns.points = left_join(sttns.points, SRTM_30.sttns, by = 'CodigoEstacion') %>% 
  left_join(df_corr, by = 'CodigoEstacion')

# scatterplot v1: r_s vs altitud
sttns.points %>% filter(!is.na(SRTM_30 )) %>% 
ggplot(aes(x = d_R_Spearman, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión coef corr Spearman r"["s"] ~ "vs altitud "),
       x = expression(r["s"]),
       y = "Altitud") +
  theme_minimal()

## panel scatterplot: r_s por mes/grupo meses vs altitud ----
df_R_Spearman_mes = left_join(df_R_Spearman_mes, SRTM_30.sttns, by = 'CodigoEstacion')
df_R_Spearman_cat_mes = left_join(df_R_Spearman_cat_mes, SRTM_30.sttns, by = 'CodigoEstacion')

ggplot(df_R_Spearman_mes, aes(x = r_s, y = SRTM_30)) +
#ggplot(df_R_Spearman_cat_mes, aes(x = r_s, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión r"["s"] ~ "por mes vs altitud"),
  #labs(title = expression("Diagrama de dispersión r"["s"] ~ "por grupo de meses vs altitud"),
       x = expression(r["s"]),
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ mes, scales = "free")
  #facet_wrap(~ cat_mes, scales = "free")

## Distribución espacial de R_s ----
shp_depto <- st_read("Departamentos202208_shp/Depto.shp")
## regiones.shp <- https://github.com/jcmunozmora/WB_agroinsumos_Colombia/tree/main/01_mapas
##Regiones_Naturales.shp <- https://github.com/SMByC/puntosdecalor.ideam.gov.co/tree/master/page/data/shapes/regions
## https://www.researchgate.net/post/Mapa-fisiografico-o-regiones-naturales-de-Colombia
#shp_regions <- st_read("Regiones_shp/regiones.shp")
## Valencia, S., et al., (2023). 
## Spatio-temporal assessment of Gridded precipitation products across topographic and climatic gradients in Colombia.
shp_regions <- st_read("regiones_naturales_colombia/shp_regiones_naturales_colombia.shp")
#mundocol <- st_read("admin00/admin00_vecinos.shp")

## grafico mapa DPTO - regiones, puntos sttns con diferencias
ggplot() + 
  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) +
  geom_sf(data = shp_depto, alpha =  0.2)

## !!! Falta agregar elementos de los mapas: Norte, escala
ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = filter(sttns.points, !is.na(d_R_Spearman)), aes(fill = d_R_Spearman, color = d_R_Spearman), size = 1) +
  scale_fill_gradient(low = "red", high = "green", name = expression(r["s"])) +# Muestra los puntos con valores de r_s
  scale_color_gradient(low = "red", high = "green", name = expression(r["s"])) +
  #theme_minimal() +  # Tema minimalista
  labs(title = expression("Distribución espacial de r"[s]) , fill = "d_R_Spearman") +  # Título y leyenda
  # es mejor dejar la escala a un lado para facilitar lectura 
  theme(legend.position = "right") 

## Violin plot con boxplot de R_s por región ----
##  modificar el parametro a FALSE evita el error al hacer el join espacial
## https://stackoverflow.com/questions/68478179/how-to-resolve-spherical-geometry-failures-when-joining-spatial-data
#sf::sf_use_s2(FALSE)
sttns.points = st_join(sttns.points, shp_regions)

ggplot() + 
  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) + 
  #geom_sf(data = filter(sttns.points, Subregion=='Andes'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Caribe'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Pacifico'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Llanos'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Amazonas'), color = "blue", size = 3)
  geom_sf(data = filter(sttns.points, is.na(layer)), color = "blue", size = 3)

## asigna la Subregión manualmente a las sttns que NO cruzan en el join espacial
sttns.points = sttns.points %>% 
  mutate(layer = ifelse(is.na(layer) & DEPARTAMENTO %in% c('Arauca', 'Vichada'), 'Orinoquia', layer)) %>% 
  mutate(layer = ifelse(is.na(layer) & DEPARTAMENTO %in% c(ADZ, 'Sucre', 'Bolivar', 'Atlantico', 'Antioquia'), 'Caribe', layer)) %>% 
  mutate(layer = ifelse(is.na(layer) & DEPARTAMENTO == 'Nariño', 'Pacifico', layer)) %>% 
  mutate(layer = ifelse(is.na(layer) & DEPARTAMENTO == 'Amazonas', 'Amazonia', layer))

sttns.points = sttns.points %>% 
  mutate(region = factor(ifelse(layer == 'Amazonia', 'Amazonía',
                                ifelse(layer == 'Andes', 'Andina',
                                       ifelse(layer == 'Pacifico', 'Pacífica',
                                              ifelse(layer == 'Orinoquia', 'Orinoquía',
                                                     layer ))))))


table(sttns.points$region)

bx_ = bx_text(subset(sttns.points, !is.na(d_R_Spearman)), serie = region, r = d_R_Spearman)

ggplot(filter(sttns.points, !is.na(d_R_Spearman))) + 
  geom_violin(aes(x = region, y = d_R_Spearman, fill = region), trim = TRUE) +  # Usamos Región como eje X y r_s como eje Y
  geom_boxplot(aes(x = region, y = d_R_Spearman, fill = region), width = 0.2) +
  scale_fill_manual(values = c("Amazonía" = "#91B495", "Andina" = "#D4B698", "Caribe" = "#F8E595",
                               "Orinoquía" = "#B5E794", "Pacífica" = "#A59CB3")) +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(x = region, label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = region, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(x = "región", y = expression(r["s"]), title = expression("Violin plot con boxplot de r"["s"] ~ "por región natural")) +  # Etiquetas de los ejes y título
  # sin leyenda y tamaño del texto de los ejes
  theme_minimal() + 
  coord_cartesian(ylim = c(0, max(bx_$df_outliers$upper_bound)))+
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 
  
#theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar las etiquetas del eje X si es necesario

## panel de scatterplot: r_s por región vs altitud ----

ggplot(filter(sttns.points, !is.na(d_R_Spearman)),
       aes(x = d_R_Spearman, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión r"["s"] ~ "por región vs altitud"),
       x = expression(r["s"]),
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ region, scales = "free")

## 9.5 df_MAE:data frame para graficar MAE ------------------------------------
df_MAE = data.frame(CodigoEstacion = 
                      sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ])

for(i in 1:len_){
  # MAE: Error Medio Absoluto de la ts original -------------------------------
  df_MAE$MAE[i] = mae(sim = pol2_list[[i]]$`chirps-v2-0`, 
                      obs = pol2_list[[i]]$sttns,  na.rm=TRUE)
  
  # d_MAE: Error Medio Absoluto de la ts diff estacionalmente -----------------
  df_MAE$d_MAE[i] = mae(sim = pol2_list[[i]]$`d_chirps.v2.0`, 
                        obs = pol2_list[[i]]$d_sttns,  na.rm=TRUE)
}

## Boxplot comparación MAE ts diff estacionalmente vs ts original -------------
# df para graficar boxpot MAE ts original vs ts transformada
df_MAE_long = df_MAE %>% dplyr::select(-c(CodigoEstacion)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "MAE") %>% 
  mutate(serie = factor(serie, labels = c("diff estacionalmente", 
                                          "original")))

# crea objeto con texto de % outliers y texto de quartiles
bx_ = bx_text(df_MAE_long, serie, MAE)

## Al comparar graficamente el MAE de la ts original vs ts diff estacionalmente
## +obs1: se observa que la medida MAE para la ts original subestima el MAE para la ts diff estacionalmente
## +obs2: la ts diff estacionalmente tiene 1.2% menos de outliers 
## +obs3: el RIQ es mayor en la ts diff estacionalmente, la variabilidad del MAE aumenta (amplitud de la caja) 

ggplot(df_MAE_long, aes(x = fct_rev(serie), y = MAE, fill = serie)) +
  geom_boxplot() +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)-20), 
            size = 4, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del MAE", x = "serie", y = "MAE") +
  theme_minimal() +
  #coord_cartesian(ylim = c(0, max(bx_$df_outliers$upper_bound)))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## Boxplot MAE por mes --------------------------------------------------------

df_cat <- function(list = pol2_list, cat_, perf){
  # Argumentos: ----------------------------------------------------------------
  # 
  # list lista con dfs con col v_cat, d_chirps.v2.0 y d_sttns
  # cat_ str nombre de la col en c/elemento list a filtrar
  # perf str con medidas de desempeño
  #----------------------------------------------------------------------------
  len_ = length(list)
  # deriva el vector de categorías de la col arg v_cat 
  v_cat <- distinct(pol2_list[[1]], !!sym(cat_)) %>% as_vector() #unique(list[[i]]$cat_mes)
  
  df_ = data.frame(matrix(numeric(len_), nc = length(v_cat), nr = len_)) %>%
    `colnames<-`(v_cat)
  
  for(i in 1:len_){
    # Calcula la medida de desempeño para c/categoría
    values_ <- sapply(v_cat, function(m){
      # Filtra los datos para la categoría 'm'
      data_mes <- filter(list[[i]], !!sym(cat_) == m)
      # Calcula la medida de desempeño de acuerdo al arg perf
      if(perf == 'MAE'){
        perf_value <- mae(sim = data_mes$d_chirps.v2.0, 
                          obs = data_mes$d_sttns, na.rm = TRUE)
      }else if(perf == 'PBIAS'){
        perf_value <- pbias(sim = data_mes$d_chirps.v2.0, 
                            obs = data_mes$d_sttns, na.rm = TRUE, dec=2)
      }else{
        perf_value <- NSE(sim = data_mes$d_chirps.v2.0, 
                          obs = data_mes$d_sttns, na.rm = TRUE)
      }
      return(perf_value)
    })
    # guarda los resultados en el listado de R_Spearman
    df_[i,] <- values_
  }
  
  # crea la col CodigoEstacion
  df_$CodigoEstacion = sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ]
  
  df_ <- df_ %>%
    pivot_longer(cols = all_of(as.vector(v_cat)), 
                 names_to = cat_, 
                 values_to = perf)
  
  df_ = df_ %>% mutate(!!cat_ := factor(!!sym(cat_), levels = v_cat))
  return(df_)
}

df_MAE_mes <- df_cat(list = pol2_list, cat_ = "mes", perf = "MAE")

bx_ = bx_text(df_MAE_mes, serie = mes, r = MAE)

ggplot(df_MAE_mes, aes(x = mes, y = MAE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  theme_minimal() +
  labs(title = "Boxplot del MAE por mes",
       x = "Mes", 
       y = "MAE") #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Boxplot MAE por cat_mes (DJF, MAM, JJA, SON) ----

df_MAE_cat_mes <- df_cat(list = pol2_list, cat_ = "cat_mes", perf = "MAE")

bx_ = bx_text(df_MAE_cat_mes, cat_mes, MAE)

ggplot(df_MAE_cat_mes, aes(x = cat_mes, y = MAE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  theme_minimal() +
  labs(title = expression("Boxplot del MAE por grupo de meses"),
       x = "meses", 
       y = "MAE") 

## Scatterplot MAE vs altitud ----

# sttns.points: trae geometry
# df_MAE: trae los MAE
sttns.points = left_join(sttns.points, df_MAE, by = 'CodigoEstacion') 

# scatterplot v1: r_s vs altitud
sttns.points %>% filter(!is.na(SRTM_30 )) %>% 
  ggplot(aes(x = d_MAE, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión MAE vs altitud "),
       x = "MAE",
       y = "Altitud") +
  theme_minimal()

## panel scatterplot: MAE por mes/grupo meses vs altitud ----
df_MAE_mes = left_join(df_MAE_mes, SRTM_30.sttns, by = 'CodigoEstacion')
df_MAE_cat_mes = left_join(df_MAE_cat_mes, SRTM_30.sttns, by = 'CodigoEstacion')

#ggplot(df_MAE_mes, aes(x = MAE, y = SRTM_30)) +
ggplot(df_MAE_cat_mes, aes(x = MAE, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  #labs(title = expression("Diagrama de dispersión MAE por mes vs altitud"),
  labs(title = expression("Diagrama de dispersión MAE por grupo de meses vs altitud"),
       x = "MAE",
       y = "Altitud") +
  theme_minimal()+
 #facet_wrap(~ mes, scales = "free")
  facet_wrap(~ cat_mes, scales = "free")

## Distribución espacial de MAE ----

## grafico mapa DPTO - regiones, puntos sttns con diferencias
#ggplot() + 
#  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) +
#  geom_sf(data = shp_depto, alpha =  0.2)

## !!! Falta agregar elementos de los mapas: Norte, escala
ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = filter(sttns.points, !is.na(d_MAE)), aes(fill = d_MAE, color = d_MAE), size = 1) +
  scale_fill_gradient(low = "green", high = "red", name = "MAE") +# Muestra los puntos con valores de r_s
  scale_color_gradient(low = "green", high = "red", name = "MAE") +
  #theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de MAE" , fill = "d_MAE") +  # Título y leyenda
  # es mejor dejar la escala a un lado para facilitar lectura 
  theme(legend.position = "right") 

## Violin plot con boxplot de MAE por región ----
bx_ = bx_text(subset(sttns.points, !is.na(d_MAE)), serie = region, r = d_MAE)

ggplot(filter(sttns.points, !is.na(d_MAE))) + 
  geom_violin(aes(x = region, y = d_MAE, fill = region), trim = TRUE) +  # Usamos Región como eje X y r_s como eje Y
  geom_boxplot(aes(x = region, y = d_MAE, fill = region), width = 0.2) +
  scale_fill_manual(values = c("Amazonía" = "#91B495", "Andina" = "#D4B698", "Caribe" = "#F8E595",
                               "Orinoquía" = "#B5E794", "Pacífica" = "#A59CB3")) +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(x = region, label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = region, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(x = "región", y = "MAE", title = expression("Violin plot con boxplot de MAE por región natural")) +  # Etiquetas de los ejes y título
  # sin leyenda y tamaño del texto de los ejes
  theme_minimal() + 
  coord_cartesian(ylim = c(0, max(bx_$df_outliers$upper_bound)))+
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

#theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar las etiquetas del eje X si es necesario

## panel de scatterplot: MAE por región vs altitud ----

ggplot(filter(sttns.points, !is.na(d_MAE)),
       aes(x = d_MAE, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión MAE por región vs altitud"),
       x = "MAE",
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ region, scales = "free")

## 9.6 df_PBIAS: data frame para graficar PBIAS -------------------------------
df_PBIAS = data.frame(CodigoEstacion = 
                        sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ])

for(i in 1:len_){
  # PBIAS: porcentaje de sesgo de la ts original ------------------------------
  df_PBIAS$PBIAS[i] = pbias(sim = pol2_list[[i]]$`chirps-v2-0`,
                            obs = pol2_list[[i]]$sttns, na.rm=TRUE, dec=1)
  
  # d_PBIAS: porcentaje de sesgo de la ts diff estacionalmente ------------------
  df_PBIAS$d_PBIAS[i] = pbias(sim = pol2_list[[i]]$`d_chirps.v2.0`,
                              obs = pol2_list[[i]]$d_sttns, na.rm=TRUE, dec=1)
}

## Boxplot comparación PBIAS ts diff estacionalmente vs ts original -------------
# df para graficar boxpot MAE ts original vs ts transformada
df_PBIAS_long = df_PBIAS %>% dplyr::select(-c(CodigoEstacion)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "PBIAS") %>% 
  mutate(serie = factor(serie, labels = c("diff estacionalmente", "original")))

# crea objeto con texto de % outliers y texto de quartiles
bx_ = bx_text(df_PBIAS_long, serie, PBIAS)

## Al comparar graficamente el PBIAS de la ts original vs ts diff estacionalmente
## +obs1: se observa que la ts diff estacionalmente tiene un PBIAS negativo, indicando que d_CHIRPS subestima d_sttns
## +obs2: la ts diff estacionalmente tiene 11.5% más de outliers 
## +obs3: el RIQ es mayor en la ts diff estacionalmente, la variabilidad del PBIAS aumenta (amplitud de la caja) 

ggplot(df_PBIAS_long, aes(x = fct_rev(serie), y = PBIAS, fill = serie)) +
  geom_boxplot() +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(lower_bound)), 
            size = 4, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2), "%")), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2), "%")), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2), "%")), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del PBIAS", x = "serie", y = "PBIAS") +
  theme_minimal() +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), abs(min(bx_$df_outliers$lower_bound))))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## Boxplot PBIAS por mes --------------------------------------------------------

df_PBIAS_mes <- df_cat(list = pol2_list, cat_ = "mes", perf = "PBIAS")
# hay valores de PBIAS NAs porque la col d_chrpsv2.0 suma cero

bx_ = bx_text(df_PBIAS_mes %>% filter(!is.na(PBIAS)), serie = mes, r = PBIAS)

ggplot(df_PBIAS_mes, aes(x = mes, y = PBIAS)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(bx_$df_outliers$lower_bound)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q1, label = paste0("Q1: ", round(Q1, 2), "%")), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q2, label = paste0("Q2: ", round(Q2, 2), "%")), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q3, label = paste0("Q3: ", round(Q3, 2), "%")), 
            vjust = -1 , size = 4) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), abs(min(bx_$df_outliers$lower_bound))))+
  # título y nombres de ejes
  theme_minimal() +
  labs(title = "Boxplot del PBIAS por mes",
       x = "Mes", 
       y = "PBIAS") #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Boxplot PBIAS por cat_mes (DJF, MAM, JJA, SON) ----

df_PBIAS_cat_mes <- df_cat(list = pol2_list, cat_ = "cat_mes", perf = "PBIAS")

bx_ = bx_text(df_PBIAS_cat_mes, cat_mes, PBIAS)

ggplot(df_PBIAS_cat_mes, aes(x = cat_mes, y = PBIAS)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(lower_bound)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q1, label = paste0("Q1: ", round(Q1, 2), "%")), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q2, label = paste0("Q2: ", round(Q2, 2), "%")), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q3, label = paste0("Q3: ", round(Q3, 2), "%")), 
            vjust = -1 , size = 4) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), abs(min(bx_$df_outliers$lower_bound))))+
  # título y nombres de ejes
  theme_minimal() +
  labs(title = expression("Boxplot del PBIAS por grupo de meses"),
       x = "Meses", 
       y = "PBIAS") 

## Scatterplot PBIAS vs altitud ----

# sttns.points: trae geometry
# df_MAE: trae los MAE
sttns.points = left_join(sttns.points, df_PBIAS, by = 'CodigoEstacion') 

# scatterplot v1: r_s vs altitud
sttns.points %>% filter(!is.na(SRTM_30 ) & abs(d_PBIAS)<1000) %>% 
  ggplot(aes(x = d_PBIAS, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión PBIAS vs altitud "),
       x = "PBIAS",
       y = "Altitud") +
  theme_minimal()

## panel scatterplot: PBIAS por mes/grupo meses vs altitud ----
df_PBIAS_mes = left_join(df_PBIAS_mes, SRTM_30.sttns, by = 'CodigoEstacion')
df_PBIAS_cat_mes = left_join(df_PBIAS_cat_mes, SRTM_30.sttns, by = 'CodigoEstacion')

#ggplot(df_PBIAS_mes, aes(x = PBIAS, y = SRTM_30)) +
ggplot(df_PBIAS_cat_mes %>%  filter(abs(PBIAS)<1000), aes(x = PBIAS, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  #labs(title = expression("Diagrama de dispersión MAE por mes vs altitud"),
  labs(title = expression("Diagrama de dispersión PBIAS por grupo de meses vs altitud"),
       x = "PBIAS",
       y = "Altitud") +
  theme_minimal()+
  #facet_wrap(~ mes, scales = "free")
  facet_wrap(~ cat_mes, scales = "free")

## Distribución espacial de PBIAS ----

## grafico mapa DPTO - regiones, puntos sttns con diferencias
#ggplot() + 
#  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) +
#  geom_sf(data = shp_depto, alpha =  0.2)

## !!! Falta agregar elementos de los mapas: Norte, escala
ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = filter(sttns.points %>%  filter(abs(d_PBIAS)<300), !is.na(d_PBIAS)), aes(fill = d_PBIAS, color = d_PBIAS), size = 1) +
  #scale_fill_gradient(low = "green", high = "red", name = "PBIAS") +# Muestra los puntos con valores de r_s
  scale_fill_gradient2(
    low = "red", mid = "green", high = "red", midpoint = 0, name = 'PBIAS') +
  #scale_color_gradient(low = "green", high = "red", name = "PBIAS") +
  scale_color_gradient2(
    low = "red", mid = "green", high = "red", midpoint = 0, name = 'PBIAS') +
  #theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de PBIAS" , fill = "PBIAS") +  # Título y leyenda
  # es mejor dejar la escala a un lado para facilitar lectura 
  theme(legend.position = "right") 

## Violin plot con boxplot de PBIAS por región ----
bx_ = bx_text(subset(sttns.points, !is.na(d_PBIAS)), serie = region, r = d_PBIAS)

#ggplot(filter(sttns.points, !is.na(d_PBIAS))) +
ggplot(filter(sttns.points, abs(d_PBIAS)<500)) + 
  geom_violin(aes(x = region, y = d_PBIAS, fill = region), trim = TRUE) +  # Usamos Región como eje X y r_s como eje Y
  geom_boxplot(aes(x = region, y = d_PBIAS, fill = region), width = 0.2) +
  scale_fill_manual(values = c("Amazonía" = "#91B495", "Andina" = "#D4B698", "Caribe" = "#F8E595",
                               "Orinoquía" = "#B5E794", "Pacífica" = "#A59CB3")) +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(x = region, label = paste0(outlier_pct, "% \n outliers"), y = min(bx_$df_outliers$lower_bound)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = region, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(x = "Región", y = "PBIAS", title = expression("Violin plot con boxplot de PBIAS por región natural")) +  # Etiquetas de los ejes y título
  # sin leyenda y tamaño del texto de los ejes
  theme_minimal() + 
  coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), max(bx_$df_outliers$upper_bound)))+
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

#theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar las etiquetas del eje X si es necesario

## panel de scatterplot: PBIAS por región vs altitud ----

#ggplot(filter(sttns.points, !is.na(d_MAE)),
filter(sttns.points, abs(d_PBIAS)<500) %>% 
ggplot(aes(x = d_PBIAS, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión PBIAS por región vs altitud"),
       x = "PBIAS",
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ region, scales = "free")

## 9.7 df_NSE: data frame para graficar NSE -----------------------------------
df_NSE = data.frame(CodigoEstacion = 
                      sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ])

for(i in 1:len_){
  # NSE: Coef Eficiencia de Nash de la ts original ----------------------------
  df_NSE$NSE[i] = NSE(sim = pol2_list[[i]]$`chirps-v2-0`,
                      obs = pol2_list[[i]]$sttns,  na.rm=TRUE)
  
  # d_NSE: Coef Eficiencia de Nash de la ts diff estacionalmente --------------
  df_NSE$d_NSE[i] = NSE(sim = pol2_list[[i]]$`d_chirps.v2.0`,
                        obs = pol2_list[[i]]$d_sttns,  na.rm=TRUE)
}

## Boxplot comparación NSE ts diff estacionalmente vs ts original -------------
# df para graficar boxpot MAE ts original vs ts transformada
df_NSE_long = df_NSE %>% dplyr::select(-c(CodigoEstacion)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "NSE") %>% 
  mutate(serie = factor(serie, labels = c("diff estacionalmente", "original")))

# crea objeto con texto de % outliers y texto de quartiles
bx_ = bx_text(df_NSE_long, serie, NSE)

## Al comparar graficamente el NSE de la ts original vs ts diff estacionalmente
## +obs1: se observa que la medida NSE para ts original subestima el NSE para la ts diff estacionalmente
## +obs2: la ts diff estacionalmente tiene 3.4% menos outliers que la ts original
## +obs3: el RIQ es mayor en la ts original, la variabilidad del NSE aumenta (amplitud de la caja) 

ggplot(df_NSE_long, aes(x = fct_rev(serie), y = NSE, fill = serie)) +
  geom_boxplot() +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 4, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del NSE", x = "serie", y = "NSE") +
  theme_minimal() +
  #scale_y_continuous(labels = percent_format(scale = 1)) +
  #coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), abs(min(bx_$df_outliers$lower_bound))))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## Boxplot NSE por mes --------------------------------------------------------

df_NSE_mes <- df_cat(list = pol2_list, cat_ = "mes", perf = "NSE")
# hay valores de PBIAS NAs porque la col d_chrpsv2.0 suma cero

bx_ = bx_text(df_NSE_mes, serie = mes, r = NSE)

ggplot(df_NSE_mes, aes(x = mes, y = NSE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = -max(upper_bound)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = mes, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  coord_cartesian(ylim = c(-max(bx_$df_outliers$upper_bound), max(bx_$df_outliers$upper_bound)))+
  # título y nombres de ejes
  theme_minimal() +
  labs(title = "Boxplot del NSE por mes",
       x = "Mes", 
       y = "NSE") #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Boxplot NSE por cat_mes (DJF, MAM, JJA, SON) ----

df_NSE_cat_mes <- df_cat(list = pol2_list, cat_ = "cat_mes", perf = "NSE")

bx_ = bx_text(df_NSE_cat_mes, cat_mes, NSE)

ggplot(df_NSE_cat_mes, aes(x = cat_mes, y = NSE)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = cat_mes, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  theme_minimal() +
  labs(title = expression("Boxplot del NSE por grupo de meses"),
       x = "Meses", 
       y = "NSE") 

## Scatterplot NSE vs altitud ----

# sttns.points: trae geometry
# df_MAE: trae los MAE
sttns.points = left_join(sttns.points, df_NSE, by = 'CodigoEstacion') 

ggplot(sttns.points, aes(x = d_NSE, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión NSE vs altitud "),
       x = "NSE",
       y = "Altitud") +
  theme_minimal()

## panel scatterplot: NSE por mes/grupo meses vs altitud ----
df_NSE_mes = left_join(df_NSE_mes, SRTM_30.sttns, by = 'CodigoEstacion')
df_NSE_cat_mes = left_join(df_NSE_cat_mes, SRTM_30.sttns, by = 'CodigoEstacion')

#ggplot(df_NSE_mes, aes(x = NSE, y = SRTM_30)) +
ggplot(df_NSE_cat_mes, aes(x = NSE, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  #labs(title = expression("Diagrama de dispersión MAE por mes vs altitud"),
  labs(title = expression("Diagrama de dispersión NSE por grupo de meses vs altitud"),
       x = "NSE",
       y = "Altitud") +
  theme_minimal()+
  #facet_wrap(~ mes, scales = "free")
  facet_wrap(~ cat_mes, scales = "free")

## Distribución espacial de NSE ----

## grafico mapa DPTO - regiones, puntos sttns con diferencias
#ggplot() + 
#  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) +
#  geom_sf(data = shp_depto, alpha =  0.2)

## !!! Falta agregar elementos de los mapas: Norte, escala
ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = filter(sttns.points , !is.na(d_NSE)), aes(fill = d_NSE, color = d_NSE), size = 1) +
  scale_fill_gradient(low = "red", high = "green", name = "NSE") +# Muestra los puntos con valores de r_s
  scale_color_gradient(low = "red", high = "green", name = "NSE") +
  #theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de NSE" , fill = "NSE") +  # Título y leyenda
  # es mejor dejar la escala a un lado para facilitar lectura 
  theme(legend.position = "right") 

## Violin plot con boxplot de NSE por región ----
bx_ = bx_text(subset(sttns.points, !is.na(d_NSE)), serie = region, r = d_NSE)

ggplot(filter(sttns.points, !is.na(d_NSE))) +
  geom_violin(aes(x = region, y = d_NSE, fill = region), trim = TRUE) +  # Usamos Región como eje X y r_s como eje Y
  geom_boxplot(aes(x = region, y = d_NSE, fill = region), width = 0.2) +
  scale_fill_manual(values = c("Amazonía" = "#91B495", "Andina" = "#D4B698", "Caribe" = "#F8E595",
                               "Orinoquía" = "#B5E794", "Pacífica" = "#A59CB3")) +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(x = region, label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = region, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = region, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  labs(x = "Región", y = "NSE", title = expression("Violin plot con boxplot de NSE por región natural")) +  # Etiquetas de los ejes y título
  # sin leyenda y tamaño del texto de los ejes
  theme_minimal() + 
  #coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), max(bx_$df_outliers$upper_bound)))+
  theme(legend.position = "none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

#theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar las etiquetas del eje X si es necesario

## panel de scatterplot: NSE por región vs altitud ----

filter(sttns.points, !is.na(d_NSE)) %>% 
  ggplot(aes(x = d_NSE, y = SRTM_30)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = expression("Diagrama de dispersión NSE por región vs altitud"),
       x = "NSE",
       y = "Altitud") +
  theme_minimal()+
  facet_wrap(~ region, scales = "free")


# 10. qmap::fitQmap ----

## 10.1 fit Qmap and do Qmap ----

d_qm_QUANT = list()
d_qm_RQUANT = list()
qm_QUANT = list()
qm_RQUANT = list()
#d_qm_SSPLIN = list()
#qm_SSPLIN = list()
qm_PTF_power = list()
d_qm_PTF_power = list()
qm_PTF_linear = list()
d_qm_PTF_linear = list()
#qm_PTF_expasympt = list()
#d_qm_PTF_expasympt = list()
qm_PTF_scale = list()
d_qm_PTF_scale = list()
qm_PTF_power.x0 = list()
d_qm_PTF_power.x0 = list()
#qm_PTF_expasympt.x0 = list()
#d_qm_PTF_expasympt.x0 = list()

for(i in 1:len_){
  # Argumentos fitQmap: -------------------------------------------------------
  # 
  # obs: ts observado sttns
  # mod: ts modelo CHIRPS v2.0
  # method: str QUANT/RQUANT No paramétricos
  # qstep: num para construir malla espaciada quantile(mod,probs=seq(0,1,by=qstep))
  # wet.day: TRUE/FALSE se usa para igualar fracción días con pcp 
  #----------------------------------------------------------------------------
  
  # reeplaza en c/iter la estimacion de la ECDF
  qm_QUANT.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
                          method="QUANT", qstep=0.01, wet.day = TRUE)
  
  d_qm_QUANT.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
                            method="QUANT", qstep=0.01, wet.day = FALSE)
  
  # en c/iter estima Q-Q vía Local Linear Least Square Reg
  qm_RQUANT.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
                           method="RQUANT", qstep=0.01, wet.day = TRUE)
  
  d_qm_RQUANT.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
                             method="RQUANT", qstep=0.01, wet.day = FALSE)
  
  # OBSERVACIÓN 1 -------------------------------------------------------------
  # Al trasladar la ts con una cte para que sea positiva y wet.day = TRUE
  # no se obtiene mejora significativa en las medidas de desempeño
  # disminuye unas decimas pero se reporta el mismo comportamiento
  # en la ts desestacionalizada el MAE se incrementa y el NSE disminuye
  #----------------------------------------------------------------------------
  
  #cte[i] = abs(min(min(d_sttns[[i]]$sttns, na.rm = TRUE),
  #              min(d_chirps.v2.0[[i]]$`chirps-v2-0`)) + 1)
  
  #d_qm_QUANT.fit <- fitQmap(d_sttns[[i]]$sttns + cte[i],
  #                          d_chirps.v2.0[[i]]$`chirps-v2-0` + cte[i],
  #                          method="QUANT", qstep=0.01, wet.day = TRUE)
  
  #d_qm_RQUANT.fit <- fitQmap(d_sttns[[i]]$sttns + cte[i],
  #                           d_chirps.v2.0[[i]]$`chirps-v2-0` + cte[i],
  #                           method="RQUANT", qstep=0.01, wet.day = TRUE)
  
  # OBSERVACIÓN 2 -------------------------------------------------------------
  # con el método SSPLIN el MAE se incrementa y el NSE disminuye
  #----------------------------------------------------------------------------
  
  #qm_SSPLIN.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
  #                         method="SSPLIN", qstep=0.01, wet.day = TRUE)
  
  #d_qm_SSPLIN.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
  #                           method="SSPLIN", qstep=0.01, wet.day = FALSE)
  
  #----------------------------------------------------------------------------
  #  Quantile mapping using parametric transformations 
  #----------------------------------------------------------------------------
  
  qm_PTF_power.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
                              transfun="power", cost="MAE", qstep=0.01, wet.day = TRUE)
  
  d_qm_PTF_power.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
                                transfun="power", cost="MAE", qstep=0.01, wet.day = FALSE)
  
  qm_PTF_linear.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
                               transfun="linear", cost="MAE", qstep=0.01, wet.day = TRUE)
  
  d_qm_PTF_linear.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
                                 transfun="linear", cost="MAE", qstep=0.01, wet.day = FALSE)
  
  #qm_PTF_expasympt.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
  #                                transfun="expasympt", cost="MAE", qstep=0.01, wet.day = TRUE)
  
  #d_qm_PTF_expasympt.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
  #                                  transfun="expasympt", cost="MAE", qstep=0.01, wet.day = FALSE)
  
  qm_PTF_scale.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
                              transfun="scale", cost="MAE", qstep=0.01, wet.day = TRUE)
  
  d_qm_PTF_scale.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
                                transfun="scale", cost="MAE", qstep=0.01, wet.day = FALSE)
  
  qm_PTF_power.x0.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
                                 transfun="power.x0", cost="MAE", qstep=0.01, wet.day = TRUE)
  
  d_qm_PTF_power.x0.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
                                   transfun="power.x0", cost="MAE", qstep=0.01, wet.day = FALSE)
  
  #qm_PTF_expasympt.x0.fit <- fitQmap(pol2_list[[i]]$sttns, pol2_list[[i]]$`chirps-v2-0`,
  #                                   transfun="expasympt.x0", cost="MAE", qstep=0.01, wet.day = TRUE)
  
  #d_qm_PTF_expasympt.x0.fit <- fitQmap(d_sttns[[i]]$sttns, d_chirps.v2.0[[i]]$`chirps-v2-0`,
  #                                     transfun="expasympt.x0", cost="MAE", qstep=0.01, wet.day = FALSE)
  
  # Argumentos doQmap: -------------------------------------------------------
  # 
  # mod: ts modelo CHIRPS v2.0
  # fobj: salida de fitQmapQUANT/fitQmapRQUANT
  # type: opc interpolación para valores que no están en la malla quantile(mod,probs=seq(0,1,by=qstep))
  #----------------------------------------------------------------------------
  qm_QUANT[[i]] <- data.frame(QUANT_ = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
                                              qm_QUANT.fit, type="linear"))
                                              #qm_QUANT.fit, type="tricub"))  
  
  
  
  d_qm_QUANT[[i]] <- data.frame(QUANT_ = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
                                                d_qm_QUANT.fit, type="linear"))
                                                #d_qm_QUANT.fit, type="tricub"))  
  
  qm_RQUANT[[i]] <- data.frame(RQUANT_ = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
                                                qm_RQUANT.fit, type="linear"))
                                                #qm_RQUANT.fit, type="tricub"))  
  
  d_qm_RQUANT[[i]] <- data.frame(RQUANT_ = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
                                                  d_qm_RQUANT.fit, type="linear"))
                                                  #d_qm_RQUANT.fit, type="tricub"))  
  
  #d_qm_QUANT[[i]] <- data.frame(QUANT_ = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0` + cte[i], 
  #                                              d_qm_QUANT.fit, type="linear") - cte[i]) 
  
  #d_qm_RQUANT[[i]] <- data.frame(RQUANT_ = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0` + cte[i], 
  #                                                d_qm_RQUANT.fit, type="linear") - cte[i]) 
  
  #qm_SSPLIN[[i]] <- data.frame(SSPLIN = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
  #                                              qm_SSPLIN.fit))
  
  #d_qm_SSPLIN[[i]] <- data.frame(SSPLIN = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
  #                                                d_qm_SSPLIN.fit))
  
  #----------------------------------------------------------------------------
  #  Quantile mapping using parametric transformations 
  #----------------------------------------------------------------------------
  
  qm_PTF_power[[i]] <- data.frame(PTF_power = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
                                                     qm_PTF_power.fit))
  
  d_qm_PTF_power[[i]] <- data.frame(PTF_power = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
                                                       d_qm_PTF_power.fit))
  
  qm_PTF_linear[[i]] <- data.frame(PTF_linear = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
                                                       qm_PTF_linear.fit))
  
  d_qm_PTF_linear[[i]] <- data.frame(PTF_linear = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
                                                         d_qm_PTF_linear.fit))
  
  #qm_PTF_expasympt[[i]] <- data.frame(PTF_expasympt = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
  #                                                           qm_PTF_expasympt.fit))
  
  #d_qm_PTF_expasympt[[i]] <- data.frame(PTF_expasympt = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
  #                                                             d_qm_PTF_expasympt.fit))
  
  qm_PTF_scale[[i]] <- data.frame(PTF_scale = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
                                                     qm_PTF_scale.fit))
  
  d_qm_PTF_scale[[i]] <- data.frame(PTF_scale = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
                                                       d_qm_PTF_scale.fit))
  
  qm_PTF_power.x0[[i]] <- data.frame(PTF_power.x0 = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
                                                           qm_PTF_power.x0.fit))
  
  d_qm_PTF_power.x0[[i]] <- data.frame(PTF_power.x0 = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
                                                             d_qm_PTF_power.x0.fit))
  
  #qm_PTF_expasympt.x0[[i]] <- data.frame(PTF_expasympt.x0 = doQmap(pol2_list[[i]]$`chirps-v2-0`, 
  #                                                                 qm_PTF_expasympt.x0.fit))
  
  #d_qm_PTF_expasympt.x0[[i]] <- data.frame(PTF_expasympt.x0 = doQmap(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
  #                                                                   d_qm_PTF_expasympt.x0.fit))
  
  # OBSERVACIÓN 3 -------------------------------------------------------------
  # con NAs en la ts no ajusta fitQmap (Error)
  # ---------------------------------------------------------------------------
  
  #qm_QUANT.fit <- fitQmap(pol2_list[[i]]$d_sttns, pol2_list[[i]]$d_chirps.v2.0,
  #                        method="QUANT", qstep=0.01, wet.day = FALSE)
  
  #qm_QUANT[[i]] <- data.frame(QUANT_ = doQmap(pol2_list[[i]]$d_chirps.v2.0, 
  #                                            #qm_QUANT.fit, type="linear"))
  #                                            qm_QUANT.fit, type="tricub"))
  
}

## 10.2 crea la cols qm_{} y d_qm{} en pol2_list -------------------------------
for(i in 1:len_){
  ## crea cols qm_{(R)QUANT} y d_qm{(R)QUANT} en lista pol2_list
  pol2_list[[i]]$qm_QUANT = qm_QUANT[[i]]$QUANT_
  pol2_list[[i]]$d_qm_QUANT = c(rep(NA, round(df_pgram$period[i])),
                                d_qm_QUANT[[i]]$QUANT_)
  
  pol2_list[[i]]$qm_RQUANT = qm_RQUANT[[i]]$RQUANT_
  pol2_list[[i]]$d_qm_RQUANT = c(rep(NA, round(df_pgram$period[i])),
                                d_qm_RQUANT[[i]]$RQUANT_)
  
  #pol2_list[[i]]$qm_SSPLIN = qm_SSPLIN[[i]]$SSPLIN
  #pol2_list[[i]]$d_qm_SSPLIN = c(rep(NA, round(df_pgram$period[i])),
  #                               d_qm_SSPLIN[[i]]$SSPLIN)
  
  # crea cols qm_PTF_{function} y d_qm_PTF_{function} -------------------------
  pol2_list[[i]]$qm_PTF_power = qm_PTF_power[[i]]$PTF_power
  pol2_list[[i]]$d_qm_PTF_power = c(rep(NA, round(df_pgram$period[i])),
                                    d_qm_PTF_power[[i]]$PTF_power)
  
  pol2_list[[i]]$qm_PTF_linear = qm_PTF_linear[[i]]$PTF_linear
  pol2_list[[i]]$d_qm_PTF_linear = c(rep(NA, round(df_pgram$period[i])),
                                     d_qm_PTF_linear[[i]]$PTF_linear)
  
  #pol2_list[[i]]$qm_PTF_expasympt = qm_PTF_expasympt[[i]]$PTF_expasympt
  #pol2_list[[i]]$d_qm_PTF_expasympt = c(rep(NA, round(df_pgram$period[i])),
  #                                      d_qm_PTF_expasympt[[i]]$PTF_expasympt)
  
  pol2_list[[i]]$qm_PTF_scale = qm_PTF_scale[[i]]$PTF_scale
  pol2_list[[i]]$d_qm_PTF_scale = c(rep(NA, round(df_pgram$period[i])),
                                    d_qm_PTF_scale[[i]]$PTF_scale)
  
  pol2_list[[i]]$qm_PTF_power.x0 = qm_PTF_power.x0[[i]]$PTF_power.x0
  pol2_list[[i]]$d_qm_PTF_power.x0 = c(rep(NA, round(df_pgram$period[i])),
                                       d_qm_PTF_power.x0[[i]]$PTF_power.x0)
  
  #pol2_list[[i]]$qm_PTF_expasympt.x0 = qm_PTF_expasympt.x0[[i]]$PTF_expasympt.x0
  #pol2_list[[i]]$d_qm_PTF_expasympt.x0 = c(rep(NA, round(df_pgram$period[i])),
  #                                         d_qm_PTF_expasympt.x0[[i]]$PTF_expasympt.x0)
  
}

## Gráfico compara QM con chirps y sttns para REMOLINO [24040060] Santander ----
set.seed(02102025)
i = sample(1:len_, 1)

pol2_list[[i]]$Date = as.Date(pol2_list[[i]]$Date)

ggplot(pol2_list[[i]], aes(x = Date)) +
  geom_line(aes(y = qm_QUANT, color = "QUANT"), linewidth = 1) +
  geom_line(aes(y = `chirps-v2-0`, color = "chirps"), linewidth = 1) +
  geom_line(aes(y = sttns, color = "sttns"), linewidth = 1) +
  scale_color_manual(values = c("QUANT" = "red", "chirps" = "black", "sttns" = "blue")) +
  labs(color = "Serie", y = 'pcp', 
       title = 'TS original - estación REMOLINO [24040060] dpto. Santander') +
  theme_minimal() +
  theme(legend.position = c(0.95, 0.9), axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

ggplot(pol2_list[[i]] , aes(x = Date)) +
  geom_line(aes(y = d_qm_QUANT, color = "qm_QUANT"), linewidth = 1) +
  geom_line(aes(y = d_chirps.v2.0, color = "d_chirps"), linewidth = 1) +
  geom_line(aes(y = d_sttns, color = "d_sttns"), linewidth = 1, na.rm = TRUE) +
  scale_color_manual(values = c("qm_QUANT" = "red", "d_chirps" = "black", "d_sttns" = "blue"))+
  labs(color = "Serie", y = 'pcp desestacionalizada', 
       title = 'TS desestacionalizada - estación REMOLINO [24040060] dpto. Santander') +
  theme_minimal() + 
  theme(legend.position = c(0.95, 0.9), axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## 10.3 df_corr: agrega cols qm_{}_R_Speaman y d_qm_{}_R_Speaman --------------
for(i in 1:len_){
  
  # QUANT R_Spearman ----------------------------------------------------------
  df_corr$qm_QUANT_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_QUANT,
                                             obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_corr$d_qm_QUANT_R_Spearman[i] = rSpearman(sim = d_qm_QUANT[[i]]$QUANT_,
                                               obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  # RQUANT R_Spearman ---------------------------------------------------------
  df_corr$qm_RQUANT_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_RQUANT,
                                              obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_corr$d_qm_RQUANT_R_Spearman[i] = rSpearman(sim = d_qm_RQUANT[[i]]$RQUANT_,
                                                obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  # SSPLIN R_Spearman ---------------------------------------------------------
  #df_corr$qm_SSPLIN_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_SSPLIN, 
  #                                              obs = pol2_list[[i]]$sttns, na.rm = TRUE)
  
  #df_corr$d_qm_SSPLIN_R_Spearman[i] = rSpearman(sim = d_qm_SSPLIN[[i]]$SSPLIN, 
  #                                              obs = d_sttns[[i]]$sttns, na.rm = TRUE)
  
  # PTF R_Spearman ------------------------------------------------------------
  
  df_corr$d_qm_PTF_power_R_Spearman[i] = rSpearman(sim = d_qm_PTF_power[[i]]$PTF_power,
                                                   obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_corr$qm_PTF_power_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_PTF_power,
                                                 obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_corr$d_qm_PTF_linear_R_Spearman[i] = rSpearman(sim = d_qm_PTF_linear[[i]]$PTF_linear,
                                                    obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_corr$qm_PTF_linear_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_PTF_linear,
                                                  obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  #df_corr$d_qm_PTF_expasympt_R_Spearman[i] = rSpearman(sim = d_qm_PTF_expasympt[[i]]$PTF_expasympt,
  #                                                     obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  #df_corr$qm_PTF_expasympt_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_PTF_expasympt,
  #                                                   obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_corr$d_qm_PTF_scale_R_Spearman[i] = rSpearman(sim = d_qm_PTF_scale[[i]]$PTF_scale,
                                                   obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_corr$qm_PTF_scale_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_PTF_scale,
                                                 obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_corr$d_qm_PTF_power.x0_R_Spearman[i] = rSpearman(sim = d_qm_PTF_power.x0[[i]]$PTF_power.x0,
                                                      obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_corr$qm_PTF_power.x0_R_Spearman[i] = rSpearman(sim = pol2_list[[i]]$qm_PTF_power.x0,
                                                    obs = pol2_list[[i]]$sttns, na.rm=TRUE)
}

## Boxplot R_Spearman de métodos BC en ts original vs diff estacionalmente ------
df_Spearman = df_corr %>% dplyr::select(-c(CodigoEstacion, R_pearson, d_R_pearson)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "r_s") %>% 
  mutate(serie = factor(gsub('pearman','',serie)))

quantile(df_corr$d_R_Spearman)
quantile(df_corr$d_qm_QUANT_R_Spearman)

bx_ = bx_text(df_Spearman, serie, r_s)

ggplot(df_Spearman, aes(x = fct_rev(serie), y = r_s, fill = serie)) +
  geom_boxplot() +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del coef de Spearman", x = "serie", y = expression(r[s])) +
  theme_minimal() +
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

# 10.4 df_MAE: agrega cols qm_{}_MAE y d_qm_{}_MAE ----------------
for(i in 1:len_){
  # QUANT MAE -----------------------------------------------------------------
  df_MAE$d_qm_QUANT_MAE[i] = mae(sim = pol2_list[[i]]$d_qm_QUANT, 
                                 obs = pol2_list[[i]]$d_sttns,  na.rm=TRUE)
  
  df_MAE$qm_QUANT_MAE[i] = mae(sim = pol2_list[[i]]$qm_QUANT, 
                               obs = pol2_list[[i]]$sttns,  na.rm=TRUE)
  
  # RQUANT MAE ----------------------------------------------------------------
  df_MAE$d_qm_RQUANT_MAE[i] = mae(sim = pol2_list[[i]]$d_qm_RQUANT, 
                                  obs = pol2_list[[i]]$d_sttns,  na.rm=TRUE)
  
  df_MAE$qm_RQUANT_MAE[i] = mae(sim = pol2_list[[i]]$qm_RQUANT, 
                                obs = pol2_list[[i]]$sttns,  na.rm=TRUE)
  
  # SSPLIN MAE ----------------------------------------------------------------
  #df_MAE$d_qm_SSPLIN_MAE[i] = mae(sim = pol2_list[[i]]$SSPLIN, 
  #                                obs = pol2_list[[i]]$d_sttns, na.rm = TRUE)
  #df_MAE$qm_SSPLIN_MAE[i] = mae(sim = pol2_list[[i]]$SSPLIN, 
  #                              obs = pol2_list[[i]]$sttns, na.rm = TRUE)
  
  # PTF MAE -------------------------------------------------------------------
  
  df_MAE$d_qm_PTF_power_MAE[i] = mae(sim = d_qm_PTF_power[[i]]$PTF_power,
                                     obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_MAE$qm_PTF_power_MAE[i] = mae(sim = pol2_list[[i]]$qm_PTF_power,
                                   obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_MAE$d_qm_PTF_linear_MAE[i] = mae(sim = d_qm_PTF_linear[[i]]$PTF_linear,
                                      obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_MAE$qm_PTF_linear_MAE[i] = mae(sim = pol2_list[[i]]$qm_PTF_linear,
                                    obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  #df_MAE$d_qm_PTF_expasympt_MAE[i] = mae(sim = d_qm_PTF_expasympt[[i]]$PTF_expasympt,
  #                                       obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  #df_MAE$qm_PTF_expasympt_MAE[i] = mae(sim = pol2_list[[i]]$qm_PTF_expasympt,
  #                                     obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_MAE$d_qm_PTF_scale_MAE[i] = mae(sim = d_qm_PTF_scale[[i]]$PTF_scale,
                                     obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_MAE$qm_PTF_scale_MAE[i] = mae(sim = pol2_list[[i]]$qm_PTF_scale,
                                   obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_MAE$d_qm_PTF_power.x0_MAE[i] = mae(sim = d_qm_PTF_power.x0[[i]]$PTF_power.x0,
                                        obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_MAE$qm_PTF_power.x0_MAE[i] = mae(sim = pol2_list[[i]]$qm_PTF_power.x0,
                                      obs = pol2_list[[i]]$sttns, na.rm=TRUE)
}

## Boxplot MAE de métodos BC en ts original vs diff estacionalmente -----------
df_MAE_long = df_MAE %>% dplyr::select(-c(CodigoEstacion)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "MAE") %>% 
  mutate(serie = factor(serie))

df_MAE = df_MAE %>% mutate(
  across(
    # crea cols pct_cambio con cols que empiezan por {qm_} y terminan por {_MAE}
    .cols = matches("^qm_.*_MAE$"),
    .fns = ~ (.x - MAE) / MAE,
    .names = "pct_cambio_{col}"),
  across(
    # crea cols pct_cambio con cols que empiezan por {d_qm_} y terminan por {_MAE}
    .cols = matches("^d_qm_.*_MAE$"),
    .fns = ~ (.x - d_MAE) / d_MAE,
    .names = "pct_cambio_{col}"))

f_pct_cambio <- function(df, col1, col2, perf){
  # Argumentos: ----------------------------------------------------------------
  # 
  # df data frame con cols col1 y col2
  # col1 num col medida desempeño de referencia  
  # col2 num col medida desempeño de Bias Correction
  # perf str medida de desempeño PBIAS/MAE/NSE
  #----------------------------------------------------------------------------
  df %>%
    { 
      if (perf == "PBIAS") {
        # si la medida de desempeño es PBIAS filtra los valores absolutos menores a 300 (PBIAS_P5 y PBIAS_P95)
        .x <- filter(., abs({{ col1 }}) < 300)
        summarise(.x, pct = (sum(abs({{ col2 }})) - sum(abs({{ col1 }}))) / sum(abs({{ col1 }})))
        } else {
          summarise(., pct = (sum({{ col2 }}) - sum({{ col1 }})) / sum({{ col1 }}))
          }
      } %>%
    # devuelve el valor en lugar de una columna
    pull(pct)
  }

# crea objeto con texto de % outliers y texto de quartiles
bx_ = bx_text(df_MAE_long, serie, MAE)

# crea objeto con texto % de cambio para c/método QM
nam_ = names(df_MAE)[str_detect(names(df_MAE), '^qm|^d_qm')]
bx_pct_cambio = data.frame(serie = nam_,  pct_ = numeric(length(nam_)))

# asigna la salida de f_pct_cambio a la columna pct_ para la fila de c/serie == i
for(i in nam_){
  if(str_detect(i, '^d_qm')){
    bx_pct_cambio$pct_[bx_pct_cambio$serie == i] = f_pct_cambio(df_MAE, d_MAE, !!sym(i), perf = "MAE")
  }else{
    bx_pct_cambio$pct_[bx_pct_cambio$serie == i] = f_pct_cambio(df_MAE, MAE, !!sym(i), perf = "MAE")
  }
}

# join df_outliers con bx_pct_cambio para traer cols max_r y upper_bound
bx_pct_cambio = inner_join(dplyr::select(bx_$df_outliers, c(serie, max_r, upper_bound)),
                           bx_pct_cambio, by = 'serie') %>% 
  mutate(pct_ = round(pct_*100,2))

ggplot(df_MAE_long, aes(x = fct_rev(serie), y = MAE, fill = serie)) +
  geom_boxplot() +
  # texto % cambio en c/método QM
  geom_text(data = bx_pct_cambio, aes(label = paste0(pct_, "% \n cambio"), y = max(max_r)), 
            size = 4, color = "black", fontface = "bold") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)-20), 
            size = 4, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del MAE", x = "serie", y = "MAE") +
  theme_minimal() +
  #coord_cartesian(ylim = c(0, max(bx_$df_outliers$upper_bound)))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

# 10.5 df_PBIAS: agrega cols qm_{}_PBIAS y d_qm_{}_PBIAS ----------------------
for(i in 1:len_){
  # QUANT PBIAS ---------------------------------------------------------------
  df_PBIAS$d_qm_QUANT_PBIAS[i] = pbias(sim = pol2_list[[i]]$d_qm_QUANT,
                                       obs = pol2_list[[i]]$d_sttns, na.rm=TRUE, dec=1)
  
  df_PBIAS$qm_QUANT_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_QUANT,
                                     obs = pol2_list[[i]]$sttns, na.rm=TRUE, dec=1)
  
  # RQUANT PBIAS --------------------------------------------------------------
  df_PBIAS$d_qm_RQUANT_PBIAS[i] = pbias(sim = pol2_list[[i]]$d_qm_RQUANT,
                                        obs = pol2_list[[i]]$d_sttns, na.rm=TRUE, dec=1)
  
  df_PBIAS$qm_RQUANT_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_RQUANT,
                                      obs = pol2_list[[i]]$sttns, na.rm=TRUE, dec=1)
  
  #df_PBIAS$d_qm_SSPLIN_PBIAS[i] = pbias(sim = pol2_list[[i]]$d_qm_SSPLIN, 
  #                                      obs = pol2_list[[i]]$d_sttns, na.rm=TRUE, dec=1)
  
  #df_PBIAS$qm_SSPLIN_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_SSPLIN, 
  #                                    obs = pol2_list[[i]]$sttns, na.rm=TRUE, dec=1)
  
  # PTF PBIAS -----------------------------------------------------------------
  
  df_PBIAS$d_qm_PTF_power_PBIAS[i] = pbias(sim = d_qm_PTF_power[[i]]$PTF_power,
                                           obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_PBIAS$qm_PTF_power_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_PTF_power,
                                         obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_PBIAS$d_qm_PTF_linear_PBIAS[i] = pbias(sim = d_qm_PTF_linear[[i]]$PTF_linear,
                                            obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_PBIAS$qm_PTF_linear_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_PTF_linear,
                                          obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  #df_PBIAS$d_qm_PTF_expasympt_PBIAS[i] = pbias(sim = d_qm_PTF_expasympt[[i]]$PTF_expasympt,
  #                                             obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  #df_PBIAS$qm_PTF_expasympt_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_PTF_expasympt,
  #                                           obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_PBIAS$d_qm_PTF_scale_PBIAS[i] = pbias(sim = d_qm_PTF_scale[[i]]$PTF_scale,
                                           obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_PBIAS$qm_PTF_scale_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_PTF_scale,
                                         obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_PBIAS$d_qm_PTF_power.x0_PBIAS[i] = pbias(sim = d_qm_PTF_power.x0[[i]]$PTF_power.x0,
                                              obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_PBIAS$qm_PTF_power.x0_PBIAS[i] = pbias(sim = pol2_list[[i]]$qm_PTF_power.x0,
                                            pol2_list[[i]]$sttns, na.rm=TRUE)
}

## Boxplot PBIAS de métodos BC en ts original vs diff estacionalmente ---------

df_PBIAS_long = df_PBIAS %>% dplyr::select(-c(CodigoEstacion)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "PBIAS") %>% 
  mutate(serie = factor(serie))

bx_ = bx_text(df_PBIAS_long, serie, PBIAS)
nam_ = names(df_PBIAS)[str_detect(names(df_PBIAS), '^qm|^d_qm')]
bx_pct_cambio = data.frame(serie = nam_,  pct_ = numeric(length(nam_)))

for(i in nam_){
  if(str_detect(i, '^d_qm')){
    bx_pct_cambio$pct_[bx_pct_cambio$serie == i] = f_pct_cambio(df_PBIAS, d_PBIAS, !!sym(i), perf = "PBIAS")
  }else{
    bx_pct_cambio$pct_[bx_pct_cambio$serie == i] = f_pct_cambio(df_PBIAS, PBIAS, !!sym(i), perf = "PBIAS")
  }
}

bx_pct_cambio = inner_join(dplyr::select(bx_$df_outliers, c(serie, max_r, upper_bound)),
                           bx_pct_cambio, by = 'serie') %>% 
  mutate(pct_ = round(pct_*100,2))

ggplot(df_PBIAS_long, aes(x = fct_rev(serie), y = PBIAS, fill = serie)) +
  geom_boxplot() +
  # texto % cambio en c/método qm
  geom_text(data = bx_pct_cambio, aes(label = paste0(pct_, "% \n cambio"), y = max(upper_bound)), 
            size = 4, color = "black", fontface = "bold") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(lower_bound)-10), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), max(bx_$df_outliers$upper_bound)))+
  # título y nombres de ejes
  labs(title = "Boxplots del PBIAS", x = "serie", y = "PBIAS") +
  theme_minimal() +
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

# 10.6 df_NSE: agrega cols qm_{}_NSE y d_qm_{}_NSE ----------------------------
for(i in 1:len_){
  # QUANT NSE -----------------------------------------------------------------
  df_NSE$d_qm_QUANT_NSE[i] = NSE(sim = pol2_list[[i]]$d_qm_QUANT,
                                 obs = pol2_list[[i]]$d_sttns, na.rm=TRUE)
  
  df_NSE$qm_QUANT_NSE[i] = NSE(sim = pol2_list[[i]]$qm_QUANT,
                               obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  # RQUANT NSE ----------------------------------------------------------------
  df_NSE$d_qm_RQUANT_NSE[i] = NSE(sim = pol2_list[[i]]$d_qm_RQUANT,
                                  obs = pol2_list[[i]]$d_sttns,  na.rm=TRUE)
  
  df_NSE$qm_RQUANT_NSE[i] = NSE(sim = pol2_list[[i]]$qm_RQUANT,
                                obs = pol2_list[[i]]$sttns,  na.rm=TRUE)
  
  #df_NSE$d_qm_SSPLIN_NSE[i] = NSE(sim = pol2_list[[i]]$d_qm_SSPLIN, 
  #                                obs = pol2_list[[i]]$d_sttns, na.rm=TRUE)
  
  #df_NSE$qm_SSPLIN_NSE[i] = NSE(sim = pol2_list[[i]]$qm_SSPLIN, 
  #                              obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  # PFT NSE -------------------------------------------------------------------
  
  df_NSE$d_qm_PTF_power_NSE[i] = NSE(sim = d_qm_PTF_power[[i]]$PTF_power,
                                     obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_NSE$qm_PTF_power_NSE[i] = NSE(sim = pol2_list[[i]]$qm_PTF_power,
                                   obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_NSE$d_qm_PTF_linear_NSE[i] = NSE(sim = d_qm_PTF_linear[[i]]$PTF_linear,
                                      obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_NSE$qm_PTF_linear_NSE[i] = NSE(sim = pol2_list[[i]]$qm_PTF_linear,
                                    obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  #df_NSE$d_qm_PTF_expasympt_NSE[i] = NSE(sim = d_qm_PTF_expasympt[[i]]$PTF_expasympt,
  #                                       obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  #df_NSE$qm_PTF_expasympt_NSE[i] = NSE(sim = pol2_list[[i]]$qm_PTF_expasympt,
  #                                     obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_NSE$d_qm_PTF_scale_NSE[i] = NSE(sim = d_qm_PTF_scale[[i]]$PTF_scale,
                                     obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_NSE$qm_PTF_scale_NSE[i] = NSE(sim = pol2_list[[i]]$qm_PTF_scale,
                                   obs = pol2_list[[i]]$sttns, na.rm=TRUE)
  
  df_NSE$d_qm_PTF_power.x0_NSE[i] = NSE(sim = d_qm_PTF_power.x0[[i]]$PTF_power.x0,
                                        obs = d_sttns[[i]]$sttns, na.rm=TRUE)
  
  df_NSE$qm_PTF_power.x0_NSE[i] = NSE(sim = pol2_list[[i]]$qm_PTF_power.x0,
                                      pol2_list[[i]]$sttns, na.rm=TRUE)
}

## Boxplot NSE de métodos BC en ts original vs diff estacionalmente -----------
df_NSE_long = df_NSE %>% dplyr::select(-c(CodigoEstacion)) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "NSE") %>% 
  mutate(serie = factor(serie))

bx_ = bx_text(df_NSE_long, serie, NSE)
nam_ = names(df_NSE)[str_detect(names(df_NSE), '^qm|^d_qm')]
bx_pct_cambio = data.frame(serie = nam_,  pct_ = numeric(length(nam_)))

for(i in nam_){
  if(str_detect(i, '^d_qm')){
    bx_pct_cambio$pct_[bx_pct_cambio$serie == i] = f_pct_cambio(df_NSE, d_NSE, !!sym(i), perf = "NSE")
  }else{
    bx_pct_cambio$pct_[bx_pct_cambio$serie == i] = f_pct_cambio(df_NSE, NSE, !!sym(i), perf = "NSE")
  }
}

bx_pct_cambio = inner_join(dplyr::select(bx_$df_outliers, c(serie, max_r, upper_bound)),
                           bx_pct_cambio, by = 'serie') %>% 
  mutate(pct_ = round(pct_*100,2))

ggplot(df_NSE_long, aes(x = fct_rev(serie), y = NSE, fill = serie)) +
  geom_boxplot() +
  # texto % cambio en c/método qm
  geom_text(data = bx_pct_cambio, aes(label = paste0(pct_, "% \n cambio"), y = max(upper_bound)), 
            size = 4, color = "black", fontface = "bold") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del NSE", x = "serie", y = "PBIAS") +
  theme_minimal() +
  #coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), 1))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

#quantile(NSE_)
#mean(NSE_) # 0.33
# indica que, Chirps tiene una buena predicción, mejor que usar 
## la media de los valores observados como predicción, sin embargo
### el promedio está lejos de 1, aún hay un margen significativo para mejorar
#hist(NSE_)

## 10.7 Ensamble Bias Correction ----------------------------------------------

# df_error: data frame join de los dfs medidas_desempeño con los métodos QM ----
df_error = left_join(df_MAE, df_PBIAS, by = c('CodigoEstacion')) %>% 
  left_join(df_NSE, by = c('CodigoEstacion'))

df_error <- df_error %>%
  # crea las cols pct_cambio_{}
  mutate(
    # aplica la fx abs a las cols que empiezan por {d} y terminan por {_PBIAS}
    across(
      .cols = matches("^d.*_PBIAS$"), .fns = ~ abs(.x), .names = "{col}_abs"),
    # aplica la función % cambio a las cols que empiezan {d_qm} y terminan por {_MAE}
    across(
      .cols = matches("^d_qm_.*_MAE$"), .fns = ~ (.x - d_MAE) / d_MAE, 
      .names = "pct_cambio_{col}"),
    # aplica la función % cambio a las cols que empiezan {d_qm} y terminan por {_PBIAS_abs}
    across(
      .cols = matches("^d_qm_.*_PBIAS_abs$"), .fns = ~ (.x - d_PBIAS_abs) / d_PBIAS_abs,
      .names = "pct_cambio_{col}"),
    # aplica la función % cambio a las cols que empiezan {d_qm} y terminan por {_NSE}
    across(
      .cols = matches("^d_qm_.*_NSE$"), .fns = ~ (.x - d_NSE) / d_NSE,
      .names = "pct_cambio_{col}")
  )

# df_error_wide: pivota df_error para obtener cols serie (método QM) y serie2 ----
df_error_wide = df_error %>% 
  # pivota cols pct_cambio_ a filas
  pivot_longer(cols = starts_with('pct_cambio_'), names_to = "serie", values_to = "pct_cambio") %>% 
  # crea la cols serie2 (medida_desempeño) y serie (método QM)
  mutate(serie2 = ifelse(str_detect(serie, '_MAE'), 'MAE',
                         ifelse(str_detect(serie, '_PBIAS'), 'PBIAS', 'NSE')),
         serie = gsub('pct_cambio_|_MAE|_PBIAS_abs|_NSE', '', serie)) %>% 
  # selecciona las cols ID, serie, serie2 y pct_cambio
  dplyr::select(c('CodigoEstacion', 'serie', 'serie2', 'pct_cambio')) %>% 
  # pivota filas de la col serie2 (medida_desempeño) a cols
  pivot_wider(id_cols = c('CodigoEstacion', 'serie'), names_from = 'serie2', values_from = 'pct_cambio'  )

# df_optim: data frame con malla de pesos para las medidas de desempeño -------
df_optim <- expand.grid(
  a1 = seq(0, 1, by = 0.01),
  a2 = seq(0, 1, by = 0.01)) %>%
  mutate(a3 = 1 - a1 - a2) %>%
  filter(a3 >= 0 & a3 <= 1)

for(i in 1:nrow(df_optim)){
  
  rank_error = df_error_wide %>% 
    # crea la col pct_cambio_ como combinación lineal de sum_{j = 1}^{3} w_j (medida_desempeño)
    mutate(pct_cambio_ = (df_optim$a1[i]*MAE) + (df_optim$a2[i]*PBIAS) - (df_optim$a3[i]*NSE) ) %>% 
    # agrupa por ID
    group_by(CodigoEstacion) %>% 
    # filtra las filas con el min pct_cambio_
    slice_min(order_by = pct_cambio_, n = 1, with_ties = FALSE)
  
  # join con rank_error para traer serie (método QM)
  df_error_ensemble = left_join(df_error, dplyr::select(rank_error, c('CodigoEstacion', 'serie')),
                                by = 'CodigoEstacion') %>% 
    # aplica función por filas
    rowwise() %>%
    # crea cols d_ensemble_{medida_desempeño} a partir del valor de la col {{serie}}_{medida_desempeño}
    mutate(d_ensemble_MAE = get(paste0(serie, "_MAE")),
           d_ensemble_PBIAS = get(paste0(serie, "_PBIAS")),
           d_ensemble_NSE = get(paste0(serie, "_NSE"))) %>%
    ungroup()
  
  # crea las cols pct_cambio_{medida_desmepño} de toda la col
  df_optim$pct_cambio_MAE[i] = f_pct_cambio(df_error_ensemble, 
                                            d_MAE, d_ensemble_MAE , perf = "MAE")
  df_optim$pct_cambio_PBIAS[i] = f_pct_cambio(df_error_ensemble, d_PBIAS, 
                                              d_ensemble_PBIAS, perf = "PBIAS")
  df_optim$pct_cambio_NSE[i] = f_pct_cambio(df_error_ensemble, d_NSE, 
                                            d_ensemble_NSE, perf = "NSE")
}

#df_optim$pct_cambio = df_optim$pct_cambio_MAE + (df_optim$pct_cambio_PBIAS) - df_optim$pct_cambio_NSE
## crea la col pct_cambio como la combinación lineal de las cols pct_cambio_{medida_desempeño}
## ideal valores negativos en pct_cambio_MAE indica disminución MAE (más cercano a cero)
## ideal valores negativos en pct_cambio_PBIAS indica disminución PBIAS_{abs} (más cercano a cero)
## ideal valores postivos en pct_cambio_NSE indica incremento NSE (más cercano a uno)
## luego para el NSE en la comb. lineal se usa el negativo
## entonces se busca el min de pct_cambio, indica que mejora el desempeño (menor error y pbias, mayor NSE) 

df_optim = df_optim %>% mutate(
  pct_cambio = pct_cambio_MAE + pct_cambio_PBIAS - pct_cambio_NSE)

## el max pct_cambio_PBIAS coincide con el min pct_cambio_MAE  
## se gana 0.007 en MAE y 0.015 en NSE mientras que en PBIAS se pierde 0.002
w_optim = slice_min(df_optim, pct_cambio, n = 50, with_ties = FALSE) %>%
  filter(min(pct_cambio_MAE) == pct_cambio_MAE &  
           max(pct_cambio_PBIAS) == pct_cambio_PBIAS &
           max(pct_cambio_NSE) == pct_cambio_NSE)

# pasa de relación 22% (pct_cambio MAE +NSE)/(pct_cambio PBIAS) a relación 13% (pct_cambio MAE +NSE)/(pct_cambio PBIAS)
# disminuye 10% ganando 4.3% en (pct_cambio MAE + NSE)
pct_optim = (w_optim$pct_cambio_MAE - w_optim$pct_cambio_NSE)/abs(w_optim$pct_cambio_PBIAS)

p_corte = seq(0, 1200, 50)
for(i in p_corte){
  df_min_pct_cambio_MAE = slice_min(df_optim, pct_cambio_MAE, n = i, with_ties = FALSE)
  df_max_pct_cambio_NSE = slice_max(df_optim, pct_cambio_NSE, n = i, with_ties = FALSE)
  
  df_min_pct_cambio = inner_join(df_min_pct_cambio_MAE, 
                                 dplyr::select(df_max_pct_cambio_NSE, c(a1, a2, a3)),
                                 by = c('a1', 'a2', 'a3')) %>% 
    slice_min(pct_cambio, n = 1, with_ties = FALSE)
  
  MAE_cambio = df_min_pct_cambio$pct_cambio_MAE - w_optim$pct_cambio_MAE
  NSE_cambio = df_min_pct_cambio$pct_cambio_NSE - w_optim$pct_cambio_NSE
  PBIAS_cambio = df_min_pct_cambio$pct_cambio_PBIAS - w_optim$pct_cambio_PBIAS
  
  pct_ = (MAE_cambio + NSE_cambio)/abs(PBIAS_cambio)
  pct_disminuye = pct_optim - pct_
  pct_aumenta = MAE_cambio + NSE_cambio
  
  #print(paste0('top ',i, ': disminuye ', round(pct_disminuye, 3), 
  #             ' ganando ', round(pct_aumenta,3),
  #             ' con relación ', round(pct_aumenta/pct_disminuye, 3) )) 
  #print(paste0('top ',i, ': pct de mejora respecto a PBIAS:', round(pct_,3))) 
  #print(paste0('top ',i, ': pct_cambio_PBIAS:', 
  #             round(df_min_pct_cambio$pct_cambio_PBIAS,2)))
}

# top 600, para seleccionar hasta que punto disminuir pct_cambio_PBIAS para ganar más pct_cambio_MAE y NSE
df_min_pct_cambio_MAE = slice_min(df_optim, pct_cambio_MAE, n = 600, with_ties = FALSE)
df_max_pct_cambio_NSE = slice_max(df_optim, pct_cambio_NSE, n = 600, with_ties = FALSE)

w_min_pct_cambio = inner_join(df_min_pct_cambio_MAE, 
                              dplyr::select(df_max_pct_cambio_NSE, c(a1, a2, a3)),
                              by = c('a1', 'a2', 'a3')) %>% 
  slice_min(pct_cambio, n = 1, with_ties = FALSE)

# df_error_ensemble: data frame para seleccionar cols ------------------------
rank_error = df_error_wide %>% 
  # w_optim a1: 0.8, a2: 0.19 y a3: 0.01
  mutate(pct_cambio_ = (w_optim$a1*MAE) + (w_optim$a2*PBIAS) - (w_optim$a3*NSE) ) %>%
  # w_min_pct_cambio a1: 0.93, a2: 0.06 y a3: 0.01
  #mutate(pct_cambio_ = (w_min_pct_cambio$a1*MAE) + (w_min_pct_cambio$a2*PBIAS) - (w_min_pct_cambio$a3*NSE)) %>%
  group_by(CodigoEstacion) %>% 
  slice_min(order_by = pct_cambio_, n = 1, with_ties = FALSE)

# join con rank_error para traer serie (método QM) 
df_error_ensemble = left_join(df_error, dplyr::select(rank_error, c('CodigoEstacion', 'serie')),
                              by = 'CodigoEstacion') %>% 
  rowwise() %>%
  mutate(d_ensemble_MAE = get(paste0(serie, "_MAE")),
         d_ensemble_PBIAS = get(paste0(serie, "_PBIAS")),
         d_ensemble_NSE = get(paste0(serie, "_NSE"))) %>%
  ungroup()

## gráfico pie del ensamble métodos QM ----------------------------------------
df_error_ensemble %>% 
  count(serie) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = serie)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            size = 9, color = "white") +
  labs(title = "Porcentaje de métodos QM en ensemble",
       fill = "Serie") +
  scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

## Boxplot MAE de ensamble métodos BC en ts original vs diff estacionalmente ----
df_MAE_long = df_error_ensemble %>% dplyr::select(MAE, d_MAE, d_ensemble_MAE) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "MAE") %>% 
  mutate(serie = factor(serie))

bx_ = bx_text(df_MAE_long, serie, MAE)
bx_pct_cambio = data.frame(serie = "d_ensemble_MAE",  pct_ = numeric(length("d_ensemble")))
bx_pct_cambio$pct_ = f_pct_cambio(df_error_ensemble, d_MAE, !!sym("d_ensemble_MAE"), perf = "MAE")

bx_pct_cambio = inner_join(dplyr::select(bx_$df_outliers, c(serie, max_r, upper_bound)),
                           bx_pct_cambio, by = 'serie') %>% 
  mutate(pct_ = round(pct_*100,2))


ggplot(df_MAE_long, aes(x = fct_rev(serie), y = MAE, fill = serie)) +
  geom_boxplot() +
  # texto % cambio en c/método qm
  geom_text(data = bx_pct_cambio, aes(label = paste0(pct_, "% \n cambio"), y = max(max_r)), 
            size = 4, color = "black", fontface = "bold") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min_r)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del MAE", x = "serie", y = "MAE") +
  theme_minimal() +
  #coord_cartesian(ylim = c(0, max(bx_$df_outliers$upper_bound)))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## Boxplot PBIAS de ensamble métodos BC en ts original vs diff estacionalmente ----
df_PBIAS_long = df_error_ensemble %>% dplyr::select(PBIAS, d_PBIAS, d_ensemble_PBIAS) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "PBIAS") %>% 
  mutate(serie = factor(serie))

bx_ = bx_text(df_PBIAS_long, serie, PBIAS)
bx_pct_cambio = data.frame(serie = "d_ensemble_PBIAS",  pct_ = numeric(length("d_")))
bx_pct_cambio$pct_ = f_pct_cambio(df_error_ensemble, d_PBIAS, !!sym("d_ensemble_PBIAS"), perf = "PBIAS")

bx_pct_cambio = inner_join(dplyr::select(bx_$df_outliers, c(serie, max_r, upper_bound)),
                           bx_pct_cambio, by = 'serie') %>% 
  mutate(pct_ = round(pct_*100,2))

ggplot(df_PBIAS_long, aes(x = fct_rev(serie), y = PBIAS, fill = serie)) +
  geom_boxplot() +
  # texto % cambio en c/método qm
  geom_text(data = bx_pct_cambio, aes(label = paste0(pct_, "% \n cambio"), y = max(upper_bound)), 
            size = 4, color = "black", fontface = "bold") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(min(bx_$df_outliers$lower_bound))), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del PBIAS", x = "serie", y = "PBIAS") +
  theme_minimal() +
  coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), abs(min(bx_$df_outliers$lower_bound))))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## Boxplot NSE de ensamble métodos BC en ts original vs diff estacionalmente ----
df_NSE_long = df_error_ensemble %>% dplyr::select(NSE, d_NSE, d_ensemble_NSE) %>% 
  pivot_longer(cols = everything(), names_to = "serie", values_to = "NSE") %>% 
  mutate(serie = factor(serie))

bx_ = bx_text(df_NSE_long, serie, NSE)
bx_pct_cambio = data.frame(serie = "d_ensemble_NSE",  pct_ = numeric(length("d_")))
bx_pct_cambio$pct_ = f_pct_cambio(df_error_ensemble, d_NSE, !!sym("d_ensemble_NSE"), perf = "NSE")

bx_pct_cambio = inner_join(dplyr::select(bx_$df_outliers, c(serie, max_r, upper_bound)),
                           bx_pct_cambio, by = 'serie') %>% 
  mutate(pct_ = round(pct_*100,2))

ggplot(df_NSE_long, aes(x = fct_rev(serie), y = NSE, fill = serie)) +
  geom_boxplot() +
  # texto % cambio en c/método qm
  geom_text(data = bx_pct_cambio, aes(label = paste0(pct_, "% \n cambio"), y = max(upper_bound)), 
            size = 4, color = "black", fontface = "bold") +
  # texto de % outliers
  geom_text(data = bx_$df_outliers, aes(label = paste0(outlier_pct, "% \n outliers"), y = min(lower_bound)), 
            size = 3, color = "black", fontface = "bold") +
  # texto quartiles y valores
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q1, label = paste0("Q1: ", round(Q1, 2))), 
            vjust = 1, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q2, label = paste0("Q2: ", round(Q2, 2))), 
            vjust = -0.5, size = 4) +
  geom_text(data = bx_$df_stats, aes(x = serie, y = Q3, label = paste0("Q3: ", round(Q3, 2))), 
            vjust = -1 , size = 4) +
  # título y nombres de ejes
  labs(title = "Boxplots del NSE", x = "serie", y = "NSE") +
  theme_minimal() +
  coord_cartesian(ylim = c(min(bx_$df_outliers$lower_bound), max(bx_$df_outliers$upper_bound)))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title=element_text(size=14),
        plot.title = element_text(size=14, face="bold")) 

## ----

for(i in 1:len_){
  # crea la col d_qm_ensamble con la serie que min las medidas de desmepño para esa serie
  pol2_list[[i]]$d_qm_ensamble = pol2_list[[i]] %>% 
    dplyr::select(df_error_ensemble$serie[i]) %>% pull()
}

df_cat2 <- function(list = pol2_list, cat_, perf){
  # Argumentos: ----------------------------------------------------------------
  # 
  # list lista con dfs con col v_cat, d_chirps.v2.0 y d_sttns
  # cat_ str nombre de la col en c/elemento list a filtrar
  # perf str con medidas de desempeño
  #----------------------------------------------------------------------------
  len_ = length(list)
  # deriva el vector de categorías de la col arg v_cat 
  v_cat <- distinct(pol2_list[[1]], !!sym(cat_)) %>% as_vector() #unique(list[[i]]$cat_mes)
  
  df_ = data.frame(matrix(numeric(len_), nc = length(v_cat), nr = len_)) %>%
    `colnames<-`(v_cat)
  
  for(i in 1:len_){
    # Calcula la medida de desempeño para c/categoría
    values_ <- sapply(v_cat, function(m){
      # Filtra los datos para la categoría 'm'
      data_mes <- filter(list[[i]], !!sym(cat_) == m)
      # Calcula la medida de desempeño de acuerdo al arg perf
      if(perf == 'MAE'){
        perf_value <- mae(sim = data_mes$d_qm_ensamble, 
                          obs = data_mes$d_sttns, na.rm = TRUE)
      }else if(perf == 'PBIAS'){
        perf_value <- pbias(sim = data_mes$d_qm_ensamble, 
                            obs = data_mes$d_sttns, na.rm = TRUE, dec=2)
      }else{
        perf_value <- NSE(sim = data_mes$d_qm_ensamble, 
                          obs = data_mes$d_sttns, na.rm = TRUE)
      }
      return(perf_value)
    })
    # guarda los resultados en el listado de R_Spearman
    df_[i,] <- values_
  }
  
  # crea la col CodigoEstacion
  df_$CodigoEstacion = sttns.points$CodigoEstacion[!sttns.points$DEPARTAMENTO %in% ADZ]
  
  df_ <- df_ %>%
    pivot_longer(cols = all_of(as.vector(v_cat)), 
                 names_to = cat_, 
                 values_to = perf)
  
  df_ = df_ %>% mutate(!!cat_ := factor(!!sym(cat_), levels = v_cat))
  return(df_)
}

df_MAE_mes2 <- df_cat2(list = pol2_list, cat_ = "mes", perf = "MAE")

df_MAE_mes_long = left_join(df_MAE_mes, rename(df_MAE_mes2, 'MAE_BC' = MAE),
                        by = c("CodigoEstacion", "mes")) %>% 
  pivot_longer(cols = c(MAE, MAE_BC),
               names_to = 'serie',
               values_to = 'MAE') %>% 
  mutate(serie = ifelse(serie == 'MAE', 'd_chirps', 'd_chirps_BC'))

#bx_ = bx_text(df_MAE_mes2, serie = mes, r = MAE)

ggplot(df_MAE_mes_long, aes(x = mes, y = MAE, fill = serie)) +
  geom_boxplot()

df_MAE_mes_polar = df_MAE_mes_long %>% group_by(mes, serie) %>%
  summarise(MAE = mean(MAE, na.rm = TRUE))

df_MAE_mes_polar$mes_num = rep(seq(1, 12, 1), each = 2)

df_MAE_mes_polar = df_MAE_mes_polar 
df_MAE_mes_polar_x= filter(df_MAE_mes_polar, mes_num == 12)
df_MAE_mes_polar_x$mes_num = 0

df_MAE_mes_polar_new <- rbind(df_MAE_mes_polar_x, df_MAE_mes_polar)

polar_plot <- ggplot(df_MAE_mes_polar, aes(x = mes_num, y = MAE, group = serie, color = serie)) +
  geom_line(linewidth = 2 ) +
  scale_x_continuous(breaks = 1:12) +
  theme_minimal()
  
MAE_max = ceiling(max(df_MAE_mes_polar$MAE) / 25) * 25

polar_plot %+%
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  labs(title = "MAE promedio por mes",x = "Mes", y = "MAE")+
  theme(legend.position = c(0.9, 0.95),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24),
        plot.title = element_text(size=24, face="bold"))

polar_plot %+%
  df_MAE_mes_polar_new +
  coord_polar(start = -pi/12) +
  scale_y_continuous(limits = c(0, MAE_max), breaks = seq(0, MAE_max, by = 25)) +
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  labs(title = "Polar plot del MAE promedio por mes",x = "Mes", y = "MAE") +
  theme(legend.position = c(0.95, 0.95),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24),
        plot.title = element_text(size=24, face="bold")) 

## Boxplot y lines PBIAS d_chirps vs d_qm_{ensamble} --------------------------
df_PBIAS_mes2 <- df_cat2(list = pol2_list, cat_ = "mes", perf = "PBIAS")

df_PBIAS_mes_long = left_join(df_PBIAS_mes, rename(df_PBIAS_mes2, 'PBIAS_BC' = PBIAS),
                            by = c("CodigoEstacion", "mes")) %>% 
  pivot_longer(cols = c(PBIAS, PBIAS_BC),
               names_to = 'serie',
               values_to = 'PBIAS') %>% 
  mutate(serie = ifelse(serie == 'PBIAS', 'd_chirps', 'd_chirps_BC'))

ggplot(df_PBIAS_mes_long, aes(x = mes, y = PBIAS, fill = serie)) +
  geom_boxplot()+
  scale_y_continuous(labels = percent_format(scale = 1), breaks = seq(-200, 200, 100)) +
  coord_cartesian(ylim = c(-300, 300))

df_PBIAS_mes_mean = df_PBIAS_mes_long %>% filter(abs(PBIAS)<1000) %>%  group_by(mes, serie) %>%
  summarise(PBIAS = mean(PBIAS, na.rm = TRUE))

ggplot(df_PBIAS_mes_mean, aes(x = mes, y = PBIAS, group = serie, color = serie)) +
  geom_line(linewidth = 2 ) +
  scale_x_discrete(labels = month.abb)+
  scale_y_continuous(breaks = seq(-50, 50, 25))+
  coord_cartesian(ylim = c(-50,50))+
  theme_minimal() +
  labs(title = "PBIAS promedio por mes",x = "Mes", y = "PBIAS")+
  theme(legend.position = c(0.9, 0.95),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24),
        plot.title = element_text(size=24, face="bold"))

## Boxplot y lines NSE d_chirps vs d_qm_{ensamble} --------------------------
df_NSE_mes2 <- df_cat2(list = pol2_list, cat_ = "mes", perf = "NSE")

df_NSE_mes_long = left_join(df_NSE_mes, rename(df_NSE_mes2, 'NSE_BC' = NSE),
                              by = c("CodigoEstacion", "mes")) %>% 
  pivot_longer(cols = c(NSE, NSE_BC),
               names_to = 'serie',
               values_to = 'NSE') %>% 
  mutate(serie = ifelse(serie == 'NSE', 'd_chirps', 'd_chirps_BC'))

ggplot(df_NSE_mes_long, aes(x = mes, y = NSE, fill = serie)) +
  geom_boxplot()+
  coord_cartesian(ylim= c(-3,1))


df_NSE_mes_mean = df_NSE_mes_long %>%  group_by(mes, serie) %>%
  summarise(NSE = mean(NSE, na.rm = TRUE))

ggplot(df_NSE_mes_mean, aes(x = mes, y = NSE, group = serie, color = serie)) +
  geom_line(linewidth = 2 ) +
  scale_x_discrete(labels = month.abb)+
  scale_y_continuous(breaks = seq(0, 1, 0.25))+
  coord_cartesian(ylim = c(0,1))+
  theme_minimal() +
  labs(title = "NSE promedio por mes",x = "Mes", y = "NSE")+
  theme(legend.position = c(0.9, 0.95),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24),
        plot.title = element_text(size=24, face="bold"))
## Violin plot con boxplot de MAE por región ----
bx_ = bx_text(subset(sttns.points, !is.na(d_MAE)), serie = region, r = d_MAE)

sttns.points = left_join(sttns.points, 
                         dplyr::select(df_error_ensemble, 
                                       c(CodigoEstacion, d_ensemble_MAE, d_ensemble_PBIAS, d_ensemble_NSE)),
                         by = 'CodigoEstacion')

df_MAE_region_long = sttns.points %>% dplyr::select(CodigoEstacion, d_MAE, d_ensemble_MAE, region) %>% 
  pivot_longer(cols = c(d_MAE, d_ensemble_MAE),
               names_to = 'serie',
               values_to = 'MAE') %>% 
  mutate(serie = ifelse(serie == 'd_MAE', 'd_chirps', 'd_chirps_BC'))

ggplot(df_MAE_region_long, aes(x = region, y = MAE, fill = serie)) +
  geom_boxplot()

df_PBIAS_region_long = sttns.points %>% dplyr::select(CodigoEstacion, d_PBIAS, d_ensemble_PBIAS, region) %>% 
  pivot_longer(cols = c(d_PBIAS, d_ensemble_PBIAS),
               names_to = 'serie',
               values_to = 'PBIAS') %>% 
  mutate(serie = ifelse(serie == 'd_PBIAS', 'd_chirps', 'd_chirps_BC'))

ggplot(df_PBIAS_region_long, aes(x = region, y = PBIAS, fill = serie)) +
  geom_boxplot()+
  scale_y_continuous(labels = percent_format(scale = 1), breaks = seq(-200, 200, 100)) +
  coord_cartesian(ylim = c(-300, 300))

df_NSE_region_long = sttns.points %>% dplyr::select(CodigoEstacion, d_NSE, d_ensemble_NSE, region) %>% 
  pivot_longer(cols = c(d_NSE, d_ensemble_NSE),
               names_to = 'serie',
               values_to = 'NSE') %>% 
  mutate(serie = ifelse(serie == 'd_NSE', 'd_chirps', 'd_chirps_BC'))

ggplot(df_NSE_region_long, aes(x = region, y = NSE, fill = serie)) +
  geom_boxplot()

# Check point 2 ----
#save.image('df_eval_metrics_ADF_QM.RData')
load('df_eval_metrics_ADF_QM.RData')
#load('df_eval_metrics.RData')

ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) +
  #geom_sf(data = filter(sttns.points, DEPARTAMENTO != ADZ), color = 'blue', size = 1) +
  geom_sf(data = filter(sttns.points, n ==480), color = 'red', size = 1) +
  #theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de sttns sin datos faltantes")   # Título y leyenda

# crea la col pct_missing
sttns.points = mutate(sttns.points, 
                      pct_missing = 1 - (n/480))

# selecciona sttns para validar métodos de imputación de datos faltantes -------
sttns_full = sttns.points %>%  filter(pct_missing ==0 & DEPARTAMENTO != ADZ) 
len_ = length(pol2_list)
lgic_81_10 = c()

for(i in 1:len_){
  pol2_list[[i]]$Date = as.Date(pol2_list[[i]]$Date)
  df_81_10 = filter(pol2_list[[i]], Date<as.Date("2011-01-01"))
  lgic_81_10[i] = sum(is.na(df_81_10$sttns))==0
}

ADZ = "Archipielago De San Andres, Providencia Y Santa Catalina"
sttns_full_81_10 = sttns.points %>% filter(DEPARTAMENTO != ADZ & CodigoEstacion %in% sapply(pol2_list[lgic_81_10], function(x) x$ID[1])) 

df_missing <- sttns.points %>% filter(pct_missing !=0 & DEPARTAMENTO != ADZ) %>% 
  dplyr::select(matches("^(19|20)"))

df_missing_81_10 <- sttns.points %>% filter(DEPARTAMENTO != ADZ & !CodigoEstacion %in% sapply(pol2_list[lgic_81_10], function(x) x$ID[1])) %>% 
  dplyr::select(matches("^(19|200|2010)"))

# identifica los indices de los NAs, para mapear c/patrón de datos faltantes
# es decir en lugar de tomar una muestra para simular los indices de datos faltantes
# se mapea c/patrón y eventualmente se toma una muestra o se usan todos los patrones
# para reemplazar en las ts completas por NAs y así estimar los valores con los métodos
# y así poder calcular las métricas de desempeño, y evaluar cuál es el mejor en Colombia
idx_missing <- apply(df_missing, 1, function(row) which(is.na(row)))
lgcl_len = sapply(idx_missing, function(x) length(x) == 0)
idx_missing2 = idx_missing[!lgcl_len]
idx_missing2 = idx_missing2[!duplicated(idx_missing2)]
length(idx_missing2) # 951 patrones de datos faltantes

idx_missing_81_10 <- apply(df_missing_81_10, 1, function(row) which(is.na(row)))
lgcl_len_81_10 = sapply(idx_missing_81_10, function(x) length(x) == 0)
idx_missing2_81_10 = idx_missing_81_10[!lgcl_len_81_10]
idx_missing2_81_10 = idx_missing2_81_10[!duplicated(idx_missing2_81_10)]
length(idx_missing2_81_10) # 902 patrones de datos faltantes

# 1 sttns con 1 patrón de datos faltantes
set.seed(01112025)
idx_sample = sample(1:length(idx_missing2), size = 1, replace = FALSE)

length(idx_missing2[[idx_sample]]) # 37 (7.7%) datos faltantes

# df_ts_full: data frame con sttns sin datos faltantes en periodo 1981-2020 ----
df_ts_full = sttns_full %>% dplyr::select(matches("^(19|20|CodigoEstacion)")) %>% 
  pivot_longer(cols = matches("^(19|20)"), names_to = "fecha", values_to = "pcp") %>% 
  as.data.frame() %>% mutate(CodigoEstacion = paste0('id_',CodigoEstacion)) %>% 
  pivot_wider(id_cols = c("fecha"), names_from = "CodigoEstacion", 
              values_from = "pcp")

df_CHIRPS_full = lapply(pol2_list[lgic_], function(x) dplyr::select(x, `chirps-v2-0`)) %>% 
  bind_cols() 

df_ts_full_81_10 = lapply(pol2_list[lgic_81_10], function(x) filter(x, Date < as.Date("2011-01-01")) %>% dplyr::select(sttns)) %>% 
  bind_cols() 

df_ts_full$fecha = as.Date(df_ts_full$fecha)

## convierte a objeto ts
ts_full = ts(dplyr::select(df_ts_full, -"fecha"), start = c(1981, 1), 
             frequency = 12, class = 'ts')

ts_CHIRPS_full = ts(df_CHIRPS_full, start = c(1981, 1), frequency = 12, class = 'ts')

ts_full_81_10 = ts(df_ts_full_81_10, start = c(1981, 1), 
                   frequency = 12, class = 'ts')

ts_ = ts_full # duplica objeto
ts_81_10 = ts_full_81_10 # duplica objeto
# !!! iterar sobre c/patrón de datos faltantes
ts_[idx_missing2[[idx_sample]],] <- NA # reemplaza valores de idx_sample por NA
#plot(ts_[,1:10])

# imputeTS Kalman (StructTS, auto.arima) ----
ts_impute = list()
#imputeTS::statsNA(ts_[,1]) # "Percentage of Missing Values:" /n "7.71%"
ts_impute[["na.kalman"]] = imputeTS::na_kalman(ts_, model = "StructTS")
#ts_impute[["na.kalman_auto.arima"]] = imputeTS::na_kalman(ts_, model = "auto.arima")

#ts_impute[["na.mean"]] = imputeTS::na_mean(ts_, option = "mean")
# cambiar la ts_ por d_ts_ diff estacionalmente para remover comp estacional
#ts_impute[["na.seadec"]] = imputeTS::na_seadec(ts_, algorithm = "kalman", model = "auto.arima") # desc basada en STL::LOESS
# direct from CHIRPS_Corregida

# métodos basados en distancias -----
## solo se puede aplicar al 40% (69% si también se considera 1-NN)
m_distance = st_distance(sttns_full$geometry) # obtiene matrix distancias
m_distance = set_units(m_distance, "km") # convierte a km   
m_distance = matrix(as.numeric(m_distance), nr = nrow(m_distance), ncol = ncol(m_distance))
#m_distance[m_distance>100 | m_distance == 0] <- NA
m_distance[m_distance>50 | m_distance == 0] <- NA # remueve distancias >50km 
table(apply(m_distance,2,function(x) sum(!is.na(x))))

m_distance_81_10 = st_distance(sttns_full_81_10$geometry) # obtiene matrix distancias
m_distance_81_10 = set_units(m_distance_81_10, "km") # convierte a km   
m_distance_81_10 = matrix(as.numeric(m_distance_81_10), nr = nrow(m_distance_81_10), ncol = ncol(m_distance_81_10))
#m_distance[m_distance>100 | m_distance == 0] <- NA
m_distance_81_10[m_distance_81_10>50 | m_distance_81_10 == 0] <- NA # remueve distancias >50km 
table(apply(m_distance_81_10,2,function(x) sum(!is.na(x))))

## obtiene idx de las posiciones no-nulas
l_idx_NN = apply(m_distance, 2 , function(x) which(!is.na(x)))
l_idx_NN_81_10 = apply(m_distance_81_10, 2 , function(x) which(!is.na(x))) 
#names(l_idx_NN) = sttns_full$CodigoEstacion
idw_ = apply(m_distance,2, function(x) x^(-2)/sum(x^-(2), na.rm = TRUE)) # reemplaza por x^-j / sum x^-j
idw_81_10 = apply(m_distance_81_10,2, function(x) x^(-2)/sum(x^-(2), na.rm = TRUE)) # reemplaza por x^-j / sum x^-j

# mean annual pcp data
NR_full = apply(dplyr::select(df_ts_full, -fecha), 2, function(x) sum(x)/40)
NR_ = apply(ts_, 2, function(x) sum(x, na.rm = T)/40) # se divide por 40 (años)

NR_full_81_10 = apply(df_ts_full_81_10, 2, function(x) sum(x)/30)
NR_81_10 = apply(ts_81_10, 2, function(x) sum(x, na.rm = T)/30) # se divide por 30 (años)

# métodos basados en correlación --------------------------------------------
# missForest
# CCW

## calcular el coef de corr con ts diff estacionalmente
## df_d_ts_full: data frame con las d_ts completas

## filtra sttns sin datos faltantes

lgic_ = c()
for(i in 1:len_){
  lgic_[i] = sum(is.na(pol2_list[[i]]$sttns))==0
}

max_ = max(df_pgram$period)
df_d_ts_full = lapply(pol2_list[lgic_], function(x) dplyr::select(x, d_sttns)) %>% 
  bind_cols() %>% mutate(t_ = 1:nrow(.)) %>% 
  # filtra las obs iniciales, dado que las ts tienen periodos diferentes
  filter(t_ > max_) %>% dplyr::select(-t_)

df_d_ts_full_81_10 = lapply(pol2_list[lgic_81_10], function(x) filter(x, Date < as.Date("2011-01-01")) %>% dplyr::select(d_sttns)) %>% 
  bind_cols() %>% mutate(t_ = 1:nrow(.)) %>% 
  # filtra las obs iniciales, dado que las ts tienen periodos diferentes
  filter(t_ > max_) %>% dplyr::select(-t_)

df_d_CHIRPS = lapply(pol2_list[lgic_], function(x) dplyr::select(x, d_qm_ensamble)) %>% 
  bind_cols()  %>% mutate(t_ = 1:nrow(.)) %>% 
  # filtra las obs iniciales, dado que las ts tienen periodos diferentes
  filter(t_ > max_) %>% dplyr::select(-t_)

df_d_CHIRPS_81_10 = lapply(pol2_list[lgic_81_10], function(x) filter(x, Date < as.Date("2011-01-01")) %>%  dplyr::select(d_qm_ensamble)) %>% 
  bind_cols()  %>% mutate(t_ = 1:nrow(.)) %>% 
  # filtra las obs iniciales, dado que las ts tienen periodos diferentes
  filter(t_ > max_) %>% dplyr::select(-t_)

df_d_missing = lapply(pol2_list[!lgic_], function(x) dplyr::select(x, d_sttns)) %>% 
  bind_cols() %>% mutate(t_ = 1:nrow(.)) %>% 
  # filtra las obs iniciales, dado que las ts tienen periodos diferentes
  filter(t_ > max_) %>% dplyr::select(-t_)

df_d_missing_81_10 = lapply(pol2_list[!lgic_81_10], function(x) filter(x, Date < as.Date("2011-01-01")) %>% dplyr::select(d_sttns)) %>% 
  bind_cols() %>% mutate(t_ = 1:nrow(.)) %>% 
  # filtra las obs iniciales, dado que las ts tienen periodos diferentes
  filter(t_ > max_) %>% dplyr::select(-t_)

## obtiene los patrones unicos de datos faltantes
idx_d_missing <- apply(df_d_missing, 2, function(x) which(is.na(x)))
lgcl_len = sapply(idx_d_missing, function(x) length(x) == 0)
idx_d_missing2 = idx_d_missing[!lgcl_len]
idx_d_missing2 = idx_d_missing2[!duplicated(idx_d_missing2)] # 953 patrones de datos faltantes

idx_d_missing_81_10 <- apply(df_d_missing_81_10, 2, function(x) which(is.na(x)))
lgcl_len_81_10 = sapply(idx_d_missing_81_10, function(x) length(x) == 0)
idx_d_missing2_81_10 = idx_d_missing_81_10[!lgcl_len_81_10]
idx_d_missing2_81_10 = idx_d_missing2_81_10[!duplicated(idx_d_missing2_81_10)] # 951 patrones de datos faltantes

table(sapply(idx_d_missing2_81_10, length))

## d_ts_full: ts con las series completas
d_ts_full = ts(df_d_ts_full, start = c(1982, 1), frequency = 12, class = 'ts')
d_CHIRPS = ts(df_d_CHIRPS, start = c(1982, 1), frequency = 12, class = 'ts')

d_ts_full_81_10 = ts(df_d_ts_full_81_10, start = c(1982, 1), frequency = 12, class = 'ts')
d_CHIRPS_81_10 = ts(df_d_CHIRPS_81_10, start = c(1982, 1), frequency = 12, class = 'ts')

## m_corr: matriz de corr Spearman para las ts diff estacionalmente
m_corr = cor(d_ts_full, method = "spearman", use = "complete.obs")
m_corr[m_corr<0.5 | m_corr == 1] <- NA # remueve coef corr < 0.5

m_corr_81_10 = cor(d_ts_full_81_10, method = "spearman", use = "complete.obs")
m_corr_81_10[m_corr_81_10<0.6 | m_corr_81_10 == 1] <- NA # remueve coef corr < 0.6

## selecciona las 3 sttn con coef corr más altos
max_ = max(apply(m_corr, 2, function(x) sum(!is.na(x))))
max_81_10 = max(apply(m_corr_81_10, 2, function(x) sum(!is.na(x))))
# mientras el max de vecinos sea mayor a 3
while(max_>3){for(i in 1:ncol(m_corr)){
  # si la columna tiene mas de 3 vecinos  
  if(sum(!is.na(m_corr[,i]))>3){
    # si el coef corr es menor o igual al minimo hacerlo NA 
    lgcl_ = m_corr[,i] <= min(m_corr[,i], na.rm = T)
    m_corr[lgcl_, i] <- NA
  }
}
  # recalcula el número max de vecinos
  max_ = max(apply(m_corr, 2, function(x) sum(!is.na(x))))
}

while(max_81_10>3){for(i in 1:ncol(m_corr_81_10)){
  # si la columna tiene mas de 3 vecinos  
  if(sum(!is.na(m_corr_81_10[,i]))>3){
    # si el coef corr es menor o igual al minimo hacerlo NA 
    lgcl_ = m_corr_81_10[,i] <= min(m_corr_81_10[,i], na.rm = T)
    m_corr_81_10[lgcl_, i] <- NA
  }
}
  # recalcula el número max de vecinos
  max_81_10 = max(apply(m_corr_81_10, 2, function(x) sum(!is.na(x))))
}

table(apply(m_corr, 2 , function(x) sum(!is.na(x))))
table(apply(m_corr_81_10, 2 , function(x) sum(!is.na(x))))

# lista con los indices de los coef corr más altos
l_idx_corr = apply(m_corr, 2 , function(x) which(!is.na(x))) # obtiene idx de las posiciones no-nulas
l_idx_corr_81_10 = apply(m_corr_81_10, 2 , function(x) which(!is.na(x))) # obtiene idx de las posiciones no-nulas
# pesos de c/coef corr
ccw_ = apply(m_corr, 2, function(x) x/sum(x, na.rm = T))
ccw_81_10 = apply(m_corr_81_10, 2, function(x) x/sum(x, na.rm = T))

d_NR_full = apply(df_d_ts_full, 2, function(x) sum(x)/39)
d_NR_full_81_10 = apply(df_d_ts_full_81_10, 2, function(x) sum(x)/29)

# sim_impute: función optimizada paraa los valores de idx_
sim_impute <- function(M_ = m_distance, M_c = m_corr, ts_ = d_ts_full, 
                       idx_ = idx_d_missing2[[i]], d_NR_full = d_NR_full,
                       ts_CHIRPS = d_CHIRPS, wdw_30 = FALSE ){
  #############################################################################
  # m_distance: matriz de distancias de las sttns con datos completos
  # m_corr: matriz de correlaciones de las sttns con datos completos
  # d_ts_full: ts sttns con datos completos
  # idx_d_missing2: lista patrones de datos faltantes (954 elementos distintos)
  # d_NR_full: normales anuales de las sttns con datos completos
  # wdw_30: bool indica sí se aplica la ventana de 30 años a las Normales anuales
  # output_: lista con la mediana de c/medidas de desempeño para el patrón idx_ en las sttns con datos completos
  #############################################################################
  
  ## l_idx_d_NN lista a partir de M_ sttns con distancias < 50km
  l_idx_d_NN = apply(M_, 2 , function(x) which(!is.na(x))) 
  
  ## l_idx_corr lista a partir de M_c sttns con coef corr Spearman > 0.5
  l_idx_corr = apply(M_c, 2 , function(x) which(!is.na(x))) # obtiene idx de las posiciones no-nulas
  
  ## idw_ lista con pesos idw x^-j / sum x^-j para sttns de M_
  idw_ = apply(M_, 2, function(x) x^(-2)/sum(x^-(2), na.rm = TRUE)) 
  
  ## ccw_ lista con pesos coef corr Spearman
  ccw_ = apply(M_c, 2, function(x) x/sum(x, na.rm = T))
  
  ## d_ts_ reemplaza valores por nulos según patrón de índices del arg idx_
  d_ts_ = ts_ # duplica objeto
  d_ts_[idx_,] <- NA # reemplaza valores de idx_sample por NA
  
  ## d_NR_ normales anuales de la sttns con datos faltantes
  if(wdw_30){
    d_NR_ = apply(d_ts_, 2, function(x) sum(x, na.rm = T)/29) # se divide por 39 (años, se removió el año 1 por diff estacionales)
  }else{
    #d_NR_ = apply(d_ts_, 2, function(x) sum(x, na.rm = T)/39) # se divide por 39 (años, se removió el año 1 por diff estacionales)
    d_NR_ = apply(d_ts_, 2, function(x) sum(x, na.rm = T)/40) # se divide por 39 (años, se removió el año 1 por diff estacionales)
  }
  
  ## l_idx_d_NN agrega 4 elementos a la lista a partir de los índices de los elementos no-nulos de M_
  len_ = length(l_idx_d_NN)
  len_2 = length(l_idx_corr)
  
  for(i in 1:len_){
    l_ = length(l_idx_d_NN[[i]]) # longitud de c/elemento de la lista
    # si hay al menos 1 elemento (sttns vecinas)
    if(l_ > 0){
      ## crea una lista con: 
      ## 1. d_ts_ cols de las d_ts vecinas completas, 
      ## 2. IDW_ pesos para c/ts vecina
      ## 3. d_NR_ razón valores normales anuales
      ## 4. AA_ pesos para promedio ajustable a si cumple criterio "within 10%"
      l_idx_d_NN[[i]] = list(d_ts_ = ts_[,l_idx_d_NN[[i]]],
                             IDW_ = idw_[!is.na(idw_[,i]),i],
                             # pesos N_x/N_i que cumplen N_i excede 10% a N_x
                             d_NR_ = ifelse((d_NR_full[l_idx_d_NN[[i]]] - d_NR_[i])/d_NR_[i] > 0.1, 
                                            d_NR_[i]/d_NR_full[l_idx_d_NN[[i]]], 0 ),
                             # promedio aritmético con sttns que cumpen estar dentro del 10% de N_x
                             AA_ = ifelse(abs((d_NR_full[l_idx_d_NN[[i]]] - d_NR_[i])/d_NR_[i]) <= 0.1, 1, 0)
      )
    }          
  }
  
  ## l_idx_corr agrega elementos a la lista a partir de los índices de elementos no-nulos de M_c
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]]) # longitud de c/elemento de la lista
    # si hay al menos 1 elemento (sttns con corr_Spearman > 0.5)
    if(l_ > 0){
      ## crea una lista con: 
      ## 1.las cols de las d_ts con mayor corr_Speraman completas, 
      ## 2. los pesos para c/ts segun el corr_Spearman
      l_idx_corr[[i]] = list(d_ts_ = ts_[,l_idx_corr[[i]]],
                             CCW_ = ccw_[!is.na(ccw_[,i]),i])
    }
  }
  
  ## l_idx_corr agrega la ts con missings a las ts de sttns con corr_Spearman > 0.5
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]])
    # si l_ > 0 indica que hay al menos 1 sttn con coef corr > 0.5
    if(l_ > 0){
      l_idx_corr[[i]]$d_ts_RF = cbind(d_ts_[,i], l_idx_corr[[i]]$d_ts_)
    }
  }
  
  ## reemplaza los pesos/razón valores normales por valores imputados
  for(i in 1:len_){
    l_ = length(l_idx_d_NN[[i]]) # longitud c/elemento de la lista
    # si l_ > 0 indica que hay al menos 1 sttn cercana
    if(l_ > 0){
      idw_ = rep(0, length(idx_))
      aa_ = rep(0, length(idx_))
      nr_ = rep(0, length(idx_))
      l_2 = length(l_idx_d_NN[[i]]$IDW_) # num de sttns vecinas
      
      for(j in 1:l_2){
        # si l_2 == 1 solo hay 1 sttn
        if(l_2 == 1){
          # si sólo hay 1 sttn vecina, el valor a imputar es el de esa sttn vecina
          # si el valor de esa sttn es NA, quedaría como NA
          idw_ = l_idx_d_NN[[i]]$d_ts_[idx_]
          aa_ = l_idx_d_NN[[i]]$d_ts_[idx_]
          nr_ = l_idx_d_NN[[i]]$d_ts_[idx_]
        }else{
          # inverse distance weight
          idw_ = idw_ + l_idx_d_NN[[i]]$IDW_[j]*replace(l_idx_d_NN[[i]]$d_ts_[idx_, j],
                                                        is.na(l_idx_d_NN[[i]]$d_ts_[idx_, j]), 0)
          # media aritmetica 
          aa_ = aa_ + (l_idx_d_NN[[i]]$AA_[j]/sum(l_idx_d_NN[[i]]$AA_))*replace(l_idx_d_NN[[i]]$d_ts_[idx_, j],
                                                                                is.na(l_idx_d_NN[[i]]$d_ts_[idx_, j]), 0)
          # Normal ratio
          nr_ = nr_ + l_idx_d_NN[[i]]$d_NR_[j]*replace(l_idx_d_NN[[i]]$d_ts_[idx_, j],
                                                       is.na(l_idx_d_NN[[i]]$d_ts_[idx_, j]), 0)
        }
      }
      ## asigna valores imputados
      l_idx_d_NN[[i]]$IDW_ = idw_ # con IDW
      l_idx_d_NN[[i]]$AA_ = aa_ # con AA
      l_idx_d_NN[[i]]$d_NR_ = nr_ # con NR
    }
  }
  
  ## reemplaza los pesos por los valores imputados
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]]) # longitud de c/elemento de la lista
    # si l_ > 0 indica que hay al menos 1 elemento (sttns con corr_Spearman > 0.5)
    if(l_ > 0){
      ccw_ = rep(0, length(idx_))
      l_2 = length(l_idx_corr[[i]]$CCW_) # num de sttns con corr_Separman > 0.5
      
      for(j in 1:l_2){
        # si l_2 == 1 solo hay 1 sttn
        if(l_2 == 1){
          # si sólo hay 1 sttn vecina, el valor a imputar es el de esa sttn
          # si el valor de esa sttn es NA, quedaría como NA
          ccw_ = l_idx_corr[[i]]$d_ts_[idx_]
        }else{
          # en caso que alguno de los valores sea NA e considera replace
          ccw_ = ccw_ + l_idx_corr[[i]]$CCW_[j]*replace(l_idx_corr[[i]]$d_ts_[idx_, j],
                                                        is.na(l_idx_corr[[i]]$d_ts_[idx_, j]), 0)
        }
      }
      # crea 2 elementos en la lista
      l_idx_corr[[i]]$CCW_ = ccw_ # valores imputados con ccW
    }
  }
  
  ## agrega el elemento mF_ con los valores imputados 
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]])
    # si l_ > 0 indica que hay al menos 1 sttn con coef corr > 0.5
    if(l_ > 0){
      l_idx_corr[[i]]$mF_ = missForest(l_idx_corr[[i]]$d_ts_RF, backend = "randomForest")$ximp[idx_, 1] # selecciona la col 1 (que contiene datos faltantes)
    }
  }
  
  ## d_ts_impute lista para guardar los valores imputados 
  d_ts_impute = list()
  len_idx = length(idx_)
  
  d_ts_impute[["na.kalman"]] = imputeTS::na_kalman(d_ts_, model = "StructTS") %>% as.data.frame()
  d_ts_impute$na.kalman = d_ts_impute$na.kalman[idx_,] # {len_idx} x {48}
  
  d_ts_impute[["sttn"]] = data.frame(matrix(ts_[idx_,], nr = len_idx)) #%>% 
    #`colnames<-`(names(d_ts_impute$na.kalman)) # {len_idx} x {48}
  
  d_ts_impute[["CHIRPS_C"]] = data.frame(matrix(ts_CHIRPS[idx_,], nr = len_idx))
  
  d_ts_impute[["CCW"]] = lapply(l_idx_corr, function(x) 
    if(length(x)>0){x$CCW_}else{rep(NA,len_idx)}) %>% 
    bind_cols() # {len_idx} x {48}
  
  d_ts_impute[["IDW"]] = lapply(l_idx_d_NN, function(x) 
    if(length(x)>0){x$IDW_}else{rep(NA,len_idx)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  d_ts_impute[["AA"]] = lapply(l_idx_d_NN, function(x) 
    if(length(x)>0){x$AA_}else{rep(NA,len_idx)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  d_ts_impute[["NR"]] = lapply(l_idx_d_NN, function(x) 
    if(length(x)>0){x$d_NR_}else{rep(NA,len_idx)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  d_ts_impute[["mF"]] = lapply(l_idx_corr, function(x) 
    if(length(x)>0){x$mF_}else{rep(NA,len_idx)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  ## métricas desempeño -----
  ## MAE 
  d_impute_MAE = list()
  ## NSE 
  d_impute_NSE = list()
  ## R_Spearman
  d_impute_R_s = list()
  
  d_impute_MAE[["na.kalman"]] = mae(sim = d_ts_impute$na.kalman,
                                    obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["IDW"]] = mae(sim = d_ts_impute$IDW,
                              obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["AA"]] = mae(sim = d_ts_impute$AA,
                             obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["NR"]] = mae(sim = d_ts_impute$NR,
                             obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["CCW"]] = mae(sim = d_ts_impute$CCW,
                              obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["mF"]] = mae(sim = d_ts_impute$mF,
                             obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["CHIRPS_C"]] = mae(sim = d_ts_impute$CHIRPS_C,
                                   obs = d_ts_impute$sttn, na.rm=TRUE)
  
  if(len_idx > 1){
    d_impute_NSE[["na.kalman"]] = NSE(sim = d_ts_impute$na.kalman,
                                      obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["IDW"]] = NSE(sim = d_ts_impute$IDW,
                                obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["AA"]] = NSE(sim = d_ts_impute$AA,
                               obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["NR"]] = NSE(sim = d_ts_impute$NR,
                               obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["CCW"]] = NSE(sim = d_ts_impute$CCW,
                                obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["mF"]] = NSE(sim = d_ts_impute$mF,
                               obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["CHIRPS_C"]] = NSE(sim = d_ts_impute$CHIRPS_C,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
  }
  
  if(len_idx > 2){
    
    d_impute_R_s[["na.kalman"]] = rSpearman(sim = d_ts_impute$na.kalman,
                                            obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["IDW"]] = rSpearman(sim = d_ts_impute$IDW,
                                      obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["AA"]] = rSpearman(sim = d_ts_impute$AA,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["NR"]] = rSpearman(sim = d_ts_impute$NR,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["CCW"]] = rSpearman(sim = d_ts_impute$CCW,
                                      obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["mF"]] = rSpearman(sim = d_ts_impute$mF,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["CHIRPS_C"]] = rSpearman(sim = d_ts_impute$CHIRPS_C,
                                           obs = d_ts_impute$sttn, na.rm=TRUE)
  }
  
  output_ = list()
  output_[["MAE"]] = lapply(d_impute_MAE, function(x) quantile(x, probs = 0.5, na.rm = TRUE) %>% as.vector())
  if(len_idx >1){
    output_[["r_s"]] = lapply(d_impute_R_s, function(x) quantile(x, probs = 0.5, na.rm = TRUE) %>% as.vector())
    output_[["NSE"]] = lapply(d_impute_NSE, function(x) quantile(x, probs = 0.5, na.rm = TRUE) %>% as.vector())
  }
  
  return(output_) 
}

lags_full = df_pgram$period[lgic_]

pad_diff <- function(x, lag_) {
  d <- diff(x, lag = lag_)
  c(rep(NA, lag_), d)     # rellenar con NAs al inicio
}

sim_impute_full <- function(M_ = m_distance, M_c = m_corr, ts_ = d_ts_full, 
                            idx_ = idx_d_missing2[[i]], d_NR_full = d_NR_full,
                            ts_CHIRPS = d_CHIRPS, lags_ = lags_full, wdw_30 = FALSE){
  #############################################################################
  # m_distance: matriz de distancias de las sttns con datos completos
  # m_corr: matriz de correlaciones de las sttns con datos completos
  # d_ts_full: ts sttns con datos completos
  # idx_d_missing2: lista patrones de datos faltantes (954 elementos distintos)
  # d_NR_full: normales anuales de las sttns con datos completos
  # wdw_30: bool indica sí se aplica la ventana de 30 años a las Normales anuales
  # output_: lista con la mediana de c/medidas de desempeño para el patrón idx_ en las sttns con datos completos
  #############################################################################
  
  ## l_idx_d_NN lista a partir de M_ sttns con distancias < 50km
  l_idx_d_NN = apply(M_, 2 , function(x) which(!is.na(x))) 
  
  ## l_idx_corr lista a partir de M_c sttns con coef corr Spearman > 0.5
  l_idx_corr = apply(M_c, 2 , function(x) which(!is.na(x))) # obtiene idx de las posiciones no-nulas
  
  ## idw_ lista con pesos idw x^-j / sum x^-j para sttns de M_
  idw_ = apply(M_, 2, function(x) x^(-2)/sum(x^-(2), na.rm = TRUE)) 
  
  ## ccw_ lista con pesos coef corr Spearman
  ccw_ = apply(M_c, 2, function(x) x/sum(x, na.rm = T))
  
  ## d_ts_ reemplaza valores por nulos según patrón de índices del arg idx_
  d_ts_ = ts_ # duplica objeto
  d_ts_[idx_,] <- NA # reemplaza valores de idx_sample por NA
  
  ## d_NR_ normales anuales de la sttns con datos faltantes
  if(wdw_30){
    d_NR_ = apply(d_ts_, 2, function(x) sum(x, na.rm = T)/29) # se divide por 39 (años, se removió el año 1 por diff estacionales)
  }else{
    #d_NR_ = apply(d_ts_, 2, function(x) sum(x, na.rm = T)/39) # se divide por 39 (años, se removió el año 1 por diff estacionales)
    d_NR_ = apply(d_ts_, 2, function(x) sum(x, na.rm = T)/40) # se divide por 39 (años, se removió el año 1 por diff estacionales)
  }
  
  ## l_idx_d_NN agrega 4 elementos a la lista a partir de los índices de los elementos no-nulos de M_
  len_ = length(l_idx_d_NN)
  len_2 = length(l_idx_corr)
  
  for(i in 1:len_){
    l_ = length(l_idx_d_NN[[i]]) # longitud de c/elemento de la lista
    # si hay al menos 1 elemento (sttns vecinas)
    if(l_ > 0){
      ## crea una lista con: 
      ## 1. d_ts_ cols de las d_ts vecinas completas, 
      ## 2. IDW_ pesos para c/ts vecina
      ## 3. d_NR_ razón valores normales anuales
      ## 4. AA_ pesos para promedio ajustable a si cumple criterio "within 10%"
      l_idx_d_NN[[i]] = list(d_ts_ = ts_[,l_idx_d_NN[[i]]],
                             IDW_ = idw_[!is.na(idw_[,i]),i],
                             # pesos N_x/N_i que cumplen N_i excede 10% a N_x
                             d_NR_ = ifelse((d_NR_full[l_idx_d_NN[[i]]] - d_NR_[i])/d_NR_[i] > 0.1, 
                                            d_NR_[i]/d_NR_full[l_idx_d_NN[[i]]], 0 ),
                             # promedio aritmético con sttns que cumpen estar dentro del 10% de N_x
                             AA_ = ifelse(abs((d_NR_full[l_idx_d_NN[[i]]] - d_NR_[i])/d_NR_[i]) <= 0.1, 1, 0)
      )
    }          
  }
  
  ## l_idx_corr agrega elementos a la lista a partir de los índices de elementos no-nulos de M_c
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]]) # longitud de c/elemento de la lista
    # si hay al menos 1 elemento (sttns con corr_Spearman > 0.5)
    if(l_ > 0){
      ## crea una lista con: 
      ## 1.las cols de las d_ts con mayor corr_Speraman completas, 
      ## 2. los pesos para c/ts segun el corr_Spearman
      l_idx_corr[[i]] = list(d_ts_ = ts_[,l_idx_corr[[i]]],
                             CCW_ = ccw_[!is.na(ccw_[,i]),i])
    }
  }
  
  ## l_idx_corr agrega la ts con missings a las ts de sttns con corr_Spearman > 0.5
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]])
    # si l_ > 0 indica que hay al menos 1 sttn con coef corr > 0.5
    if(l_ > 0){
      l_idx_corr[[i]]$d_ts_RF = cbind(d_ts_[,i], l_idx_corr[[i]]$d_ts_)
    }
  }
  
  ## reemplaza los pesos/razón valores normales por valores imputados
  for(i in 1:len_){
    l_ = length(l_idx_d_NN[[i]]) # longitud c/elemento de la lista
    # si l_ > 0 indica que hay al menos 1 sttn cercana
    if(l_ > 0){
      idw_ = rep(0, length(idx_))
      aa_ = rep(0, length(idx_))
      nr_ = rep(0, length(idx_))
      l_2 = length(l_idx_d_NN[[i]]$IDW_) # num de sttns vecinas
      
      for(j in 1:l_2){
        # si l_2 == 1 solo hay 1 sttn
        if(l_2 == 1){
          # si sólo hay 1 sttn vecina, el valor a imputar es el de esa sttn vecina
          # si el valor de esa sttn es NA, quedaría como NA
          idw_ = l_idx_d_NN[[i]]$d_ts_#[idx_]
          aa_ = l_idx_d_NN[[i]]$d_ts_#[idx_]
          nr_ = l_idx_d_NN[[i]]$d_ts_#[idx_]
        }else{
          # inverse distance weight
          #idw_ = idw_ + l_idx_d_NN[[i]]$IDW_[j]*replace(l_idx_d_NN[[i]]$d_ts_[idx_, j],
          #                                              is.na(l_idx_d_NN[[i]]$d_ts_[idx_, j]), 0)
          idw_ = idw_ + l_idx_d_NN[[i]]$IDW_[j]*replace(l_idx_d_NN[[i]]$d_ts_[, j],
                                                        is.na(l_idx_d_NN[[i]]$d_ts_[, j]), 0)
          # media aritmetica 
          #aa_ = aa_ + (l_idx_d_NN[[i]]$AA_[j]/sum(l_idx_d_NN[[i]]$AA_))*replace(l_idx_d_NN[[i]]$d_ts_[idx_, j],
          #                                                                      is.na(l_idx_d_NN[[i]]$d_ts_[idx_, j]), 0)
          aa_ = aa_ + (l_idx_d_NN[[i]]$AA_[j]/sum(l_idx_d_NN[[i]]$AA_))*replace(l_idx_d_NN[[i]]$d_ts_[, j],
                                                                                is.na(l_idx_d_NN[[i]]$d_ts_[, j]), 0)
          # Normal ratio
          #nr_ = nr_ + l_idx_d_NN[[i]]$d_NR_[j]*replace(l_idx_d_NN[[i]]$d_ts_[idx_, j],
          #                                             is.na(l_idx_d_NN[[i]]$d_ts_[idx_, j]), 0)
          nr_ = nr_ + l_idx_d_NN[[i]]$d_NR_[j]*replace(l_idx_d_NN[[i]]$d_ts_[, j],
                                                       is.na(l_idx_d_NN[[i]]$d_ts_[, j]), 0)
        }
      }
      ## asigna valores imputados
      l_idx_d_NN[[i]]$IDW_ = idw_ # con IDW
      l_idx_d_NN[[i]]$AA_ = aa_ # con AA
      l_idx_d_NN[[i]]$d_NR_ = nr_ # con NR
    }
  }
  
  ## reemplaza los pesos por los valores imputados
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]]) # longitud de c/elemento de la lista
    # si l_ > 0 indica que hay al menos 1 elemento (sttns con corr_Spearman > 0.5)
    if(l_ > 0){
      ccw_ = rep(0, length(idx_))
      l_2 = length(l_idx_corr[[i]]$CCW_) # num de sttns con corr_Separman > 0.5
      
      for(j in 1:l_2){
        # si l_2 == 1 solo hay 1 sttn
        if(l_2 == 1){
          # si sólo hay 1 sttn vecina, el valor a imputar es el de esa sttn
          # si el valor de esa sttn es NA, quedaría como NA
          ccw_ = l_idx_corr[[i]]$d_ts_#[idx_]
        }else{
          # en caso que alguno de los valores sea NA e considera replace
          #ccw_ = ccw_ + l_idx_corr[[i]]$CCW_[j]*replace(l_idx_corr[[i]]$d_ts_[idx_, j],
          #                                              is.na(l_idx_corr[[i]]$d_ts_[idx_, j]), 0)
          ccw_ = ccw_ + l_idx_corr[[i]]$CCW_[j]*replace(l_idx_corr[[i]]$d_ts_[, j],
                                                        is.na(l_idx_corr[[i]]$d_ts_[, j]), 0)
        }
      }
      # crea 2 elementos en la lista
      l_idx_corr[[i]]$CCW_ = ccw_ # valores imputados con ccW
    }
  }
  
  ## agrega el elemento mF_ con los valores imputados 
  for(i in 1:len_2){
    l_ = length(l_idx_corr[[i]])
    # si l_ > 0 indica que hay al menos 1 sttn con coef corr > 0.5
    if(l_ > 0){
      #l_idx_corr[[i]]$mF_ = missForest(l_idx_corr[[i]]$d_ts_RF, backend = "randomForest")$ximp[idx_, 1] # selecciona la col 1 (que contiene datos faltantes)
      l_idx_corr[[i]]$mF_ = missForest(l_idx_corr[[i]]$d_ts_RF, backend = "randomForest")$ximp[,1] # selecciona la col 1 (que contiene datos faltantes)
    }
  }
  
  ## d_ts_impute lista para guardar los valores imputados 
  d_ts_impute = list()
  len_idx = length(idx_)
  len_ts = nrow(d_ts_)
  
  d_ts_impute[["na.kalman"]] = imputeTS::na_kalman(d_ts_, model = "StructTS") %>% as.data.frame()
  #d_ts_impute$na.kalman = d_ts_impute$na.kalman #d_ts_impute$na.kalman[idx_,] # {len_idx} x {48}
  
  #d_ts_impute[["sttn"]] = data.frame(matrix(ts_[idx_,], nr = len_idx)) #%>% 
  d_ts_impute[["sttn"]] = data.frame(matrix(ts_, nr = len_ts)) #%>% 
  #`colnames<-`(names(d_ts_impute$na.kalman)) # {len_idx} x {48}
  
  #d_ts_impute[["CHIRPS_C"]] = data.frame(matrix(ts_CHIRPS[idx_,], nr = len_idx))
  d_ts_impute[["CHIRPS_C"]] = data.frame(matrix(ts_CHIRPS, nr = len_ts))
  
  d_ts_impute[["CCW"]] = lapply(l_idx_corr, function(x) 
    if(length(x)>0){x$CCW_}else{rep(NA,len_ts)}) %>% 
    bind_cols() # {len_idx} x {48}
  
  d_ts_impute[["IDW"]] = lapply(l_idx_d_NN, function(x) 
    if(length(x)>0){x$IDW_}else{rep(NA,len_ts)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  d_ts_impute[["AA"]] = lapply(l_idx_d_NN, function(x) 
    if(length(x)>0){x$AA_}else{rep(NA,len_ts)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  d_ts_impute[["NR"]] = lapply(l_idx_d_NN, function(x) 
    if(length(x)>0){x$d_NR_}else{rep(NA,len_ts)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  d_ts_impute[["mF"]] = lapply(l_idx_corr, function(x) 
    if(length(x)>0){x$mF_}else{rep(NA,len_ts)}) %>% 
    bind_cols() # {len_idx} x {48} 
  
  # diff estacional ----
  d_ts_impute[["sttn"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["sttn"]]), 
                                               lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  
  d_ts_impute[["na.kalman"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["na.kalman"]]), 
                                                    lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  
  d_ts_impute[["CHIRPS_C"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["CHIRPS_C"]]), 
                                                   lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  
  d_ts_impute[["CCW"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["CCW"]]), 
                                              lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  
  d_ts_impute[["IDW"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["IDW"]]), 
                                              lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  
  d_ts_impute[["AA"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["AA"]]), 
                                             lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  
  d_ts_impute[["NR"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["NR"]]), 
                                             lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  
  d_ts_impute[["mF"]] = as.data.frame(mapply(FUN = pad_diff, x = as.data.frame(d_ts_impute[["mF"]]), 
                                             lag_ = lags_, SIMPLIFY = TRUE))[idx_,]
  ## métricas desempeño -----
  ## MAE 
  d_impute_MAE = list()
  ## NSE 
  d_impute_NSE = list()
  ## R_Spearman
  d_impute_R_s = list()
  
  d_impute_MAE[["na.kalman"]] = mae(sim = d_ts_impute$na.kalman,
                                    obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["IDW"]] = mae(sim = d_ts_impute$IDW,
                              obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["AA"]] = mae(sim = d_ts_impute$AA,
                             obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["NR"]] = mae(sim = d_ts_impute$NR,
                             obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["CCW"]] = mae(sim = d_ts_impute$CCW,
                              obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["mF"]] = mae(sim = d_ts_impute$mF,
                             obs = d_ts_impute$sttn, na.rm=TRUE)
  
  d_impute_MAE[["CHIRPS_C"]] = mae(sim = d_ts_impute$CHIRPS_C,
                                   obs = d_ts_impute$sttn, na.rm=TRUE)
  
  if(len_idx > 1){
    d_impute_NSE[["na.kalman"]] = NSE(sim = d_ts_impute$na.kalman,
                                      obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["IDW"]] = NSE(sim = d_ts_impute$IDW,
                                obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["AA"]] = NSE(sim = d_ts_impute$AA,
                               obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["NR"]] = NSE(sim = d_ts_impute$NR,
                               obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["CCW"]] = NSE(sim = d_ts_impute$CCW,
                                obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["mF"]] = NSE(sim = d_ts_impute$mF,
                               obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_NSE[["CHIRPS_C"]] = NSE(sim = d_ts_impute$CHIRPS_C,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
  }
  
  if(len_idx > 2){
    
    d_impute_R_s[["na.kalman"]] = rSpearman(sim = d_ts_impute$na.kalman,
                                            obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["IDW"]] = rSpearman(sim = d_ts_impute$IDW,
                                      obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["AA"]] = rSpearman(sim = d_ts_impute$AA,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["NR"]] = rSpearman(sim = d_ts_impute$NR,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["CCW"]] = rSpearman(sim = d_ts_impute$CCW,
                                      obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["mF"]] = rSpearman(sim = d_ts_impute$mF,
                                     obs = d_ts_impute$sttn, na.rm=TRUE)
    
    d_impute_R_s[["CHIRPS_C"]] = rSpearman(sim = d_ts_impute$CHIRPS_C,
                                           obs = d_ts_impute$sttn, na.rm=TRUE)
  }
  
  output_ = list()
  output_[["MAE"]] = lapply(d_impute_MAE, function(x) quantile(x, probs = 0.5, na.rm = TRUE) %>% as.vector())
  if(len_idx >1){
    output_[["r_s"]] = lapply(d_impute_R_s, function(x) quantile(x, probs = 0.5, na.rm = TRUE) %>% as.vector())
    output_[["NSE"]] = lapply(d_impute_NSE, function(x) quantile(x, probs = 0.5, na.rm = TRUE) %>% as.vector())
  }
  
  return(output_) 
}

## procesamiento en paralelo
plan(multisession, workers = parallel::detectCores() - 1)

#perf_missing = list()
perf_missing_full = list()
#perf_missing_81_10 = list()

perf_missing <- future_lapply(
  1:length(idx_d_missing2),
  function(i){
    sim_impute(m_distance, m_corr, d_ts_full, 
               idx_d_missing2[[i]], d_NR_full, d_CHIRPS, FALSE)
  }
)

perf_missing_full <- future_lapply(
  1:length(idx_missing2),
  function(i){
    sim_impute_full(m_distance, m_corr, ts_full, 
               idx_missing2[[i]], NR_full, ts_CHIRPS_full, lags_full, FALSE)
  }
)

perf_missing_81_10 <- future_lapply(
  1:length(idx_d_missing2_81_10),
  function(i){
    sim_impute(m_distance_81_10, m_corr_81_10, d_ts_full_81_10, 
               idx_d_missing2_81_10[[i]], d_NR_full_81_10, d_CHIRPS_81_10, 
               wdw_30 = TRUE)
  }
)

#saveRDS(perf_missing_full, file = 'performance_missing_full.RDS')
#perf_missing_full = readRDS('performance_missing_full.RDS')

table(sapply(idx_d_missing2_81_10, length))
idx_d_missing2_81_10[882]

length(perf_missing[[1]]) # 3 r_s, MAE, NSE
length(perf_missing[[1]][[1]]) #6 methods

df_MAE_sim = perf_missing[[1]][[1]] %>%  bind_cols()
df_R_Spearman_sim = perf_missing[[1]][[2]] %>%  bind_cols()
df_NSE_sim = perf_missing[[1]][[3]] %>%  bind_cols()

for(i in 2:length(perf_missing)){
  df_MAE_sim = bind_rows(df_MAE_sim, perf_missing[[i]][[1]] %>%  bind_cols())
  if(length(perf_missing[[i]])>1){
    df_R_Spearman_sim = bind_rows(df_R_Spearman_sim, perf_missing[[i]][[2]] %>%  bind_cols())
    df_NSE_sim = bind_rows(df_NSE_sim, perf_missing[[i]][[3]] %>%  bind_cols())
  }
}

df_MAE_sim_81_10 = perf_missing_81_10[[1]][[1]] %>%  bind_cols()
df_R_Spearman_sim_81_10 = perf_missing_81_10[[1]][[2]] %>%  bind_cols()
df_NSE_sim_81_10 = perf_missing_81_10[[1]][[3]] %>%  bind_cols()

for(i in 2:length(perf_missing_81_10)){
  df_MAE_sim_81_10 = bind_rows(df_MAE_sim_81_10, perf_missing_81_10[[i]][[1]] %>%  bind_cols())
  if(length(perf_missing_81_10[[i]])>1){
    df_R_Spearman_sim_81_10 = bind_rows(df_R_Spearman_sim_81_10, perf_missing_81_10[[i]][[2]] %>%  bind_cols())
    df_NSE_sim_81_10 = bind_rows(df_NSE_sim_81_10, perf_missing_81_10[[i]][[3]] %>%  bind_cols())
  }
}

df_MAE_sim_full = perf_missing_full[[1]][[1]] %>%  bind_cols()
df_R_Spearman_sim_full = perf_missing_full[[1]][[2]] %>%  bind_cols()
df_NSE_sim_full = perf_missing_full[[1]][[3]] %>%  bind_cols()

for(i in 2:length(perf_missing_full)){
  df_MAE_sim_full = bind_rows(df_MAE_sim_full, perf_missing_full[[i]][[1]] %>%  bind_cols())
  if(length(perf_missing_full[[i]])>1){
    df_R_Spearman_sim_full = bind_rows(df_R_Spearman_sim_full, perf_missing_full[[i]][[2]] %>%  bind_cols())
    df_NSE_sim_full = bind_rows(df_NSE_sim_full, perf_missing_full[[i]][[3]] %>%  bind_cols())
  }
}


which(df_R_Spearman_sim$CHIRPS_C >0.99)
lgc_len = sapply(idx_d_missing2, function(x) length(x)==2)
#lgc_len = sapply(idx_d_missing2, function(x) length(x)>=3 & length(x)<=4 )
idx_d_missing3 = idx_d_missing2[lgc_len]

idx_d_missing3[[52]] %>%  length()
idx_d_missing3[[415]] %>%  length()

# patrón de 4 datos faltantes
quantile(rSpearman(sim = d_CHIRPS[idx_d_missing3[[52]],], 
          obs = d_ts_full[idx_d_missing3[[52]], ], 
          na.rm = TRUE), probs = 0.5, na.rm = TRUE)

# patrón de 3 datos faltantes
quantile(rSpearman(sim = d_CHIRPS[idx_d_missing3[[415]],], 
                   obs = d_ts_full[idx_d_missing3[[415]], ], 
                   na.rm = TRUE), probs = 0.5, na.rm = TRUE)

d_impute_R_s_long = df_R_Spearman_sim %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "r_s") %>% 
  mutate(Method = str_c('d_',Method))

d_impute_R_s_long_81_10 = df_R_Spearman_sim_81_10 %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "r_s") %>% 
  mutate(Method = factor(str_c('d_',Method), 
                         levels = c("d_AA", "d_CCW", "d_CHIRPS_C", "d_IDW", 
                                    "d_mF", "d_na.kalman", "d_NR"),
                         ordered = TRUE))

d_impute_R_s_long_full = df_R_Spearman_sim_full %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "r_s") %>% 
  bind_rows(d_impute_R_s_long) %>% 
  mutate(Method = factor(Method, 
                         levels = c("AA", "CCW", "CHIRPS_C", "IDW", "mF", "na.kalman", "NR",
                                    "d_AA", "d_CCW", "d_CHIRPS_C", "d_IDW", "d_mF", "d_na.kalman", "d_NR"),
                         ordered = TRUE))

ggplot(d_impute_R_s_long_full, aes(x = Method, y = r_s, fill = Method)) +
#ggplot(d_impute_R_s_long_81_10, aes(x = Method, y = r_s, fill = Method)) +  
  geom_boxplot() +
  geom_hline(yintercept = 0.7, color = "red", linetype = "dashed", size = 1) +
  # título y nombres de ejes
  labs(title = "Boxplots del r_s", x = "Método imputación", y = "r_s") +
  theme_minimal() +
  scale_fill_manual(values = c(
    "AA" = '#fcd9cd', "CCW" = "#cfebf3", "CHIRPS_C" = "#d9ead3",
    "IDW" = "#fce5cd", "mF" = "#cfd9f3", "na.kalman" = '#f4cccc',
    "NR" = '#fcf1cd',
    "d_AA" = '#fcd9cd', "d_CCW" = "#cfebf3", "d_CHIRPS_C" = "#d9ead3",
    "d_IDW" = "#fce5cd", "d_mF" = "#cfd9f3", "d_na.kalman" = '#f4cccc',
    "d_NR" = '#fcf1cd')) +
  #coord_cartesian(ylim = c(0, max(bx_$df_outliers$upper_bound)))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24),
        plot.title = element_text(size=24, face="bold")) 


#MAE
d_impute_MAE_long = df_MAE_sim %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "MAE") %>% 
  mutate(Method = str_c('d_',Method))

d_impute_MAE_long_81_10 = df_MAE_sim_81_10 %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "MAE") %>% 
  mutate(Method = factor(str_c('d_',Method), 
                         levels = c("d_AA", "d_CCW", "d_CHIRPS_C", "d_IDW", 
                                    "d_mF", "d_na.kalman", "d_NR"),
                         ordered = TRUE))

d_impute_MAE_long_full = df_MAE_sim_full %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "MAE") %>% 
  bind_rows(d_impute_MAE_long) %>% 
  mutate(Method = factor(Method, 
                         levels = c("AA", "CCW", "CHIRPS_C", "IDW", "mF", "na.kalman", "NR",
                                    "d_AA", "d_CCW", "d_CHIRPS_C", "d_IDW", "d_mF", "d_na.kalman", "d_NR"),
                         ordered = TRUE))

ggplot(d_impute_MAE_long_full, aes(x = Method, y = MAE, fill = Method)) +
#ggplot(d_impute_MAE_long_81_10, aes(x = Method, y = MAE, fill = Method)) +
  geom_boxplot() +
  # título y nombres de ejes
  labs(title = "Boxplots del MAE", x = "Método imputación", y = "MAE") +
  theme_minimal() +
  scale_fill_manual(values = c(
    "AA" = '#fcd9cd', "CCW" = "#cfebf3", "CHIRPS_C" = "#d9ead3",
    "IDW" = "#fce5cd", "mF" = "#cfd9f3", "na.kalman" = '#f4cccc',
    "NR" = '#fcf1cd',
    "d_AA" = '#fcd9cd', "d_CCW" = "#cfebf3", "d_CHIRPS_C" = "#d9ead3",
    "d_IDW" = "#fce5cd", "d_mF" = "#cfd9f3", "d_na.kalman" = '#f4cccc',
    "d_NR" = '#fcf1cd')) +
  #coord_cartesian(ylim = c(0, max(bx_$df_outliers$upper_bound)))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24),
        plot.title = element_text(size=24, face="bold")) 

#NSE
d_impute_NSE_long = df_NSE_sim %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "NSE") %>% 
  mutate(Method = str_c('d_',Method))

d_impute_NSE_long_81_10 = df_NSE_sim_81_10 %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "NSE") %>% 
  mutate(Method = factor(str_c('d_',Method), 
                         levels = c("d_AA", "d_CCW", "d_CHIRPS_C", "d_IDW", 
                                    "d_mF", "d_na.kalman", "d_NR"),
                         ordered = TRUE))

d_impute_NSE_long_full = df_NSE_sim_full %>% 
  pivot_longer(cols = everything(), names_to = "Method", values_to = "NSE") %>% 
  bind_rows(d_impute_NSE_long) %>% 
  mutate(Method = factor(Method, 
                         levels = c("AA", "CCW", "CHIRPS_C", "IDW", "mF", "na.kalman", "NR",
                                    "d_AA", "d_CCW", "d_CHIRPS_C", "d_IDW", "d_mF", "d_na.kalman", "d_NR"),
                         ordered = TRUE))

ggplot(d_impute_NSE_long_full, aes(x = Method, y = NSE, fill = Method)) +
#ggplot(d_impute_NSE_long_81_10, aes(x = Method, y = NSE, fill = Method)) +
  geom_boxplot() +
  # título y nombres de ejes
  labs(title = "Boxplots del NSE", x = "Método imputación", y = "NSE") +
  theme_minimal() +
  scale_fill_manual(values = c(
    "AA" = '#fcd9cd', "CCW" = "#cfebf3", "CHIRPS_C" = "#d9ead3",
    "IDW" = "#fce5cd", "mF" = "#cfd9f3", "na.kalman" = '#f4cccc',
    "NR" = '#fcf1cd',
    "d_AA" = '#fcd9cd', "d_CCW" = "#cfebf3", "d_CHIRPS_C" = "#d9ead3",
    "d_IDW" = "#fce5cd", "d_mF" = "#cfd9f3", "d_na.kalman" = '#f4cccc',
    "d_NR" = '#fcf1cd')) +
  coord_cartesian(ylim = c(-2, 1))+
  # sin leyenda y tamaño del texto de los ejes
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24),
        plot.title = element_text(size=24, face="bold"))

# Check point 3 ----
#save.image('df_eval_metrics_missing.RData')
load('df_eval_metrics_missing.RData')

ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) +
  #geom_sf(data = filter(sttns.points, DEPARTAMENTO != ADZ), color = 'blue', size = 1) +
  geom_sf(data = filter(sttns.points, n ==480), color = 'red', size = 1) +
  geom_sf(data = sttns_sinNN2, color = 'blue', size = 1) +
  #theme_minimal() +  # Tema minimalista
  labs(title = "Distribución espacial de sttns sin datos faltantes")   


# df_missing: data frame con 30 años continuos ----
df_missing <- sttns.points %>%
  #filter(region =='Andina', n < 480) %>%
  #filter(region =='Caribe', n < 480) %>%
  #filter(region =='Pacífica', n < 480) %>%
  #filter(region =='Orinoquía', n < 480) %>%
  #filter(region =='Orinoquía', n >= 470) %>%
  #filter(region =='Amazonía', n < 480) %>%
  #filter(region =='Amazonía', n >= 467) %>% 
  dplyr::select(matches("^(19|20)"))
  #dplyr::select(matches("^(199|20)")) %>% dplyr::select(-matches("^1990"))
  #1986 - 2015 
  #dplyr::select(matches("^(19|20||CodigoEstacion|NombreEstacion)")) %>% dplyr::select(-matches("^(1981|1982|1983|1984|1985|2016|2017|2018|2019|2020)")) %>% 
  # 1991 - 2020 
  #dplyr::select(matches("^(199|20|CodigoEstacion|NombreEstacion)")) %>% dplyr::select(-matches("^1990")) %>% 
  # 1981 - 2010 
  #dplyr::select(matches("^(19|200|2010|CodigoEstacion|NombreEstacion)")) %>% 
  #pivot_longer(cols = matches("^(19|20)"),
  #pivot_longer(cols = matches("^(19|20)"),
  #pivot_longer(cols = matches("^(19|200|2010)"),
  #             names_to = 'mes', values_to = 'pcp') %>% 
  #filter(!is.na(pcp)) %>% 
  #count(CodigoEstacion, NombreEstacion)

cod_sttn_86_15 = filter(df_missing, n == 360) %>% distinct(CodigoEstacion) %>% mutate(yy86_15 = 1)
cod_sttn_91_20 = filter(df_missing, n == 360) %>% distinct(CodigoEstacion) %>% mutate(yy91_20 = 1)
cod_sttn_81_10 = filter(df_missing, n == 360) %>% distinct(CodigoEstacion) %>% mutate(yy81_10 = 1)

df_completo = full_join(cod_sttn_81_10, cod_sttn_91_20, by = 'CodigoEstacion') %>% 
  full_join(cod_sttn_86_15, by = 'CodigoEstacion') %>% 
  mutate(across(where(is.numeric), ~ replace_na(., 0)))

df_completo %>% mutate(s_ = yy81_10+yy91_20+yy86_15) %>% count(s_)

#df_completo %>% mutate(s_ = yy81_10+yy91_20+yy86_15) %>% filter(s_ ==1) %>%
df_completo %>% mutate(s_ = yy81_10+yy91_20+yy86_15) %>% filter(s_ ==2) %>% 
  dplyr::select('CodigoEstacion', 'yy81_10', 'yy86_15', 'yy91_20') %>% 
  pivot_longer(cols = c('yy81_10', 'yy91_20', 'yy86_15'),
               names_to = "yy", values_to = 'value') %>% filter(value ==1) %>%  
               #names_to = "yy", values_to = 'value') %>% filter(value ==1) %>% 
  group_by(CodigoEstacion) %>%  summarise(yy = paste0(min(yy),'-',max(yy)) ) %>%
  count(yy)

cod_sttn_30 = filter(df_missing, n == 360) %>% distinct(CodigoEstacion)
sttns_30 = filter(sttns.points, CodigoEstacion %in% cod_sttn_30$CodigoEstacion)

sttns_30 = filter(sttns.points, CodigoEstacion %in% df_completo$CodigoEstacion)
sttns_30 %>% count(region)

idx_missing <- apply(df_missing, 1, function(row) which(is.na(row)))
v_idx = c()
#v_idx_Andina = c()
#v_idx_Caribe = c()
#v_idx_Pacifica = c()
#v_idx_Orinoquia = c()
#v_idx_Amazonia = c()
v_idx_30 = c()

for(i in 1:length(idx_missing)){
  v_idx = c(v_idx, as.vector(idx_missing[[i]]))
  #v_idx_Andina = c(v_idx_Andina, as.vector(idx_missing[[i]]))
  #v_idx_Caribe = c(v_idx_Caribe, as.vector(idx_missing[[i]]))
  #v_idx_Pacifica = c(v_idx_Pacifica, as.vector(idx_missing[[i]]))
  #v_idx_Orinoquia = c(v_idx_Orinoquia, as.vector(idx_missing[[i]]))
  #v_idx_Amazonia = c(v_idx_Amazonia, as.vector(idx_missing[[i]]))
  #v_idx_30[i] = length(idx_missing[[i]])==0  
}

sum(v_idx_30) # 68 (1991-2020) #90 (1981-2010) # 80 (1986-2015)
# intersectando los tres conjuntos

plot(density(v_idx_Amazonia), lty = 1, lwd = 2,col = '#91B495')
lines(density(v_idx_Orinoquia), lty = 1, lwd = 2,col = '#B5E794')
lines(density(v_idx_Pacifica), lty = 1, lwd = 2,col = '#A59CB3')
lines(density(v_idx_Caribe), lty = 1, lwd = 2, col = '#F8E595')
lines(density(v_idx_Andina), lty = 1, lwd = 2, col = '#D4B698')
lines(density(v_idx), lty = 1, lwd = 2, col = 'black')

scale_fill_manual(values = c("Amazonía" = "#91B495", "Andina" = "#D4B698", "Caribe" = "#F8E595",
                             "Orinoquía" = "#B5E794", "Pacífica" = "#A59CB3")) +
  

# A1. Compara ALTITUD_IDEAM con SRTM30 ----
SRTM_30.sttns = cbind(SRTM_30.sttns, 'ALTITUD' = sttns$ALTITUD,
                      'DEPARTAMENTO' = sttns$DEPARTAMENTO,
                      'CodigoEstacion' = sttns$CodigoEstacion) 

SRTM_30.sttns = SRTM_30.sttns %>%  mutate(dif_ = ALTITUD - SRTM_30_Col1)

## Histograma del sesgo: Altitud_IDEAM - SRTM30
ggplot(SRTM_30.sttns, aes(x = dif_)) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 breaks = seq(-3000,2000, by =100)) +
  stat_bin(breaks = seq(-3000,2000, by =100),
           geom = "text",
           aes(label = after_stat(count)), # usa la frecuencia como etiqueta
           vjust = -0.5, # posición vertical del texto
           size =10 # aumenta el tamaño de las etiquetas
  )

quantile(SRTM_30.sttns$dif_, probs = seq(0,1, by =0.1))
## criterio1: abs(sesgo) > 60 (valor de los P10 y P90)
## criterio2: abs(PBIAS) > 25% (desempeño regular)
ALTTD_ = SRTM_30.sttns %>%
  mutate(PBIAS = 100*((SRTM_30_Col1- ALTITUD) / ALTITUD)) %>% 
  filter(abs(dif_)>60 & abs(PBIAS)>25)

ALTTD_2 = sttns.points[ALTTD_$ID,]
# join para traer PBIAS y graficar los puntos con PBIAS alto
ALTTD_2 = inner_join(ALTTD_2, select(ALTTD_,CodigoEstacion, PBIAS), by ='CodigoEstacion')

shp_depto <- st_read("Departamentos202208_shp/Depto.shp")
shp_regions <- st_read("Regiones_shp/regiones.shp")

## chequeando las diferencias altas para el DPTO de Santander
#SNTDR = filter(SRTM_30.sttns, abs(dif_)>500 & DEPARTAMENTO=='Santander')
#SNTDR = sttns.points[SNTDR$ID,]

## grafico mapa DPTO - regiones, puntos sttns con diferencias
ggplot() + 
  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1.5) +
  geom_sf(data = shp_depto, alpha =0.2) +
  # gráfica las sttns con sobreestimación (ALTTD_IDEAM << SRTM30)
  geom_sf(data = subset(ALTTD_2,PBIAS>0 ), col = 'red', lwd = 0.2)+
  # gráfica las sttns con subestimación (ALTTD_IDEAM >> SRTM30) 
  geom_sf(data = subset(ALTTD_2,PBIAS<0 ), col = 'blue', lwd = 0.2) 

## SRTM_30.sttns_: df para graficas histogramas
## pivota cols a filas 
SRTM_30.sttns_ = SRTM_30.sttns %>% rename('SRTM_30'=SRTM_30_Col1, 'Altitud_IDEAM'=ALTITUD) %>% 
  pivot_longer(cols = c("Altitud_IDEAM", "SRTM_30"), 
               names_to = 'Serie', values_to = 'Altitud')

## histogramas de altitud
ggplot(SRTM_30.sttns_, aes(x =Altitud , fill=Serie)) + 
  geom_histogram(alpha = 0.3, position = "identity", boundary = 0, 
                 breaks = seq(0,4500, by =100)) +
  scale_x_continuous(breaks = seq(0, 4500, by = 400)) +
  scale_fill_manual(values = c("Altitud_IDEAM" = "orange", "SRTM_30" = "blue"))+
  labs(title = "Histogramas de altitud de las estaciones metereológicas", 
       x = "Altitud (m)", y = "Count", fill = "Altitud") +
  theme(axis.text=element_text(size=28),
        axis.text.x = element_text(angle =90, hjust = 1),
        plot.title = element_text(size=28,face="bold"),
        axis.title=element_text(size=28),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28)
  )

# A2: evalua estacionalidad a ts CHRPS (No R. H0 - concluye no-estacionaria) ----
## A2.1 obtiene orden AR usando ts CHIRPS 
## con diferencia ordinaria (diff_ = TRUE)
m_order = get_ar_order(x = pol2_list, diff_ = TRUE)

# los NAs de las ts sttns ocasionan error a ajustar el modelo AR (por eso usa na.approx)
#na_count_ = c()
#for(i in 1:length(pol2_list)){
#  na_count_[i] = sum(is.na(pol2_list[[i]]$sttns))
#}
# !!! para ejecutar m_order_sttns, modificar (descomentar/comentar) get_ar_order
# m_order_sttns = get_ar_order(x = pol2_list, len = length(pol2_list), diff_ = TRUE)

table(m_order)

## sin diferencia ordinaria (diff_ = FALSE)
m_order2 = get_ar_order(x = pol2_list, len = length(pol2_list), diff_ = FALSE)
table(m_order2)

## A2.2 gráficos de las 3 estaciones que arrojan error al ajustar un modelo AR con diff ordinaria
no_order = which(m_order==-1)
ID_no_order = c()
k = 1
for(i in no_order){
  ID_no_order[k] = unique(pol2_list[[i]]$ID)
  k <- k+1
}

#plot(pol3_list[[7]]$chirps.v2.0, type = 'l', main = unique(pol3_list[[7]]$ID))
ggplot() + 
  geom_sf(data = shp_depto) +
  geom_sf(data = sttns.points[sttns.points$CodigoEstacion %in% ID_no_order,], 
          col = 'red', lwd = 0.2) +
  geom_sf_text(aes(label = CodigoEstacion), 
               data = sttns.points[sttns.points$CodigoEstacion %in% ID_no_order,], 
               size = 3) 

## A2.3 prueba aTSA::adf.test a ts CHIRPS
## con m_order (arg lag prueba ADF) orden AR obtenido diff ordinaria
## type1: no drift and no trend
df_ADF = ADF_test(pol2_list, m_order, type = 'type1')

#df_ADF = ADF_test(pol2_list, m_order, type = 'type1') %>% 
#  # type2: drift but no trend
#  bind_cols(ADF_test(pol2_list, m_order, type = 'type2')) %>% 
#  # type3: drift and linear trend 
#  bind_cols(ADF_test(pol2_list, m_order, type = 'type3'))

table(df_ADF$RD_t1) # No R. H0: concluye Raíz unitaria (no-estacionaria)

## con m_order2 (arg lag prueba ADF) orden AR obtenido sin diff ordinalmente
df_ADF2 = ADF_test(pol2_list, m_order2, type = 'type1')

table(df_ADF2$RD_t1) # No R. H0: concluye Raíz unitaria (no-estacionaria)

## A2.4 prueba urca::ur.df a ts CHIRPS

df_ur = data.frame(RD_t1 = rep('', length(pol2_list)))
m_order_ = ifelse(m_order<0, 12, m_order)

for(i in 1:length(pol2_list)){
  ADF_ur = urca::ur.df(pol2_list[[i]]$`chirps-v2-0`, type  = "none", 
                       lags = m_order_[i], selectlags = "AIC")
  #lags = m_order2[i], selectlags = "AIC")
  df_ur$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur$RD_t1) # No R. H0: concluye Raíz unitaria (no-estacionaria)
# test urca::ur.df con m_order2 el arg orden AR obtenido sin diff ordinalmente
#table(df_ur$RD_t1) # No R. H0: concluye Raíz unitaria (no-estacionaria)

# A3 Descomposición STL usando librería base (No R. H0 - concluye no-estacionaria) ----
# t.window = 21 (fpp3) auto
#nextodd <- function(x) {
#  x <- round(x)
#  if (x%%2 == 0) 
#    x <- x + 1
#  as.integer(x)}

#nextodd(ceiling(1.5 * frequency(y)/(1 - 1.5/(10 * as.integer(length(y)) + 1))))

stl_chirps.v2.0 = list()
for(i in 1:len_){
  y <- tsibble::tsibble(pol2_list[[i]] %>% 
                          mutate(Date = as.Date(Date, format = '%Y-%m-%d')) 
                        %>% dplyr::select(c("Date" ,  "chirps-v2-0")), index = Date)
  
  y = y %>% mutate(Date = yearmonth(Date)) %>% 
    as_tsibble(index = Date)
  # t.window = 21 (fpp3) auto
  #nextodd(ceiling(1.5 * frequency(y)/(1 - 1.5/(10 * as.integer(length(y)) + 1))))
  fit1 <- stl(y, s.window = "periodic")
  # t.window = 13 (fpp2)
  #fit11<- stl(y, t.window = 13, s.window = "periodic")
  # Selecting a shorter trend window, avoid produces a trend-cycle component that is too rigid
  #fit12<- stl(y, t.window = 7, s.window = "periodic")
  stl_chirps.v2.0[[i]] =  data.frame(chirps_col = y$`chirps-v2-0` - fit1$time.series[, "seasonal"]) %>% 
    rename(`chirps-v2-0` = chirps_col)
  #stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit11$time.series[, "seasonal"]
  #stl_chirps.v2.0[[i]] = y$chirps.v2.0 - fit12$time.series[, "seasonal"]
}

## A3.1 obtiene orden AR usando ts CHIRPS 
## con diferencia ordinaria (diff_ = TRUE)
m_order5 = get_ar_order(x = stl_chirps.v2.0, diff_ = TRUE)
table(m_order5)

## sin diferencia ordinaria (diff_ = TRUE)
m_order6 = get_ar_order(x = stl_chirps.v2.0, diff_ = FALSE)
table(m_order6)
m_order6 = ifelse(m_order6==0, 12, m_order6)

## A3.2 prueba aTSA::adf.test a ts stl_CHIRPS
## con m_order5 (arg lag prueba ADF) orden AR obtenido diff ordinaria
## type1: no drift and no trend
df_ADF5 = ADF_test(stl_chirps.v2.0, m_order5, type = 'type1')
table(df_ADF5$RD_t1) # 1004 No R. H0: concluye hay Raíz unitaria (no-estacionaria)

## con m_order6 (arg lag prueba ADF) orden AR obtenido sin diff ordinalmente
df_ADF6 = ADF_test(stl_chirps.v2.0, m_order6, type = 'type1')
table(df_ADF6$RD_t1) # 604 No R. H0: concluye hay Raíz unitaria (no-estacionaria)

## A3.3 prueba urca::ur.df a ts stl_CHIRPS
df_ur5 = data.frame(RD_t1 = rep('', len_))
for(i in 1:len_){
  ADF_ur = urca::ur.df(stl_chirps.v2.0[[i]]$`chirps-v2-0`, type  = "none", 
                       lags = m_order5[i], selectlags = "AIC")
  #lags = m_order6[i], selectlags = "AIC")
  df_ur5$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur5$RD_t1) # 1004 No R. H0: concluye hay Raíz unitaria (no-estacionaria)
# test urca::ur.df con m_order6 el arg orden AR obtenido sin diff ordinalmente
#table(df_ur5$RD_t1) # 856 No R. H0: concluye hay Raíz unitaria (no-estacionaria)

# A4 Descomposición STL usando librería feast (No R. H0 - concluye no-estacionaria) ----
f_stl_chirps.v2.0 = list()
for(i in 1:len_){
  y <- tsibble::tsibble(pol2_list[[i]] %>% 
                          mutate(Date = as.Date(Date, format = '%Y-%m-%d')) 
                        %>% dplyr::select(c("Date" ,  "chirps-v2-0")), index = Date)
  
  y = y %>% mutate(Date = yearmonth(Date)) %>% 
    as_tsibble(index = Date)
  
  # The default setting for monthly data is trend(window=21)
  # t.window = 13 (fpp2)
  # Selecting a shorter trend window, avoid produces a trend-cycle component that is too rigid
  fit2 = y |>
    model(
      #STL(`chirps.v2.0` ~ trend(window = 13) +
      STL(`chirps-v2-0` ~ 
            season(window = "periodic"),
          robust = TRUE)) |>
    components() 
  f_stl_chirps.v2.0[[i]] = data.frame(chirps_col = y$`chirps-v2-0` - fit2$season_year) %>% 
    rename(`chirps-v2-0` = chirps_col)
}

## A4.1 obtiene orden AR usando ts f_stl_CHIRPS 
## con diferencia ordinaria (diff_ = TRUE)
m_order7 = get_ar_order(x = f_stl_chirps.v2.0, diff_ = TRUE)
table(m_order7)

## sin diferencia ordinaria (diff_ = TRUE)
m_order8 = get_ar_order(x = f_stl_chirps.v2.0, diff_ = FALSE)
table(m_order8)
m_order8 = ifelse(m_order8==0, 12, m_order8)

## A4.2 prueba aTSA::adf.test a ts f_stl_CHIRPS
## con m_order7 (arg lag prueba ADF) orden AR obtenido diff ordinaria
## type1: no drift and no trend
df_ADF7 = ADF_test(f_stl_chirps.v2.0, m_order7, type = 'type1')
table(df_ADF7$RD_t1) # 1004 No R. H0: concluye hay Raíz unitaria (no-estacionaria)

## con m_order8 (arg lag prueba ADF) orden AR obtenido sin diff ordinalmente
df_ADF8 = ADF_test(f_stl_chirps.v2.0, m_order8, type = 'type1')
table(df_ADF8$RD_t1) # 596 No R. H0: concluye hay Raíz unitaria (no-estacionaria)

## A4.3 prueba urca::ur.df a ts f_stl_CHIRPS
df_ur9 = data.frame(RD_t1 = rep('', len_))
m_order7 = ifelse(m_order7==-1, 12, m_order7)
for(i in 1:len_){
  ADF_ur = urca::ur.df(f_stl_chirps.v2.0[[i]]$`chirps-v2-0`, type  = "none", 
                       lags = m_order7[i], selectlags = "AIC")
  #lags = m_order8[i], selectlags = "AIC")
  df_ur9$RD_t1[i] = interp_urdf_FM(ADF_ur)
}

table(df_ur9$RD_t1) # 1004 No R. H0: concluye hay Raíz unitaria (no-estacionaria)
# test urca::ur.df con m_order8 el arg orden AR obtenido sin diff ordinalmente
#table(df_ur9$RD_t1) # 841 No R. H0: concluye hay Raíz unitaria (no-estacionaria)

## A4.4 grafica compara ts desestacionalizadas
## diff estacional (rojo) y removiendo comp. estacional con LOESS (azul, gris)
set.seed(24112024)
s = sample(1:len_, size = 1)

min_length <- min(length(stl_chirps.v2.0[[s]]$`chirps-v2-0`), length(d_chirps.v2.0[[s]]$`chirps-v2-0`))
s1_adj <- window(stl_chirps.v2.0[[s]]$`chirps-v2-0`, end = time(stl_chirps.v2.0[[s]]$`chirps-v2-0`)[min_length])
y_adj = window(pol2_list[[s]]$`chirps-v2-0`, end = time(pol2_list[[s]]$`chirps-v2-0`)[min_length])
s2_adj <- d_chirps.v2.0[[s]]$`chirps-v2-0`
s3_adj <- window(f_stl_chirps.v2.0[[s]]$`chirps-v2-0`, end = time(f_stl_chirps.v2.0[[s]]$`chirps-v2-0`)[min_length])

plot(s2_adj, type = "l", col = "red", lwd = 2, lty = 2, xlab = "Index", ylim  = c(-200,700))
lines(1:length(s2_adj), s1_adj, type = "l", col = "blue", lwd = 2, lty =2)
lines(1:length(s2_adj), y_adj, type = "l", col = "black", lwd = 2)
lines(1:length(s2_adj), s3_adj, type = "l", col = "gray", lwd = 2, lty =2)

legend("topleft", legend = c("serie original CHIRPS", "sin comp. estacional - stats::stl", 
                             "sin comp. estacional - feast::STL", "dif estacional"), 
       col = c("black", "blue", "gray", "red"), lwd = 2, lty = c(1,2,2,2), bty = "n")

# A5. Gráfico coef corr spearman cruzado ---------------------------------------
ccf_spearman <- function(x, y, max_lag = 10) {
  # Argumentos: ----------------------------------------------------------------
  # 
  # x df con ts num 
  # y df con ts num con misma dim que x
  # max_lag num rezago máximos
  #----------------------------------------------------------------------------
  x_rank <- x
  y_rank <- y
  
  lags <- seq(-max_lag, max_lag)
  cor_values <- sapply(lags, function(lag) {
    if (lag < 0) {
      # (1 - lag) adelanta la serie y_rank, corr(x_t, y_{t+lag}) para t = 1,...,(n_x - lag)
      cor(x_rank[1:(length(x) + lag)], y_rank[(1 - lag):length(y)], 
          method = "spearman", use = "complete.obs")
    } else if (lag > 0) {
      # (1 + lag): adelanta la serie x_rank, corr(x_{t+lag}, y_{t}) para t = 1,...,(n_y - lag)
      cor(x_rank[(1 + lag):length(x)], y_rank[1:(length(y) - lag)], 
          method = "spearman", use = "complete.obs")
    } else {
      # lag = 0 corr instantánea corr(x_t, y_t) para t = 1,..,n
      cor(x_rank, y_rank, method = "spearman", use = "complete.obs")
    }
  })
  df_cor <- as.data.frame(t(cor_values))
  colnames(df_cor) <- paste0("Lag_", lags)
  
  return(df_cor)
}

d_ccf_R_Spearman = list()
for(i in 1:len_){
  d_ccf_R_Spearman[[i]] = ccf_spearman(d_chirps.v2.0[[i]]$`chirps-v2-0`, 
                                       d_sttns[[i]]$sttns, max_lag = 12)
}

df_d_ccf_R_Spearman = d_ccf_R_Spearman %>% bind_rows()
#df_d_ccf_R_Spearman = df_d_ccf_R_Spearman %>% bind_cols(R_s = d_R_Spearman)
#table(df_d_ccf_R_Spearman$Lag_0 == df_d_ccf_R_Spearman$R_s)

# Pivotear el data frame a formato largo
df_d_ccf_R_Spearman <- df_d_ccf_R_Spearman %>%
  pivot_longer(cols = everything(), 
               names_to = "Lag", 
               values_to = "Spearman_Correlation")

df_d_ccf_R_Spearman = df_d_ccf_R_Spearman %>% 
  mutate(Lag = factor(Lag, levels = paste0("Lag_", -12:12)))

# Crear el gráfico de boxplot ccf R Spearman 
ggplot(df_d_ccf_R_Spearman, aes(x = Lag, y = Spearman_Correlation)) +
  geom_boxplot(fill = "lightblue", color = "black") + 
  theme_minimal() +
  labs(title = "Distribución de la Correlación de Spearman por Lag",
       x = "Lag",
       y = "Coeficiente de Spearman") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotar etiquetas en eje X
