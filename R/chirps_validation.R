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
shp_regions <- st_read("Regiones_shp/regiones.shp")
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
sf::sf_use_s2(FALSE)
sttns.points = st_join(sttns.points, shp_regions)

ggplot() + 
  geom_sf(data = shp_regions, col = 'darkgreen', alpha =0.2, lwd =1) + 
  #geom_sf(data = filter(sttns.points, Subregion=='Andes'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Caribe'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Pacifico'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Llanos'), color = "blue", size = 3)
  #geom_sf(data = filter(sttns.points, Subregion=='Amazonas'), color = "blue", size = 3)
  geom_sf(data = filter(sttns.points, is.na(Subregion)), color = "blue", size = 3)

## asigna la Subregión manualmente a las sttns que NO cruzan en el join espacial
sttns.points = sttns.points %>% 
  mutate(Subregion = ifelse(is.na(Subregion) & DEPARTAMENTO == 'Arauca', 'Llanos', Subregion),
         Subregion = ifelse(is.na(Subregion) & DEPARTAMENTO %in% c(ADZ, 'Sucre'), 'Caribe', Subregion),
         Subregion = ifelse(is.na(Subregion) & DEPARTAMENTO == 'Choco', 'Pacifico', Subregion),
         Subregion = ifelse(is.na(Subregion) & DEPARTAMENTO == 'Amazonas', 'Amazonas', Subregion))

prop.table(table(sttns.points$Subregion))

sttns.points = sttns.points %>% 
  mutate(region = factor(ifelse(Subregion == 'Amazonas', 'Amazonía',
                                ifelse(Subregion == 'Andes', 'Andina',
                                       ifelse(Subregion == 'Pacifico', 'Pacífica',
                                              ifelse(Subregion == 'Llanos', 'Orinoquía',
                                                     Subregion ))))))


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
f_MAE = function(df, obs, mod){
  MAE = df %>%
    mutate(abs_error = abs({{mod}} - {{obs}})) %>%
    summarise(MAE = mean(abs_error, na.rm = TRUE)) %>%
    pull(MAE)
  return(MAE)
}

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

# NSE (Nash-Sutcliffe Efficiency) ---------------------------------------------
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

quantile(NSE_)
mean(NSE_) # 0.33
# indica que, Chirps tiene una buena predicción, mejor que usar 
## la media de los valores observados como predicción, sin embargo
### el promedio está lejos de 1, aún hay un margen significativo para mejorar
hist(NSE_)

# Ensamble Bias Correction ---------------------------------------------------
