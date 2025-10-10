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
library(qmap)

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

# 2. Carga datos elevación de datos abiertos ----
SRTM_30 <- rast("Servicio-159/SRTM30/SRTM_30_Col1.tif")
print(SRTM_30)

SRTM_30 <- project(SRTM_30, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
crs(SRTM_30) 

# 3. Convierte a objeto sf ----
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
# 4. Extrae la altitud de los puntos de las estaciones ----
SRTM_30.sttns <- terra::extract(SRTM_30, sttns.points)
v_ = which(is.na(SRTM_30.sttns$SRTM_30_Col1)) # filas sin elevación

for (i in v_){
  SRTM_30.sttns$SRTM_30_Col1[i] = terra::extract(geodata::elevation_3s(
    lon = pcp_col$LONGITUD[i], lat = pcp_col$LATITUD[i], path = './tmpr'), 
    cbind(pcp_col$LONGITUD[i], pcp_col$LATITUD[i])
    )[[1]]
}

# 5. Función para extraer los pixeles chirps más cercanos ----
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
                      %>%  select(c("Date" ,  "chirps-v2-0")), index = Date)

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

View(d_chirps.v2.0[[i]])

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


