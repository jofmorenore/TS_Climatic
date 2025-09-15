library(readr)
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(lubridate)   # para trabajar con fechas
library(scales)      # para formatear los breaks en el eje X

# 1. Crea df_sttns ----
setwd('~/TS_climatic/IDEAM')
files = list.files(pattern = '.csv')

## carga los datos de las estaciones del IDEAM
lista = list()
for(i in files){
  lista[[i]] = read_csv(i)
}

## concatena por filas
### la PK es la combinación CodigoEstacion y Fecha
df_sttns = lista %>% bind_rows()
df_sttns$Fecha = as.Date(df_sttns$Fecha)

#nrow(filter(df_sttns, Fecha < as.Date("2021-01-01")))

### Num de estaciones: 1476
tbl_sttns =  filter(df_sttns, Fecha < as.Date("2021-01-01")) %>%
  count(CodigoEstacion, NombreEstacion)

nrow(tbl_sttns)
###     0%    25%    50%    75%   100% 
###   7.00 411.75 452.00 471.00 480.00 
quantile(tbl_sttns$n)
### Completitud >90%: 938 estaciones
c90 = round((40*12)*0.9)
c65 = round((40*12)*0.65)
nrow(filter(tbl_sttns, n>c90))

table(df_sttns$NivelAprobacion)

# 2. Agrega columnas de CNE_IDEAM a tbl_sttns ----

CNE_IDEAM <- read_excel("CNE_IDEAM.xls")

col_IDEAM = c("CODIGO", "CATEGORIA", "TECNOLOGIA", "ESTADO", "FECHA_INSTALACION", 
              "FECHA_SUSPENSION", "ALTITUD", "LATITUD", "LONGITUD", "DEPARTAMENTO", 
              "MUNICIPIO", "OBSERVACION")

### Agrega el DPTO y las otras cols
tbl_sttns2 = CNE_IDEAM %>% select(all_of(col_IDEAM)) %>% 
  inner_join(tbl_sttns, by = c('CODIGO'='CodigoEstacion'))

# 3. Criterios de selección de estaciones ----

## C1. Completitud 
### >90% de los datos (Urrea et al., 2016)
### min de 25 datos en las series de los meses calendario (Pedraza et al. 2018)
### OMM (2011)

### DPTOS de la región oriente con poca densidad de sttns (ver A2)
### c65 (completitud 65%) debido a (publicación del IDEAM) donde menciona "hasta el 65%"
### "Amazonas" "Arauca"   "Caqueta"  "Casanare" "Guainía"  "Guaviare" "Vaupes" "Vichada"
DPTO_65 = c("Amazonas","Arauca","Caqueta","Casanare","Guainía","Guaviare","Vaupes","Vichada")

### Filtra las estaciones que complen criterio de completitud 
tbl_sttns3 = filter(tbl_sttns2, ifelse(DEPARTAMENTO %in% DPTO_65, n>c65, n>c90))

sort(prop.table(table(tbl_sttns3$CATEGORIA)))
sort(prop.table(table(tbl_sttns3$TECNOLOGIA)))
sort(prop.table(table(tbl_sttns3$ESTADO)))

df_sttns2 = inner_join(df_sttns, rename(tbl_sttns3, "CodigoEstacion"="CODIGO"), 
                       by = c("CodigoEstacion", "NombreEstacion"))

# carga df con SUPERFICIE y REGION por DPTO
#KM2 = read_excel("Superficie_DPTO.xlsx")
# join para traer la columna REGION 
#df_sttns2 = inner_join(df_sttns2, select(KM2, -SUPERFICIE_KM2), by = 'DEPARTAMENTO')

# pivotar para que cada fila sea una sttn
df_sttns2_train = df_sttns2 %>% arrange(Fecha) %>% filter(Fecha < as.Date("2021-01-01")) %>% 
  pivot_wider(
    id_cols = c(CodigoEstacion,NombreEstacion,DEPARTAMENTO,n,ALTITUD,LATITUD,LONGITUD),
    names_from = Fecha, 
    values_from = Valor
  )

arrow::write_parquet(df_sttns2_train, "sttns_pcp_col_train.parquet")

df_sttns2_test = df_sttns2 %>% arrange(Fecha) %>% filter(Fecha >= as.Date("2021-01-01")) %>% 
  pivot_wider(
    id_cols = c(CodigoEstacion, NombreEstacion,DEPARTAMENTO,REGION,ALTITUD, LATITUD, LONGITUD),
    names_from = Fecha, 
    values_from = Valor
  )

arrow::write_parquet(df_sttns2_test, "sttns_pcp_col_test.parquet")

#sttns_df = arrow::read_parquet("sttns_pcp_col_train.parquet")

# A1. Histograma datos faltantes 1981-01 a 1985-12 ----
## df_sttns2_: df para graficar histograma del recuento de datos faltantes en 1981-1985
## filtro para restringir al periodo 1981 - 1985
## pivota las cols de meses (para que queden los NAs en las filas)
## crea la col gap, Fecha y sttn 

df_sttns2_ = filter(df_sttns2, Fecha < as.Date("1986-01-01")) %>% 
  pivot_wider(
    id_cols = c(CodigoEstacion,NombreEstacion,DEPARTAMENTO,n,ALTITUD,LATITUD,LONGITUD),
    names_from = Fecha, 
    values_from = Valor
  )%>% 
  pivot_longer(cols = matches("^(19|20)"), # columnas que comienzan por 19.. o 20..
               names_to  = "Fecha",
               values_to = "pcp") %>% 
  mutate(gap = is.na(pcp),
         Fecha = as.Date(Fecha),
         sttn = paste0(DEPARTAMENTO,' - ',NombreEstacion, ' - ', round((n/480)*100), "%")) %>% 
  arrange(desc(n))

## agrupa por cols, y resume en t_gap la suma del número de NAs
df_sttns2_gr = df_sttns2_ %>% group_by(CodigoEstacion,NombreEstacion,DEPARTAMENTO,n) %>% 
  summarise(t_gap = sum(gap))

## Histograma del recuento de faltantes en el periodo 1981-1985 
ggplot(df_sttns2_gr, aes(x = t_gap)) +
  geom_histogram(
    breaks = seq(0, 60, by =12),   # define los cortes
    fill   = "lightgray",
    color  = "black"
  ) +
  stat_bin(
    breaks = seq(0, 60, by =12),
    geom = "text",
    aes(label = after_stat(count)), # usa la frecuencia como etiqueta
    vjust = -0.5, # posición vertical del texto
    size =10 # aumenta el tamaño de las etiquetas
  ) +
  labs(
    title = "Histograma del recuento de datos faltantes en el periodo 1981-01 a 1985-12",
    x = "Datos faltantes",
    y = "Frecuencia"
  ) + scale_x_continuous(breaks = seq(0,60, by =12)) +
  scale_y_continuous(limits = c(0,1000)) +
  theme(axis.text=element_text(size=28),
        plot.title = element_text(size=28,face="bold"),
        axis.title=element_text(size=28))

# A2 Gráfico datos faltantes ----
# df_c65:  df de las sttns que no cumplen el criterio c90 pero se incluyen
## pivota las cols de meses (para que queden los NAs en las filas)
## crea la col gap (que se grafica como punto negro), Fecha (eje X) y sttn (eje Y)
## ordena por n (número de obs: sttn-mes)
df_c65 = filter(df_sttns2_train, n<c90) %>%  
  pivot_longer(cols = matches("^(19|20)"), # columnas que comienzan por 19.. o 20..
               names_to  = "Fecha",
               values_to = "pcp") %>% 
  mutate(gap = is.na(pcp),
         Fecha = as.Date(Fecha),
         sttn = paste0(DEPARTAMENTO,' - ',NombreEstacion, ' - ', round((n/480)*100), "%")) %>% 
  arrange(desc(n))


# df_sample: muestra de 94 sttns que cumplen c90

set.seed(01092025)
df_sample = filter(df_sttns2_train, n>=c90) %>%  
  filter(CodigoEstacion %in% sample(CodigoEstacion,94,replace = F)) %>% 
  pivot_longer(cols = matches("^(19|20)"), # columnas que comienzan por 19.. o 20..
               names_to  = "Fecha",
               values_to = "pcp") %>% 
  mutate(gap = is.na(pcp),
         Fecha = as.Date(Fecha),
         sttn = paste0(DEPARTAMENTO,' - ',NombreEstacion, ' - ', round((n/480)*100), "%")) %>% 
  arrange(desc(n))

## crea objeto para ordenar las sttns de mayor a menor completitud de datos
orden_estaciones <- df_c65 %>%
  group_by(sttn) %>%
  summarise(n_gaps = sum(gap), .groups = "drop") %>%
  arrange(desc(n_gaps)) %>%
  pull(sttn)

orden_estaciones <- df_sample %>%
  group_by(sttn) %>%
  summarise(n_gaps = sum(gap), .groups = "drop") %>%
  arrange(desc(n_gaps)) %>%
  pull(sttn)

## redefine la col sttn como factor
df_c65 <- df_c65 %>%
  mutate(sttn = factor(sttn, levels = orden_estaciones))

df_sample <- df_sample %>%
  mutate(sttn = factor(sttn, levels = orden_estaciones))

## Gráfico datos faltantes en sttns que no cumplen c90
ggplot(df_c65, aes(x = Fecha, y = sttn)) +
  # Línea base (continua) donde hay datos
  geom_line(data = filter(df_c65, !gap),
            aes(color = sttn), size = 3) +
  # Puntos negros donde hay NA
  geom_point(data = filter(df_c65, gap),
             color = "black", size = 2) +
  # xticks del eje X cada 5 años
  scale_x_date(
    breaks = seq(as.Date("1980-01-01"), as.Date("2020-12-01"), by = "5 years"),
    date_labels = "%Y"
  ) +
  labs(title = "Gaps in Observed and CHIRPS rainfall data",
       x = "Date",
       y = "Stations",
       color = "Station") +
  theme_minimal() +
  # remueve la leyenda
  theme(legend.position = "none")

## Gráfico datos faltantes en muestra 10%
ggplot(df_sample, aes(x = Fecha, y = sttn)) +
  # Línea base (continua) donde hay datos
  geom_line(data = filter(df_sample, !gap),
            aes(color = sttn), size = 3) +
  # Puntos negros donde hay NA
  geom_point(data = filter(df_sample, gap),
             color = "black", size = 2) +
  # xticks del eje X cada 5 años
  scale_x_date(
    breaks = seq(as.Date("1980-01-01"), as.Date("2020-12-01"), by = "5 years"),
    date_labels = "%Y"
  ) +
  labs(title = "Gaps in Observed and CHIRPS rainfall data",
       x = "Date",
       y = "Stations",
       color = "Station") +
  theme_minimal() +
  # remueve la leyenda
  theme(legend.position = "none")

#tbl_sttns_ = CNE_IDEAM %>% select(all_of(col_IDEAM)) %>% 
#  inner_join(df_sttns %>%  count(CodigoEstacion, NombreEstacion), 
#             by = c('CODIGO'='CodigoEstacion')) %>% 
#  inner_join(select(KM2, -SUPERFICIE_KM2), by = 'DEPARTAMENTO') %>% 
#  count(REGION)

### A3. determina regiones con baja densidad de sttns ----

### Carga superficie de c/DPTO para calcular densidad sttns x kM2 en c/DPTO
KM2 = read_excel("Superficie_DPTO.xlsx")

### tabla resumen conteo sttns aplicando c90 y c65, calcula pct_90  y pct_65 y dnsn_90 y dnsd_65
tbl_dpto_sttns =  count(tbl_sttns2, DEPARTAMENTO) %>% 
  left_join(filter(tbl_sttns2, n>c90) %>% count(DEPARTAMENTO) %>% rename('n>90'='n'), 
            by = "DEPARTAMENTO") %>% 
  left_join(filter(tbl_sttns2, n>c65) %>% count(DEPARTAMENTO) %>% rename('n>65'='n'), 
            by = "DEPARTAMENTO") %>% 
  mutate(`n>90` = ifelse(is.na(`n>90`), 0,`n>90`)) %>% 
  mutate(pct_90 = `n>90`/n) %>%
  mutate(pct_65 = `n>65`/n) %>% 
  left_join(KM2, by = "DEPARTAMENTO") %>% 
  mutate(dnsd_1k_km2 = ifelse(SUPERFICIE_KM2 > 1000,(n/SUPERFICIE_KM2)*1000, (n/SUPERFICIE_KM2)*10),
         dnsd_90 =  ifelse(SUPERFICIE_KM2 > 1000, (`n>90`/SUPERFICIE_KM2)*1000, (`n>90`/SUPERFICIE_KM2)*10),
         dnsd_65 =  ifelse(SUPERFICIE_KM2 > 1000, (`n>65`/SUPERFICIE_KM2)*1000, (`n>65`/SUPERFICIE_KM2)*10),
         pct_chang = (dnsd_90-dnsd_1k_km2)/dnsd_1k_km2)

## criterio densidad de sttns usar Q1 0.64 sttns x 10k KM2
## criterio pct sttns usar cte 0.4
## criterio sttns usar Q2 39 sttns x DPTO
cte_pct = 0.4
Q1_dnsd = quantile(tbl_dpto_sttns$dnsd_1k_km2)[2]
Q2_sttns = quantile(tbl_dpto_sttns$n)[3]

# filtra los DPTOS que cumplen los 3 criterios
### "Amazonas" "Arauca"   "Caqueta"  "Casanare" "Guainía"  "Guaviare" "Vaupes" "Vichada"
tbl_ = filter(tbl_dpto_sttns, pct_90 < cte_pct & dnsd_90 < Q1_dnsd & n< Q2_sttns)
tbl_$DEPARTAMENTO

# chequear densidad por Región  Natural
tbl_2 = tbl_dpto_sttns %>% group_by(REGION) %>%  
  summarise(SUPERFICIE_KM2 = sum(SUPERFICIE_KM2),
            n = sum(n),
            `n>90` = sum(`n>90`),
            `n>65` = sum(`n>65`),
            dnsd = n/SUPERFICIE_KM2,
            dnsd_90 = `n>90`/SUPERFICIE_KM2*1000,
            dnsd_65 = `n>65`/SUPERFICIE_KM2*1000,
            pct_change = (`n>90`-n)/n)

### distribución del NivelAprobacion x sttns ----

df_sttns2 = df_sttns %>% count(CodigoEstacion, NombreEstacion)

df_sttns3 = df_sttns %>% count(CodigoEstacion, NombreEstacion, NivelAprobacion) %>%
  pivot_wider(
    names_from = NivelAprobacion,  
    values_from = n,              
    values_fill = 0               
  )

s_ = sum(df_sttns3$Definitivo) + sum(df_sttns3$Preliminar) + sum(df_sttns3$`En revisión`)
sum(df_sttns3$Definitivo)/s_ 
sum(df_sttns3$Preliminar)/s_
sum(df_sttns3$`En revisión`)/s_

# compara los datos IDEAM con pcp_col ----
### Carga 846 estaciones iniciales
pcp_col <- read_excel("~/TS_climatic/pcp_col.xlsx")
# hay 728 estaciones que están en el diccionario, por tanto hay datos de estas estaciones
pcp_col2 = select(pcp_col, "ID", "X", "Y") %>%  
  inner_join(tbl_sttns3, by = c('ID'='CODIGO'))

dput(names(pcp_col))

pcp_col2 = select(pcp_col, "ID", "X", "Y") %>%  
  inner_join(CNE_IDEAM2, by = c('ID'='CODIGO'))

sort(prop.table(table(pcp_col2$CATEGORIA)))
sort(prop.table(table(pcp_col2$TECNOLOGIA)))
sort(prop.table(table(pcp_col2$ESTADO)))

pcp_col3 = select(pcp_col, "ID", "X", "Y") %>%  
  inner_join(df_sttns3, by = c('ID'='CodigoEstacion')) %>% 
  arrange(desc(`En revisión`), desc(Preliminar), desc(Definitivo)) 

s_ = (sum(pcp_col3$Definitivo)+sum(pcp_col3$Preliminar)+sum(pcp_col3$`En revisión`))
sum(pcp_col3$Definitivo)/s_
sum(pcp_col3$Preliminar)/s_
sum(pcp_col3$`En revisión`)/s_


df_sttns$Fecha = as.Date(df_sttns$Fecha)
# pivota columnas a filas
df_long <- pcp_col %>%
  pivot_longer(
    cols = matches("^19|^20"), 
    names_to = "FECHA",        
    values_to = "VALOR")

df_long$FECHA = as.Date(paste0(df_long$FECHA, '-01'), format = '%Y-%m-%d')

df_long2 = df_long %>% inner_join(df_sttns, by = c('ID'='CodigoEstacion', 'FECHA'='Fecha'))
df_long2 = df_long2 %>% inner_join(select(CNE_IDEAM2, c(CODIGO, DEPARTAMENTO, LATITUD, LONGITUD)), 
                                   by = c('ID'='CODIGO')) 

prop.table(table(df_long2$VALOR==df_long2$Valor)) # 96.5% de coincidencias
df_long2$MSE = (df_long2$VALOR - df_long2$Valor)^2
df_MSE = df_long2 %>% group_by(ID) %>% summarise(RMSE = sqrt(mean(MSE)))
mean(df_MSE$RMSE) # 60.87 mm es el promedio del RMSE
quantile(df_MSE$RMSE) # Q1: 4.06, Q2: 33.2 y Q3: 114.9
# todas las 797 tiene al menos un dato que difiere
View(df_long2 %>%  filter(VALOR!=Valor) %>% count(ID))
View(df_long2 %>%  filter(ID == 21015030))
View(df_long2 %>%  filter(ID == 16020140 & VALOR!=Valor))
