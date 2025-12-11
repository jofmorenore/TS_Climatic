
library(dplyr)
library(lubridate)
library(factoextra)
library(cluster)
library(sf)
library(gridExtra)
library(tidyr)
library(gridExtra)
library(grid)

library (geofd)
library(sf)
library(fda)
library(maps)
library(dplyr)
library(ggplot2)

setwd('~/TS_climatic/')
pol2_list <- readRDS('sttns_validation_chirps.RDS')

# obtiene media multi-anual
list_pcp_annual_mean = list()
geometry_ = list()

for(i in 1:length(pol2_list)){
  list_pcp_annual_mean[[i]] = dplyr::select(pol2_list[[i]], ID, Date, sttns_fill) %>% 
    group_by(ID, year(Date)) %>% 
    summarise(annual_pcp = sum(sttns_fill)) %>% 
    summarise(mean(annual_pcp))
  geometry_[[i]] = unique(pol2_list[[i]]$`geometry-x`)
}

df_pcp_annual_mean = bind_rows(list_pcp_annual_mean)
for(i in 1:length(geometry_)){
  df_pcp_annual_mean$X[i] = geometry_[[i]][[1]][1]
  df_pcp_annual_mean$Y[i] = geometry_[[i]][[1]][2]
}

# crea variable escalada
df_pcp_annual_mean$pcp_scaled = scale(df_pcp_annual_mean$`mean(annual_pcp)`) 
write.table(df_pcp_annual_mean, '~/df_pcp_annual_mean.csv', sep = ',', col.names = T, row.names = F)
# 1. K-means -------------------------------------------------------------------

# Wss
fnb1 <- fviz_nbclust(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, method = "wss") +
  geom_hline(yintercept = km_list[[4]]$tot.withinss, linetype = 2) 

  # de 3  a 4 grupos
# Silhouette
fnb2 <- fviz_nbclust(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, 
                     method = "silhouette", k.max = 20) # 2 o 3 grupos
# gap_stat
# compute gap statistic
set.seed(12)
gap_stat <- clusGap(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, nstart = 25,
                    K.max = 25, B = 50)
print(gap_stat, method = "firstmax")

set.seed(12)
fnb3 <- fviz_nbclust(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, method = "gap_stat", k.max = 20) + 
  geom_hline(yintercept = 0.81, linetype = 2) +
  coord_cartesian(xlim = c(1.75,20)) # 2, 13 o 16 grupos

grid.arrange(fnb1, fnb2, fnb3, ncol = 3)

df_ = dplyr::select(df_pcp_annual_mean, pcp_scaled)
row.names(df_) = df_pcp_annual_mean$ID

n_optim = c(2,3,4,13,16)
km_list = list()

for(k in n_optim){
  km_list[[k]] = kmeans(df_, k, nstart = 25)
}

print(km_list[[4]])

df_pcp_annual_mean = df_pcp_annual_mean %>% 
  mutate(clust2 = factor(km_list[[2]]$cluster),
         clust3 = factor(km_list[[3]]$cluster),
         clust4 = factor(km_list[[4]]$cluster),
         clust13 = factor(km_list[[13]]$cluster))

ggplot(df_pcp_annual_mean, aes(X, Y, color = clust2)) +
  geom_point(alpha = 0.25)

df_km_pcp_annual_mean = st_as_sf(df_pcp_annual_mean, 
                                 coords = c('X', 'Y'),
                                 crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

shp_depto <- st_read("Departamentos202208_shp/Depto.shp")

map2 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_km_pcp_annual_mean, 'Cluster'=clust2), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))

map4 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_km_pcp_annual_mean, 'Cluster'=clust4), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))

map13 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_km_pcp_annual_mean, 'Cluster'=clust13) , 
          aes(fill = Cluster, pch= Cluster, color = Cluster), size = 1, stroke = 1) +
  scale_shape_manual(values = c(0:5, 0:6))+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.8))

fviz_cluster(km_list[[4]], 
             data = dplyr::select(df_pcp_annual_mean, X, Y))

bx2 <- ggplot(rename(df_km_pcp_annual_mean, 'Cluster'=clust2), 
              aes(x = Cluster, y = `mean(annual_pcp)`, fill = Cluster)) +
  geom_boxplot() + 
  scale_y_continuous(breaks = seq(0, 13000, 2000)) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))+
  labs(y = "Precipitación anual media")

bx4 <- ggplot(rename(df_km_pcp_annual_mean, 'Cluster'=clust4), 
              aes(x = Cluster, y = `mean(annual_pcp)`, fill = Cluster)) +
  geom_boxplot() + 
  scale_y_continuous(breaks = seq(0, 13000, 2000)) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))+
  labs(y = "Precipitación anual media")

bx13 <- ggplot(rename(df_km_pcp_annual_mean, 'Cluster'=clust13), 
              aes(x = Cluster, y = `mean(annual_pcp)`, fill = Cluster)) +
  geom_boxplot() + 
  scale_y_continuous(breaks = seq(0, 13000, 2000)) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.8))+
  labs(y = "Precipitación anual media")

grid.arrange(map2, bx2, ncol = 2)
table(df_km_pcp_annual_mean$clust4, df_km_pcp_annual_mean$clust2)
grid.arrange(map4, bx4, ncol = 2)
grid.arrange(map13, bx13, ncol = 2)

pct2_ <- df_km_pcp_annual_mean %>% 
  count(clust2) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust2)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))  

pct4_ <- df_km_pcp_annual_mean %>% 
  count(clust4) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust4)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            size = 5, color = "white") +
  labs(title = "Porcentaje de estaciones metereológicas en cada cluster",
       fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))  

pct13_ <- df_km_pcp_annual_mean %>% 
  count(clust13) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust13)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))  

grid.arrange(pct2_, pct4_, pct13_, ncol = 3)

# 2. Random Forest -----

library(MASS)
library(cluster)
library(survival)
library(randomForest)
library(Hmisc)

## Rand: función ...
if (exists("Rand") ) rm(Rand)
Rand <- function(tab,adjust=T) {
  
  ##########################################################################
  # The function computes the (adjusted) Rand index between two partitions #
  # Copyright Steve Horvath and Luohua Jiang, UCLA, 2003                   #
  ##########################################################################
  
  # helper function
  choosenew <- function(n,k) {
    n <- c(n); out1 <- rep(0,length(n));
    for (i in c(1:length(n)) ){
      if ( n[i]<k ) {out1[i] <- 0}
      else {out1[i] <- choose(n[i],k) }
    }
    out1
  }
  
  a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
  n <- nrow(tab)
  for (i in 1:n) {
    for(j in 1:n) {
      a <- a+choosenew(tab[i,j],2)
      nj <- sum(tab[,j])
      c <- c+choosenew(nj,2)
    }
    ni <- sum(tab[i,])
    b <- b+choosenew(ni,2)
    nn <- nn+ni
  }
  if(adjust==T) {
    d <- choosenew(nn,2)
    adrand <- (a-(b*c/n)/d)/(0.5*(b+c/n)-(b*c/n)/d)
    adrand
  } else {
    b <- b-a
    c <- c/n-a
    d <- choosenew(nn,2)-a-b-c
    rand <- (a+d)/(a+b+c+d)
    rand
  }
}

## pamNew: función ...
if (exists("pamNew") ) rm(pamNew)
pamNew <- function (x, k, diss1 = inherits(x, "dist"), metric1 = "euclidean"){
  
  #############################################################################
  # PAM: Partitioning Around Medoids
  # @x: matriz o dist
  # @k: num de clusters
  # @diss1: logical si x es un objeto de dist
  # @metric1: dist usada cuando x es una matriz
  #############################################################################
  
  if (diss1)
  { # si x es un objeto dist
    if (!is.null(attr(x, "Labels"))) { original.row.names <- attr(x, "Labels")}
    names(x) <- as.character(c(1:attr(x, "Size")))
  } 
  else
  { # si x es una matriz
    if(!is.null(dimnames(x)[[1]])) { original.row.names <- dimnames(x)[[1]]}
    row.names(x) <- as.character(c(1:dim(x)[[1]]))
  }
  # ejecuta PAM estándar
  pam1 <- pam(x,k,diss=diss1, metric=metric1)
  
  label2 <- pam1$clustering # asignación incial de clusters
  # extrae y ordena silueta
  silinfo1 <- pam1$silinfo$widths #| 1 clust asignado | 2 'neighbor clust' | 3 ancho silueta
  index1 <- as.numeric(as.character(row.names(silinfo1)))
  silinfo2 <- silinfo1[order(index1),]
  # reasigna segun silueta
  labelnew <- ifelse(silinfo2[,3]<0, silinfo2[,2], silinfo2[,1])
  names(labelnew) <- original.row.names
  out <- list()
  out[['label']] = labelnew # vector de salida con las nuevas etiquetas del cluster    
  out[['avg_silhouette_width']] = mean(silinfo2[,3])
  return(out)
  #############################################################################################################
  # A new pam clustering function which corrects the clustering membership based on the sillhouette strength. #
  # The clustering membership of an observation with a negative sillhouette strength is reassigned to its     #
  # neighboring cluster.                                                                                      #
  # The inputs of the function are similar to the original 'pam' function.                                    #
  # The function returns a vector of clustering labels.                                                       #
  # Copyright 2003 Tao Shi and Steve Horvath (last modified 10/31/03)                                         #
  #############################################################################################################
}

## collect.garbage: función ...
if (exists("collect.garbage") ) rm(collect.garbage)
collect.garbage <- function(){
  ## The following function collects garbage until the memory is clean.
  ## Usage: 1. immediately call this function after you call a function or
  ##        2. rm()
  while (gc()[2,4] != gc()[2,4]){}
}

## RFdist: función ...
if (exists("RFdist") ) rm(RFdist)
# datRF = mtrx 
# mtry1 = 4 # num de variables en c/división
# no.tree = 500 # num de árboles
# no.rep = 50 # num de repeticiones/ num de bosques a promediar
# addcl1 = T 
# addcl2 = F 
# imp = T 
# oob.prox1 = T
# 
# rm(datRF, mtry1, no.tree, no.rep, addcl1, addcl2, imp, oob.prox1)
# rm(sample1, g1, dat, nrow1, ncol1, RFproxAddcl1, RFproxAddcl2,
#    RFprox1Conver, RFprox2Conver, RFimportance1, RFimportance2,
#    RFerrrate1, RFerrrate2, rep1,i, index1)
# 
 RFdist <- function(datRF, mtry1, no.tree, no.rep, addcl1=T, addcl2=T, imp=T, 
                    oob.prox1=T, proxConver=F, seed = 123) {
  
  ####################################################################
  # datRF: data.frame o matrix n x p (sin col respuesta?)
  # mtry1: int núm. de variables probadas en cada división
  # no.tree: int núm. de árboles (ntree) que se pasan a randomForest
  # no.rep: int núm. de repeticiones / número de bosques a promediar. 
  # addcl1: logical (default TRUE) ejecuta el procedimiento con el método synthetic1 
  ##        (re-muestreo de cada variable por permutación — muestreo con reemplazo de las columnas).
  # addcl2: logical (default TRUE) ejecuta el procedimiento con el método synthetic2 
  ##        (generación de variables aleatorias uniformes en el rango observado por columna).
  # imp: logical (default TRUE) acumula/importa medidas de importancia de variables 
  ##      (RF$importance) promediadas sobre repeticiones.
  # oob.prox1: logical (default TRUE) calcula proximities OOB.
  # proxConver: logical (default FALSE) la función acumula y devuelve 
  ##            estadísticas de convergencia de la proximidad a lo largo de las repeticiones.
  ####################################################################
  
  synthetic1 <- function(dat) {
    # Construye un dataset sintético donde la mitad de las filas son las originales y 
    # la otra mitad son versiones permutadas por columna.
    sample1 <- function(X)   { sample(X, replace=T) } 
    g1      <- function(dat) { apply(dat,2,sample1) }
    nrow1 <- dim(dat)[[1]];
    yy <- rep(c(1,2),c(nrow1,nrow1) );
    data.frame(cbind(yy,
                     rbind(data.frame(dat),
                           data.frame(g1(dat)))))
  }
  
  synthetic2 <- function(dat) {
    # Construye dataset sintético donde la mitad son las originales y la otra mitad son 
    # valores aleatorios runif distribuidos uniformemente entre el mínimo y máximo 
    # observados de cada columna.
    sample2 <- function(X)   { runif(length(X), min=min(X), max =max(X)) }
    g2      <- function(dat) { apply(dat,2,sample2) }
    nrow1 <- dim(dat)[[1]];
    yy <- rep(c(1,2),c(nrow1,nrow1) );
    data.frame(cbind(yy,rbind(dat,data.frame(g2(dat)))))
  }
  
  cleandist <- function(x) { 
    # convierte a objeto dist, sustituye cualquier valor ≤ 0 por un valor muy pequeño (1e-10) 
    ## (evita distancias cero/negativas)
    x1 <- as.dist(x)
    x1[x1<=0] <- 0.0000000001
    as.matrix(x1)
  }
  
  nrow1 <- dim(datRF)[[1]]# num filas
  ncol1 <- dim(datRF)[[2]] # num cols
  RFproxAddcl1 <- matrix(0,nrow=nrow1,ncol=nrow1) # matrix n x n
  RFproxAddcl2 <- matrix(0,nrow=nrow1,ncol=nrow1) # matrix n x n
  RFprox1Conver <- cbind(1:no.rep,matrix(0,(no.rep),3)) # matrix no.rep x p
  RFprox2Conver <- cbind(1:no.rep,matrix(0,(no.rep),3)) # matrix no.rep x p
  # se asume 4 cols de importancia; esto asume la estructura que devuelve randomForest::importance
  RFimportance1 <- matrix(0, nrow=ncol1, ncol=4) # matrix p x 4
  RFimportance2 <- matrix(0, nrow=ncol1, ncol=4) # matrix p x 4
  RFerrrate1 <- 0
  RFerrrate2 <- 0
  rep1 <- rep(666,2*nrow1) # vector para reordenar ídx del objeto proximity devuelto por randomForest
  
  #addcl1
  if (addcl1) { # (si addcl1 = TRUE)
    for (i in c(0:no.rep)) { 
      set.seed(seed)
      index1 <- sample(c(1:(2*nrow1))) # crea idx (orden) aleatorio 
      rep1[index1] <-  c(1:(2*nrow1)) 
      datRFsyn <- synthetic1(datRF)
      datRFsyn <- synthetic1(datRF)[index1,] # asigna df sintético reordenado con index1
      yy <- datRFsyn[,1] # respuesta
      # ajusta un RandomForest
      RF1 <- randomForest(factor(yy)~.,data=datRFsyn[,-1], 
                          ntree=no.tree, 
                          oob.prox=oob.prox1, 
                          proximity=TRUE, # matriz de proximidad
                          do.trace=F,
                          mtry=mtry1,
                          importance=imp) 
      collect.garbage()
      RF1prox <- RF1$proximity[rep1,rep1] # reordena matriz de proximidad 
      if (i > 0) { 
        if (i > 1){
          # calcula métricas de convergencia
          # dif entre media acum de prox con la nueva replica y la media anterior
          xx <- ((RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]))/i) - (RFproxAddcl1/(i-1))
          yy <- mean( c(as.dist((RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]))/i))) 
          RFprox1Conver[i,2] <- max(abs(c(as.dist(xx)))) # max dif absoluta entre pares
          RFprox1Conver[i,3] <- mean((c(as.dist(xx)))^2) # meadia de cuadrados de dif 
          RFprox1Conver[i,4] <- yy # prox promedio
        }
        RFproxAddcl1 <- RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]) 
        if(imp) { RFimportance1 <- RFimportance1+ 1/no.rep*(RF1$importance) }
        RFerrrate1 <- RFerrrate1+ 1/no.rep*(RF1$err.rate[no.tree])
      }
    }
  }
  
  # addcl2
  if (addcl2) {  #(si addcl2 = TRUE)
    for (i in c(0:no.rep)) {
      set.seed(seed)
      # mismo esquema que addcl1
      index1 <- sample(c(1:(2*nrow1))) 
      rep1[index1] <-  c(1:(2*nrow1)) 
      # con synthetic2 usa val uniformes entre min y max por col en lugar de pertumtaciones por col
      datRFsyn <- synthetic2(datRF)[index1,] 
      yy <- datRFsyn[,1] 
      RF2 <- randomForest(factor(yy)~.,data=datRFsyn[,-1], 
                          ntree=no.tree, 
                          oob.prox=oob.prox1, 
                          proximity=TRUE,
                          do.trace=F,
                          mtry=mtry1,
                          importance=imp) 
      collect.garbage()
      RF2prox <- RF2$proximity[rep1,rep1]
      if (i > 0) { 
        if (i > 1){
          xx <- ((RFproxAddcl2 + (RF2prox[c(1:nrow1),c(1:nrow1)]))/i) - (RFproxAddcl2/(i-1))
          yy <- mean( c(as.dist((RFproxAddcl2 + (RF2prox[c(1:nrow1),c(1:nrow1)]))/i))) 
          RFprox2Conver[i,2] <- max(abs(c(as.dist(xx))))
          RFprox2Conver[i,3] <- mean((c(as.dist(xx)))^2)
          RFprox2Conver[i,4] <- yy
        }
        RFproxAddcl2 <- RFproxAddcl2 + (RF2prox[c(1:nrow1),c(1:nrow1)]) 
        if(imp) { RFimportance2 <- RFimportance2+ 1/no.rep*(RF2$importance)}
        RFerrrate2 <- RFerrrate2+ 1/no.rep*(RF2$err.rate[no.tree])
      }
    }
  }
  # promedia distancia: RFproxAddcl1/no.rep
  # convierte a distancia sqrt(1 - prox)
  # cleandist: sustituye cualquier valor ≤ 0 por un valor muy pequeño (1e-10) 
  distRFAddcl1 <- cleandist(sqrt(1-RFproxAddcl1/no.rep))
  distRFAddcl2 <- cleandist(sqrt(1-RFproxAddcl2/no.rep))
  
  distRF <- list(cl1=NULL, err1=NULL, imp1=NULL, prox1Conver=NULL, 
                 cl2=NULL, err2=NULL, imp2=NULL, prox2Conver=NULL)
  
  if(addcl1) {
    distRF$cl1 <- distRFAddcl1 # matriz de dist n x n 
    distRF$err1 <- RFerrrate1 # num tasa de error oob promedio 
    if(imp) distRF$imp1 <- RFimportance1 # matriz de importancia de var promediada p x 4
    if(proxConver) distRF$prox1Conver <- RFprox1Conver # matriz de convergencia i: 1,..., no.rep x j: iter | max abs change | mean squared change | mean prox
  }
  if(addcl2) {
    distRF$cl2 <- distRFAddcl2
    distRF$err2 <- RFerrrate2
    if(imp) distRF$imp2 <- RFimportance2
    if(proxConver) distRF$prox2Conver <- RFprox2Conver
  } 
  
  ####################################################################
  # Unsupervised randomForest function                               #
  # Return a list "distRF" containing some of the following 6 fields #
  #  depending on the options specified:                             #
  #  (1) cl1:  addcl1 distance (sqrt(1-RF.proxAddcl1))               #
  #  (2) err1: error rate                                            #
  #  (3) imp1: variable importance for addcl1                        #
  #  (4) prox1Conver: a matrix containing two convergence meausres   #
  #                   for addcl1 proximity                           #
  #                   a). max( abs( c(aveprox(N))-c(aveprox(N-1))))  #
  #                   b). mean((c(aveprox(N))-c(aveprox(N-1)))^2)    #
  #                   where N is number of forests (no.rep).         #
  #  (5) cl2, (6) err2, (7)imp2 and (8) prox2Conver for addcl2       #
  # Copyright Steve Horvath and Tao Shi (2004)                       #
  ####################################################################
  distRF
}

# Cluster using Random Forest ---------------------------------------------

# matriz de normal climática mensual
list_pcp_monthly_normal = list()
geometry_ = list()

for(i in 1:length(pol2_list)){
  list_pcp_monthly_normal[[i]] = dplyr::select(pol2_list[[i]], ID, Date, sttns_fill) %>% 
    group_by(ID, month(Date)) %>% 
    summarise(normal_monthly_pcp = mean(sttns_fill)) %>% 
    pivot_wider(id_cols = ID, names_from = 'month(Date)', values_from ='normal_monthly_pcp')
  geometry_[[i]] = unique(pol2_list[[i]]$`geometry-x`)
}

df_pcp_monthly_normal = bind_rows(list_pcp_monthly_normal)

for(i in 1:length(geometry_)){
  df_pcp_monthly_normal$X[i] = geometry_[[i]][[1]][1]
  df_pcp_monthly_normal$Y[i] = geometry_[[i]][[1]][2]
}

mtrx <- as.matrix(df_pcp_monthly_normal[,2:(ncol(df_pcp_monthly_normal)-2)])
nfrs <- 50
ntrs <- 500 # 1000
# mtry1: sqrt(p) -> sqrt(12) = 3.4
rfds <- RFdist(mtrx, 
               mtry1 = 4, # num de variables en c/división
               ntrs, # num de árboles
               nfrs, # num de repeticiones/ num de bosques a promediar
               addcl1 = T, 
               addcl2 = F, 
               imp = T, 
               oob.prox1 = T)

#lbrf <- pamNew(rfds$cl1, 5)

lbrf = list()
for(k in 2:16){
  lbrf[[k]] <- pamNew(rfds$cl1, k)
}

avg_silhouette_width = c()
for(i in 2:length(lbrf)){
  avg_silhouette_width = c(avg_silhouette_width,
                           lbrf[[i]]$avg_silhouette_width)
}

df_RF_ = data.frame(k = 2:16, avg_silhouette_width)

ggplot(df_RF_, aes(x = k, y = avg_silhouette_width))+
  geom_line() +
  geom_point()+
  scale_x_continuous(breaks = seq(1, 16, 1)) +
  labs(x = 'Number of clusters k',
       y = 'Average Silhouette width')+
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=24))
  
table(lbrf[[2]]$label)
table(lbrf[[4]]$label)

df_pcp_monthly_normal = df_pcp_monthly_normal %>% ungroup() %>%  
  mutate(clust2_RF = factor(lbrf[[2]]$label),
         clust4_RF = factor(lbrf[[4]]$label))

df_RF_pivot_longer = df_pcp_monthly_normal %>% 
  pivot_longer(cols = 2:13, names_to = 'mes', values_to = 'pcp_monthly_mean') %>% 
  mutate(mes = factor(mes, levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))

df_RF_pcp_monthly_normal = st_as_sf(df_pcp_monthly_normal, 
                                    coords=c('X', 'Y'), 
                                    crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

map2_RF <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_RF_pcp_monthly_normal, 'Cluster'=clust2_RF), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))

map4_RF <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_RF_pcp_monthly_normal, 'Cluster'=clust4_RF), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))

bx2_RF <- ggplot(rename(df_RF_pivot_longer, 'Cluster'=clust2_RF), 
              aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))+
  labs(y = "Normal precipitación mensual")

bx4_RF <- ggplot(rename(df_RF_pivot_longer, 'Cluster'=clust4_RF), 
              aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))+
  labs(y = "Normal precipitación mensual")

grid.arrange(map2_RF, bx2_RF, ncol = 2)
grid.arrange(map4_RF, bx4_RF, ncol = 2)

pct2_RF <- df_pcp_monthly_normal %>% 
  count(clust2_RF) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust2_RF)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))  

pct4_RF <- df_pcp_monthly_normal %>% 
  count(clust4_RF) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust4_RF)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))  

grid.arrange(pct2_RF, pct4_RF, ncol = 2,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 16, fontface = "bold")
             ))

# 3. Hierachical clustering of spatially correlated fd -------------------------

data(CanadianWeather, package="fda")
str(CanadianWeather, max.level = 1)

# como Longitud de Canada es West, es negativo  
CanadianWeather$coordinates[,2]=-CanadianWeather$coordinates[,2]
head(CanadianWeather$coordinates)

CanadianWeather_planar = matrix(NA, nr=nrow(CanadianWeather$coordinates), nc=2)

for(i in 1:nrow(CanadianWeather$coordinates)){
  CanadianWeather_planar[i,] = as.data.frame(CanadianWeather$coordinates)[i,] %>%
    st_as_sf(.,coords=c("W.longitude", "N.latitude"), crs = 4326) %>%
    st_transform(crs = convertlatlong2UTM(.)) %>%
    st_coordinates()}

CanadianWeather_planar = as.data.frame(CanadianWeather_planar) %>% `colnames<-`(c("W.longitude", "N.latitude"))

# Smooth temp values using a Fourier basis with 65 functions 
nbasis <- 65
Temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
day.range <- c(1, 365)
Temp.basis <- create.fourier.basis(rangeval = day.range, nbasis = nbasis)
Temp.fd <- smooth.basis(argvals = 1:365, y = Temp, fdParobj = Temp.basis)$fd
plot(Temp.fd)

#daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65, period = 365,
#                                   axes=list('axesIntervals'))
#Temp.fd <- with(CanadianWeather, smooth.basisPar(day.5,
#                                                 dailyAv[,,'Temperature.C'], daybasis65)$fd)
#
#plot(Temp.fd)

# Remover la tendencia espacial con el modelo de regresión funcional 
## X_i(t) = \alpha(t) + \beta_1(t) Longitude_i + \beta_2(t) Latitude_i + e_i (t)

coord = as.data.frame(CanadianWeather$coordinates)
names(coord)
fRegress(CanadianWeather$coordinates)

TempRgn.f <- fRegress(Temp.fd ~ N.latitude + W.longitude , coord)
#TempRgn.f2 <- fRegress(Temp.fd ~ N.latitude + W.longitude , CanadianWeather_planar)
#fdobj.res = TempRgn.f$yfdobj-TempRgn.f$yhatfdobj
fdobj.res = Temp.fd-TempRgn.f$yhatfdobj
plot(fdobj.res)
#summary(TempRgn.f)
fdobj.res

# evalua los residuales
day_grid <- 1:365
res_M <- eval.fd(day_grid, fdobj.res)


# Fit a spherical model to the estimated trace-variogram by using the OLS technique
coords.cero <- data.frame(Lon = -64.06, Lat = 45.79)

okfd.res <- okfd(new.coords = coords.cero, coords=coord,
                 data=res_M, smooth.type='fourier', nbasis=65, argvals=day.5,
                 fix.nugget=FALSE, fix.kappa=FALSE)

okfd.res$trace.vari.array[[1]] # 8034.03 y 22.08
okfd.res$trace.vari.array[[4]] # 8034.03 y 22.08

plot(okfd.res)

M <- fourierpen(fdobj.res$basis,  Lfdobj=0)

res.fd <- smooth.basis(argvals = 1:365, y = res_M, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res$coefs), res.fd, M)

new.emp.trace.vari1 <- trace.variog(coords=coord,
                                    L2norm=L2norm, bin=FALSE)

new.emp.trace.vari2 <- trace.variog(coords=CanadianWeather_planar,
                                    L2norm=L2norm, bin=FALSE)

fit1 = geofd::fit.tracevariog(new.emp.trace.vari1, models = "spherical",
                              sigma2.0 = 7769, phi.0 = 21,
                              fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                              max.dist.variogram=NULL)

fit2 = geofd::fit.tracevariog(new.emp.trace.vari2, models = "spherical",
                              sigma2.0 = 9000, phi.0 = 2200,
                              fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                              max.dist.variogram=NULL)


fit1$best$cov.pars # 8033, 22.08
fit2$best$cov.pars # 6975, 22.08


plot(new.emp.trace.vari1$u, new.emp.trace.vari1$v)


# Ejemplo geofd: An R Package for Function-Valued Geostatistical Prediction ----
data(maritimes.coords)
data(maritimes.data)

head(maritimes.coords)
head(maritimes.data[,1:4], n=5)
# data set is smoothed by using a B-splines basis with 65 functions without penalization

n <- dim(maritimes.data)[1]
argvals<-seq(1,n, by=1)

# parameters for smoothing the data
s<-35
rangeval <- range(argvals)
norder <- 4
nbasis <- 65
bspl.basis <- create.bspline.basis(rangeval, nbasis, norder)
lambda <-0
datafdPar <- fdPar(bspl.basis, Lfdobj=2, lambda)
smfd <- smooth.basis(argvals,maritimes.data,datafdPar)
datafd <- smfd$fd

# calculate the L2 norm between the smoothed curves

# s: number of sites where curves are observed, 
# datafd: a functional data object representing a smoothed data set
# M: a symmetric matrix of order equal to the number of basis functions defined by the B-splines basis object

M <- bsplinepen(bspl.basis,Lfdobj=0)
L2norm <- l2.norm(s, datafd, M)

# coords: the geographical coordinates in decimal degrees
# L2norm: a matrix whose values are the L2 norm between all pair of smoothed functions (an output from the function l2.norm)
# bin: which is a logical argument indicating whether the output is a binned variogram, 
# maxdist: a numerical value defining the maximum distance for calculating the trace-variogram.
# uvec, breaks and nugget.tolerance: defined as in the function variog of the package geoR.

dista=max(dist(maritimes.coords))*0.9
tracev=trace.variog(maritimes.coords, L2norm, bin=FALSE,
                    max.dist=dista,uvec="default",breaks="default",nugget.tolerance)

# fit a theoretical model to the estimated trace-variogram

# tracev: estimations of the trace-variogram function (an output of the function trace.variog)
# model: list vwith the models that we want to fit
# some initial values for the parameters in these models

models=fit.tracevariog(tracev, models=c("spherical","exponential",
                                        "gaussian","matern"),sigma2.0=2000, phi.0=4, fix.nugget=FALSE,
                       nugget=0, fix.kappa=TRUE, kappa=1, max.dist.variogram=dista)

models$fitted[[1]] # spherical, cov.pars: 2112.12   5.66     
models$fitted[[2]] # exponential, cov.pars: 3333.32   5.12    
models$fitted[[3]] # gaussian, cov.pars: 1810.90    2.51
models$fitted[[4]] # matern, cov.pars: 2346.85 1.83
