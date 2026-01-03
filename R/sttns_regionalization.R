
library(dplyr)
library(lubridate)
library(factoextra)
library(cluster)
library(sf)
library(gridExtra)
library(tidyr)
library(gridExtra)
library(grid)

library(MASS)
library(cluster)
library(survival)
library(randomForest)
library(Hmisc)

library(dtwclust)
library(doParallel)
library(stringr)
library(ggplot2)
library(clue)

library (geofd)
library(sf)
library(sp)
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

shp_depto <- st_read("Departamentos202208_shp/Depto.shp")

# crea variable escalada
df_pcp_annual_mean$pcp_scaled = scale(df_pcp_annual_mean$`mean(annual_pcp)`) 

c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#FB61D7")

#scale_shape_manual(values=c(3, 16, 17))+
#scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+

#write.table(df_pcp_annual_mean, '~/df_pcp_annual_mean.csv', sep = ',', col.names = T, row.names = F)
# 1. K-means -------------------------------------------------------------------

# Wss
fnb1 <- fviz_nbclust(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, method = "wss") +
  geom_hline(yintercept = km_list[[4]]$tot.withinss, linetype = 2) +
  my_theme_bx()+
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16))

  # de 3  a 4 grupos
# Silhouette
fnb2 <- fviz_nbclust(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, 
                     method = "silhouette", k.max = 20) + 
  my_theme_bx()+
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16)) # 2 o 3 grupos
# gap_stat
# compute gap statistic
set.seed(12)
gap_stat <- clusGap(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, nstart = 25,
                    K.max = 25, B = 50)
print(gap_stat, method = "firstmax")

set.seed(12)
fnb3 <- fviz_nbclust(df_pcp_annual_mean$pcp_scaled, FUN = kmeans, method = "gap_stat", k.max = 20) + 
  geom_hline(yintercept = 0.81, linetype = 2) +
  coord_cartesian(xlim = c(1.75,20))+
  my_theme_bx()+
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16)) # 2, 13 o 16 grupos

grid.arrange(fnb1, fnb2, fnb3, ncol = 3)

df_ = dplyr::select(df_pcp_annual_mean, pcp_scaled)
row.names(df_) = df_pcp_annual_mean$ID

n_optim = c(2,3,4,13,16)
km_list = list()

for(k in n_optim){
  km_list[[k]] = kmeans(df_, k, nstart = 25)
}


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

my_theme <- function() {
  theme(
    legend.position = c(0.9, 0.9),
    axis.text = element_text(size = 24),
    axis.text.x = element_text(size = 24, hjust = 1, angle = 90),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )
}

map2 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_km_pcp_annual_mean, 'Cluster'=clust2), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(16, 17))+
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+
  theme_minimal()+
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))

map4 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_km_pcp_annual_mean, 'Cluster'=clust4), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(16, 17, 18,4))+
  scale_color_manual(values=c("#00BFC4", "#F8766D", "#A3A500", "#C77CFF" ))+
  theme_minimal()+
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))

map13 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_km_pcp_annual_mean, 'Cluster'=clust13) , 
          aes(fill = Cluster, pch= Cluster, color = Cluster), size = 1, stroke = 1) +
  scale_shape_manual(values = c(8, 6, 18, 4, 
                                15, 16, 
                                3, 10, 
                                17, 0, 1, 2, 5))+
  scale_color_manual(values=c("#C77CFF","#Ffa500","#A3A500", 
                              "#D12405", "#39B600", "#00BFC4",
                             "#FF62BC", "#000000","#F8766D", 
                             "#C0C0C0", "#619CFF", "#C2E5EB", "#CD7F32"))+
  theme_minimal()+
  my_theme()+
  theme(legend.position = c(0.9, 0.8))+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2),
         shape = guide_legend(ncol = 2))


fviz_cluster(km_list[[4]], 
             data = dplyr::select(df_pcp_annual_mean, X, Y))

my_theme_bx <- function(){
  theme(
    legend.position = c(0.9, 0.9),
    axis.text = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )
}

bx2 <- ggplot(rename(df_km_pcp_annual_mean, 'Cluster'=clust2), 
              aes(x = Cluster, y = `mean(annual_pcp)`, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#00BFC4", "#F8766D"))+
  scale_y_continuous(breaks = seq(0, 13000, 2000)) +
  theme_minimal()+
  labs(y = "Precipitación anual media") +
  my_theme_bx()

bx4 <- ggplot(rename(df_km_pcp_annual_mean, 'Cluster'=clust4), 
              aes(x = Cluster, y = `mean(annual_pcp)`, fill = Cluster)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#A3A500", "#C77CFF"))+
  scale_y_continuous(breaks = seq(0, 13000, 2000)) +
  theme_minimal()+
  labs(y = "Precipitación anual media")+
  my_theme_bx()

bx13 <- ggplot(rename(df_km_pcp_annual_mean, 'Cluster'=clust13), 
              aes(x = Cluster, y = `mean(annual_pcp)`, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#C77CFF","#Ffa500","#A3A500", 
                             "#D12405", "#39B600", "#00BFC4",
                             "#FF62BC", "#000000","#F8766D", 
                             "#C0C0C0", "#619CFF", "#C2E5EB", "#CD7F32"))+
  scale_y_continuous(breaks = seq(0, 13000, 2000)) +
  theme_minimal()+
  labs(y = "Precipitación anual media")+
  my_theme_bx()+
  theme(legend.position = c(0.9, 0.7))+
  guides(color = guide_legend(ncol = 2),
         shape = guide_legend(ncol = 2))
  
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
            fontface = "bold", size = 5, color = "white") +
  scale_fill_manual(values=c("#00BFC4", "#F8766D"))+
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24))


pct4_ <- df_km_pcp_annual_mean %>% 
  count(clust4) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust4)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface = "bold", size = 5, color = "white") +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#A3A500", "#C77CFF"))+
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24))

pct13_ <- df_km_pcp_annual_mean %>% 
  count(clust13) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust13)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            size = 4, color = "white") +
  scale_fill_manual(values=c("#C77CFF","#Ffa500","#A3A500", 
                             "#D12405", "#39B600", "#00BFC4",
                             "#FF62BC", "#000000","#F8766D", 
                             "#C0C0C0", "#619CFF", "#C2E5EB", "#CD7F32"))+
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24))

grid.arrange(pct2_, pct4_, pct13_, ncol = 3,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

# 2. Random Forest -----

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
    data.frame(cbind(yy,
                     rbind(data.frame(dat),
                           data.frame(g2(dat)))))
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
mtrx = scale(mtrx) 

# RandomForest se puede ahcer son una sola columna? R. No arroja error en  factor(yy) ~ ., data = datRFsyn[, 
#mtrx_1 <- as.matrix(df_pcp_annual_mean[,5])

#nfrs <- 50
#ntrs <- 500 # 1000
#mtry1 <- 4 

## 2.1 addcl1: re-muestreo de cada variable por permutación — muestreo con reemplazo de las columnas -----

nfrs_ = 100 #c(10, 50, 100)#seq(10, 100, 20)
ntrs_ = c(100,500,1000,2000) #seq(100, 1000, 200)
mtry1_ = c(3,4,6,8) #seq(2, 10, 2)

# mtry1: sqrt(p) -> sqrt(12) = 3.4
out_ = list()

for(i in mtry1_){
  for(j in ntrs_){
    for(h in nfrs_){
      rfds_ <- RFdist(mtrx, 
                      mtry1 = i, # num de variables en c/división
                      no.tree = j, # num de árboles
                      no.rep = h, # num de repeticiones/ num de bosques a promediar
                      addcl1 = T, 
                      addcl2 = F, 
                      imp = T, 
                      oob.prox1 = T)
      
      avg_ = c()
      for(k in 2:16){
        lbrf <- pamNew(rfds_$cl1, k)
        avg_ = c(avg_, lbrf$avg_silhouette_width)
      }
      
      nam_ = paste(i,j,h, sep ='|')  
      out_[[nam_]] = avg_
    }
  }
}

# el num de repeticiones no afecta la métrica del aveg silhouette
df_avg_silhouette2 = bind_rows(out_) %>% `colnames<-`(paste0('p_', names(.))) %>%
  mutate(k = 2:16) %>% 
  pivot_longer(cols = starts_with('p_'),, names_to = "params", values_to = "avg_")

ggplot(mutate(df_avg_silhouette2, 
                            params = str_replace(str_replace(params, '\\|100$',''), 'p_', ''),
                            n_tree = as.numeric(str_replace(params, pattern = "^[^|]*\\|",'')),
                            mtry1 = paste0('mtry1 = ',str_replace(params, pattern = "\\|.*",''))) %>% 
                       filter(n_tree>100) %>% 
                       filter(mtry1 %in% c('mtry1 = 3', 'mtry1 = 4')), 
                     aes(x = k, y = avg_, color = factor(n_tree))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = 'Number of clusters k',
    y = 'Average Silhouette width',
    color = "Parámetros \nn.rep=100 \nn.tree:"
  ) +
  theme_minimal() +
  facet_wrap(~ mtry1, scales = "fixed") +
  my_theme()+ 
  theme(strip.text = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(vjust=0, angle =0))

rfds <- RFdist(mtrx, 
               mtry1 = 4, # num de variables en c/división
               1000, # num de árboles
               100, # num de repeticiones/ num de bosques a promediar
               addcl1 = T, 
               addcl2 = F, 
               imp = T, 
               oob.prox1 = T)

lbrf_2 <- pamNew(rfds$cl1, 2)
lbrf_4 <- pamNew(rfds$cl1, 4)

df_pcp_monthly_normal2 = df_pcp_monthly_normal %>% ungroup() %>%  
  mutate(clust2_RF = factor(lbrf_2$label),
         clust4_RF = factor(lbrf_4$label))

df_RF_pivot_longer = df_pcp_monthly_normal2 %>% 
  pivot_longer(cols = 2:13, names_to = 'mes', values_to = 'pcp_monthly_mean') %>% 
  mutate(mes = factor(mes, levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))

df_RF_pcp_monthly_normal = st_as_sf(df_pcp_monthly_normal2, 
                                    coords=c('X', 'Y'), 
                                    crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

map2_RF <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_RF_pcp_monthly_normal, 'Cluster'=clust2_RF), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(17, 16))+
  scale_color_manual(values=c("#F8766D","#00BFC4"))+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9)) + 
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))

map4_RF <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_RF_pcp_monthly_normal, 'Cluster'=clust4_RF), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(17, 18, 16, 4))+
  scale_color_manual(values=c("#F8766D", "#A3A500","#00BFC4", "#C77CFF" ))+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))+ 
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))

bx2_RF <- ggplot(rename(df_RF_pivot_longer, 'Cluster'=clust2_RF), 
                 aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  labs(y = "Normal precipitación mensual") +
  my_theme_bx()

bx4_RF <- ggplot(rename(df_RF_pivot_longer, 'Cluster'=clust4_RF), 
                 aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#A3A500","#00BFC4", "#C77CFF" ))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()

grid.arrange(map2_RF, bx2_RF, ncol = 2)
grid.arrange(map4_RF, bx4_RF, ncol = 2)

join_cluster = left_join(df_pcp_annual_mean,
                         dplyr::select(df_pcp_monthly_normal2, ID, clust4_RF, clust2_RF),
                         by = c('ID'))

table(join_cluster$clust2, join_cluster$clust2_RF)
table(join_cluster$clust4, join_cluster$clust4_RF)

pct2_RF <- df_pcp_monthly_normal2 %>% 
  count(clust2_RF) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust2_RF)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface = "bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

pct4_RF <- df_pcp_monthly_normal2 %>% 
  count(clust4_RF) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust4_RF)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#F8766D", "#A3A500","#00BFC4", "#C77CFF" ))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface = "bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

grid.arrange(pct2_RF, pct4_RF, ncol = 2,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

## 2.2 addcl2: generación de v.a. uniformes en el rango observado por columna ----

nfrs_ = 100 #c(10, 50, 100)#seq(10, 100, 20)
ntrs_ = c(500,1000,2000) #seq(100, 1000, 200)
mtry1_ = c(3,4) #seq(2, 10, 2)

# mtry1: sqrt(p) -> sqrt(12) = 3.4
out_2 = list()

for(i in mtry1_){
  for(j in ntrs_){
    for(h in nfrs_){
      rfds_2 <- RFdist(mtrx, 
                       mtry1 = i, # num de variables en c/división
                       no.tree = j, # num de árboles
                       no.rep = h, # num de repeticiones/ num de bosques a promediar
                       addcl1 = F, 
                       addcl2 = T, 
                       imp = T, 
                       oob.prox1 = T)
      
      avg_ = c()
      for(k in 2:16){
        lbrf <- pamNew(rfds_2$cl2, k)
        avg_ = c(avg_, lbrf$avg_silhouette_width)
      }
      
      nam_ = paste(i,j,h, sep ='|')  
      out_2[[nam_]] = avg_
    }
  }
}

# el num de repeticiones no afecta la métrica del aveg silhouette
df_avg_silhouette3 = bind_rows(out_2) %>% `colnames<-`(paste0('p_', names(.))) %>%
  mutate(k = 2:16) %>% 
  pivot_longer(cols = starts_with('p_'),, names_to = "params", values_to = "avg_")

ggplot(mutate(df_avg_silhouette3, 
              params = str_replace(str_replace(params, '\\|100$',''), 'p_', ''),
              n_tree = as.numeric(str_replace(params, pattern = "^[^|]*\\|",'')),
              mtry1 = paste0('mtry1 = ',str_replace(params, pattern = "\\|.*",''))) %>% 
         filter(n_tree>100) %>% 
         filter(mtry1 %in% c('mtry1 = 3', 'mtry1 = 4')), 
       aes(x = k, y = avg_, color = factor(n_tree))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = 'Number of clusters k',
    y = 'Average Silhouette width',
    color = "Parámetros \nn.rep=100 \nn.tree:"
  ) +
  theme_minimal() +
  facet_wrap(~ mtry1, scales = "fixed") +
  my_theme()+ 
  theme(strip.text = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(vjust=0, angle =0))

rfds2 <- RFdist(mtrx, 
                mtry1 = 4, # num de variables en c/división
                2000, # num de árboles
                100, # num de repeticiones/ num de bosques a promediar
                addcl1 = F, 
                addcl2 = T, 
                imp = T, 
                oob.prox1 = T)

lbrf_2 <- pamNew(rfds2$cl2, 2)
lbrf_3 <- pamNew(rfds2$cl2, 3)

df_pcp_monthly_normal3 = df_pcp_monthly_normal %>% ungroup() %>%  
  mutate(clust2_RF = factor(lbrf_2$label),
         clust3_RF = factor(lbrf_3$label))

join_cluster2 = left_join(df_pcp_annual_mean,
                         dplyr::select(df_pcp_monthly_normal3, ID, clust3_RF, clust2_RF),
                         by = c('ID'))

df_RF_pivot_longer2 = df_pcp_monthly_normal3 %>% 
  pivot_longer(cols = 2:13, names_to = 'mes', values_to = 'pcp_monthly_mean') %>% 
  mutate(mes = factor(mes, levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))

df_RF_pcp_monthly_normal3 = st_as_sf(df_pcp_monthly_normal3, 
                                    coords=c('X', 'Y'), 
                                    crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


map2_RF2 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_RF_pcp_monthly_normal3, 'Cluster'=clust2_RF), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(17, 16))+
  scale_color_manual(values=c("#F8766D","#00BFC4"))+
  theme_minimal()+
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))

map3_RF2 <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_RF_pcp_monthly_normal3, 'Cluster'=clust3_RF), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(17, 16, 4))+
  scale_color_manual(values=c("#F8766D", "#00BFC4", "#C77CFF"))+
  theme_minimal()+
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))


bx2_RF2 <- ggplot(rename(df_RF_pivot_longer2, 'Cluster'=clust2_RF), 
                 aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()

bx3_RF2 <- ggplot(rename(df_RF_pivot_longer2, 'Cluster'=clust3_RF), 
                 aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "#C77CFF"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9))+
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()

grid.arrange(map2_RF2, bx2_RF2, ncol = 2)
grid.arrange(map3_RF2, bx3_RF2, ncol = 2)

pct2_RF2 <- df_pcp_monthly_normal3 %>% 
  count(clust2_RF) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust2_RF)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface = 'bold', size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

pct3_RF2 <- df_pcp_monthly_normal3 %>% 
  count(clust3_RF) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust3_RF)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "#C77CFF"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface = 'bold', size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

grid.arrange(pct2_RF2, pct3_RF2, ncol = 2,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

# 3. K-means normales mensuales -----------------------------------------------

# Wss
fnb1 <- fviz_nbclust(mtrx, FUN = kmeans, method = "wss") +
  geom_hline(yintercept = km_list[[5]]$tot.withinss, linetype = 2) +
  my_theme_bx()

# de 3  a 4 grupos
# Silhouette
fnb2 <- fviz_nbclust(mtrx, FUN = kmeans, 
                     method = "silhouette", k.max = 20)+
  my_theme() # 2 o 3 grupos

# gap_stat
# compute gap statistic
set.seed(12)
gap_stat <- clusGap(mtrx, FUN = kmeans, nstart = 25,
                    K.max = 25, B = 50)
print(gap_stat, method = "firstmax")

set.seed(12)
fnb3 <- fviz_nbclust(mtrx, FUN = kmeans, method = "gap_stat", k.max = 25) #+ 
  #geom_hline(yintercept = 0.81, linetype = 2) +
  #coord_cartesian(xlim = c(1.75,20)) # 2, 13 o 16 grupos

grid.arrange(fnb1, fnb2, ncol = 2)

n_optim = c(2,3,4,5)
km_list = list()

for(k in n_optim){
  km_list[[k]] = kmeans(mtrx, k, nstart = 25)
}

df_pcp_monthly_normal4 = df_pcp_monthly_normal %>% ungroup() %>% 
  mutate(clust2 = factor(km_list[[2]]$cluster),
         clust3 = factor(km_list[[3]]$cluster),
         clust4 = factor(km_list[[4]]$cluster),
         clust5 = factor(km_list[[5]]$cluster))

ggplot(df_pcp_monthly_normal4, aes(X, Y, color = clust2)) +
  geom_point(alpha = 0.25)

row.names(mtrx) <- df_pcp_monthly_normal$ID
fviz_cluster(km_list[[2]], data = mtrx)

pca_res <- prcomp(mtrx, scale = TRUE, center = TRUE)
loadings <- pca_res$rotation
eigenvalues <- pca_res$sdev^2
cor_var_PC <- sweep(loadings, 2, sqrt(eigenvalues), FUN = "*")
cor_var_PC[, 1:2]

contrib <- loadings^2
contrib <- sweep(contrib, 2, colSums(contrib), FUN = "/")
contrib <- contrib * 100
contrib_PC1_PC2 <- contrib[, 1:2]

PC_matrix <- pca_res$x
plot_data <- PC_matrix[, 1:2]
head(plot_data)

ggplot(data.frame(plot_data), aes(PC1, PC2, color = factor(km_list[[2]]$cluster))) +
  geom_point() +
  theme_minimal()

df_km_pivot_longer <- df_pcp_monthly_normal4 %>% 
  pivot_longer(cols = 2:13, names_to = 'mes', values_to = 'pcp_monthly_mean') %>% 
  mutate(mes = factor(mes, levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))

df_pcp_monthly_normal4 = st_as_sf(df_pcp_monthly_normal4, 
                                  coords = c('X', 'Y'),
                                  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


map2_km <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_pcp_monthly_normal4, 'Cluster'=clust2), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(16, 17))+
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+
  theme_minimal()+
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))

map5_km <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_pcp_monthly_normal4, 'Cluster'=clust5), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(4,17,15, 18, 16))+
  scale_color_manual(values=c("#C77CFF","#F8766D", "#FB61D7", "#A3A500","#00BFC4"))+
  theme_minimal()+
  my_theme()+
  theme(legend.position = c(0.9, 0.85))+
  guides(color = guide_legend(override.aes = list(size = 4)))

bx2_km <- ggplot(rename(df_km_pivot_longer, 'Cluster'=clust2), 
                 aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#00BFC4", "#F8766D"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()
  
bx5_km <- ggplot(rename(df_km_pivot_longer, 'Cluster'=clust5), 
                 aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#C77CFF","#F8766D", "#FB61D7", "#A3A500","#00BFC4"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal()+
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()+
  theme(legend.position = c(0.9, 0.85))

grid.arrange(map2_km, bx2_km, ncol = 2)
grid.arrange(map5_km, bx5_km, ncol = 2)

pct2_km2 <- df_pcp_monthly_normal4 %>% 
  count(clust2) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust2)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#00BFC4", "#F8766D"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface ="bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

pct5_km2 <- df_pcp_monthly_normal4 %>% 
  count(clust5) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust5)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#C77CFF","#F8766D", "#FB61D7", "#A3A500","#00BFC4"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface ="bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

grid.arrange(pct2_RF, pct4_RF, ncol = 2,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

grid.arrange(pct2_km2, pct5_km2, ncol = 2,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

names(df_pcp_monthly_normal2) # 2 o 5 grupos
names(df_pcp_annual_mean) # 2, 4 o 13 grupos

join_km = left_join(dplyr::select(df_pcp_annual_mean, ID, clust2, clust4) ,
                    dplyr::select(df_pcp_monthly_normal2, ID, clust2, clust5) %>% 
                      rename('clust2_km2'=clust2, 'clust5_km2'= clust5),
                    by ='ID') 

table(join_km$clust2, join_km$clust2_km2)
table(join_km$clust4, join_km$clust5_km2)

# 4. librería dtwclust ----

# series are z-normalized by means of the zscore function.
mtrx.norm <- zscore(mtrx)
##centroid= "pam": Partition around medoids (PAM). This basically means that 
## the cluster centroids are always one of the time series in the data.

## 4.1 Simple partitional clustering with Euclidean distance and PAM centroids -----
# sugiere 2 grupos
#distance = "L2", centroid = "pam", seed = 3247, control = partitional_control(nrep = 10L)

workers<-makeCluster(8L)

invisible(clusterEvalQ(workers,library("dtwclust")))
# Registerthebackend;thisstepMUSTbedone
registerDoParallel(workers)

clust.p_L2_pam <- tsclust(mtrx.norm, 
                          type="partitional", 
                          k=2L:13L, 
                          distance="L2", window.size = 20L,
                          centroid="pam",
                          seed = 42)#,
                          #control = partitional_control(nrep = 8L))

stopCluster(workers)
# Gobacktosequentialcomputation
registerDoSEQ()

## Cluster evaluation | Internal - Crisp partitions
## CVI            -   Minimized or Maximized   -   Considerations
## Silhouette     -   Maximized   -   Requires the whole cross-distance matrix.
## Score function -   Maximized  - Calculates a global centroid.
## Calinski-Harabasz - Maximized - Calculates a global centroid.
## Davies-Bouldin -   Minimized    - Calculates distances to the computed cluster centroids.
## Modified Davies-Bouldin (DB*) - Minimized - Calculates distances to the computed cluster centroids.
## Dunn           -   Maximized   -   Requires the whole cross-distance matrix.
## COP            -   Minimized   -   Requires the whole cross-distance matrix.

names(clust.p_L2_pam) <-paste0("k_", rep(2:13, each = 1))
cvi_p_L2_pam = sapply(clust.p_L2_pam, cvi, type ="internal")

options(scipen = 999)
cvi_min = function(M){
  ####################################################################
  # M: matrix
  ####################################################################
  df <- M[row.names(M) %in% c('DB', 'DBstar', 'COP'),] %>% 
    as.data.frame() %>% mutate(cvi = row.names(.)) %>% 
    pivot_longer(cols = starts_with('k_'), names_to = "n_clust", values_to = "index") %>% 
    mutate(n_clust = as.numeric(str_replace(n_clust, 'k_', '')))
  return(df)
}

cvi_max = function(M){
  ####################################################################
  # M: matrix
  ####################################################################
  df <- M[row.names(M) %in% c('Sil', 'SF', 'CH', 'D'),] %>% 
    as.data.frame() %>% mutate(cvi = row.names(.)) %>% 
    pivot_longer(cols = starts_with('k_'), names_to = "n_clust", values_to = "index") %>% 
    mutate(n_clust = as.numeric(str_replace(n_clust, 'k_', '')))
  return(df)
}

cvi_min_ = cvi_min(cvi_p_L2_pam)

min_cvi_ = ggplot(cvi_min_, aes(x = n_clust, y = index)) + 
  geom_line() + 
  facet_wrap(~ cvi, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = seq(2, 13, 1)) +
  theme_minimal() + 
  labs(y = '', x = 'Number of clusters k ')+ 
  theme(strip.text = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey", size = 0.25))

cvi_max_ = cvi_max(cvi_p_L2_pam)

View(cvi_max_)
cvi_max_ = cvi_p_L2_pam[row.names(cvi_p_L2_pam) %in% c('Sil', 'SF', 'CH', 'D'),] %>% 
  as.data.frame() %>% mutate(cvi = row.names(.)) %>% 
  pivot_longer(cols = starts_with('k_'), names_to = "n_clust", values_to = "index") %>% 
  mutate(n_clust = as.numeric(str_replace(n_clust, 'k_', ''))) #%>% 
  #pivot_wider(names_from = cvi, values_from = index)

max_cvi_ = ggplot(cvi_max_, aes(x = n_clust, y = index)) + 
  geom_line() +
  facet_wrap(~ cvi, scales = "free") +
  scale_x_continuous(breaks = seq(2, 13, 1)) +
  theme_minimal() +
  labs(y = '', x = 'Number of clusters k ')+ 
  theme(strip.text = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey", size = 0.25))

grid.arrange(min_cvi_, max_cvi_, ncol = 2,
             top = textGrob(
               "Índices de validación de clusters",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

# !!! No compila debido a las repeticiones
#names(clust.p_L2_pam2) <-paste0("k_", rep(2:5, each = 8))
#sapply(clust.p_L2_pam2, cvi, type ="internal")

workers<-makeCluster(8L)

invisible(clusterEvalQ(workers,library("dtwclust")))
# Registerthebackend;thisstepMUSTbedone
registerDoParallel(workers)

clust2.p_L2_pam <- tsclust(mtrx.norm, 
                          type="partitional", 
                          k=2L, 
                          distance="L2", window.size = 20L,
                          centroid="pam",
                          seed = 42,
                          control = partitional_control(nrep = 4L))
stopCluster(workers)
# Gobacktosequentialcomputation
registerDoSEQ()

names(clust2.p_L2_pam) <-paste0("r_", 1L:4L)
clust2.p_L2_pam_ <- cl_ensemble(list = clust2.p_L2_pam)
cl_dissimilarity(clust2.p_L2_pam_)

# !!! No funciona
#table(Medoid = cl_class_ids(cl_medoid(clust2.p_L2_pam)),
#      "True Classes"= rep(c(2L, 1L), each= 4L))

## 4.2 Simple partitional clustering with dtw distance and PAM centroids ----
## sugiere 2 o 5 grupos
# type="partitional", k=2L:17L, distance="dtw", centroid="pam"

workers<-makeCluster(8L)

invisible(clusterEvalQ(workers,library("dtwclust")))
# Registerthebackend;thisstepMUSTbedone
registerDoParallel(workers)

# con distance="dtw" In distfun(x = x, centroids = centroids, window.size = window.size) :
# Package 'bigmemory' is not available, cannot parallelize computation with 'dtw'. Use options(dtwclust_suggest_bigmemory = FALSE) to avoid this warning.
options(dtwclust_suggest_bigmemory = FALSE)

clust.p_dtw_pam <- tsclust(mtrx.norm, 
                           type="partitional", 
                           k=2L:13L, 
                           distance="dtw", window.size = 20L,
                           centroid="pam",
                           seed = 42,
                           control = partitional_control(nrep = 25L))

names(clust.p_dtw_pam) <-paste0("k_", rep(2:13, each = 25))

cvi_p_dtw_pam = sapply(clust.p_dtw_pam, cvi, type ="internal")

stopCluster(workers)
# Gobacktosequentialcomputation
registerDoSEQ()

cvi_min_2 = cvi_p_dtw_pam[row.names(cvi_p_dtw_pam) %in% c('DB', 'DBstar', 'COP'),] %>% 
  as.data.frame() %>% `colnames<-`(paste0(names(.), '_',1:ncol(.))) %>% 
  mutate(cvi = row.names(.)) %>%  `colnames<-`(str_replace(names(.), "_[^_]*$",'')) %>% 
  pivot_longer(cols = starts_with('k_'), names_to = "n_clust", values_to = "index") %>% 
  mutate(n_clust = as.numeric(str_replace(n_clust, 'k_', '')))

min_cvi_2 <- ggplot(cvi_min_2, aes(x = n_clust, y = index, group = n_clust)) + 
  geom_boxplot() + 
  facet_wrap(~ cvi, scales = "free", ncol = 2) + 
  scale_x_continuous(breaks = seq(2, 13, 1)) +
  theme_minimal()+ 
  labs(y = '', x = 'Number of clusters k ')+ 
  theme(strip.text = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey", size = 0.25))

cvi_max_2 = cvi_p_dtw_pam[row.names(cvi_p_dtw_pam) %in% c('Sil', 'SF', 'CH', 'D'),] %>% 
  as.data.frame() %>% `colnames<-`(paste0(names(.), '_',1:ncol(.))) %>% 
  mutate(cvi = row.names(.)) %>%  `colnames<-`(str_replace(names(.), "_[^_]*$",'')) %>% 
  pivot_longer(cols = starts_with('k_'), names_to = "n_clust", values_to = "index") %>% 
  mutate(n_clust = as.numeric(str_replace(n_clust, 'k_', '')))

group_by(cvi_max_2, cvi, n_clust) %>% summarise(quantile(index, probs = 0.5)) %>%  View()

max_cvi_2 <- ggplot(cvi_max_2, aes(x = n_clust, y = index, group = n_clust)) + 
  geom_boxplot() + 
  facet_wrap(~ cvi, scales = "free") +
  scale_x_continuous(breaks = seq(2, 13, 1)) +
  theme_minimal() + 
  labs(y = '', x = 'Number of clusters k ')+ 
  theme(strip.text = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey", size = 0.25))

grid.arrange(min_cvi_2, max_cvi_2, ncol = 2,
             top = textGrob(
               "Índices de validación de clusters",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

workers<-makeCluster(8L)

invisible(clusterEvalQ(workers,library("dtwclust")))
# Registerthebackend;thisstepMUSTbedone
registerDoParallel(workers)

clust4.p_dtw_pam <- tsclust(mtrx.norm, 
                            type="partitional", 
                            k=4L, 
                            distance="dtw", window.size = 20L,
                            centroid="pam",
                            seed = 42,
                            control = partitional_control(nrep = 25L))
stopCluster(workers)
# Gobacktosequentialcomputation
registerDoSEQ()

names(clust4.p_dtw_pam) <-paste0("r_", 1L:25L)
clust4.p_dtw_pam_ <- cl_ensemble(list = clust4.p_dtw_pam)
cl_dissimilarity(clust4.p_dtw_pam)

g1_p_dtw_pam <- plot(clust4.p_dtw_pam$r_8, type = "sc")
g2_p_dtw_pam <- plot(clust4.p_dtw_pam$r_8, type = "centroids")

grid.arrange(g1_p_dtw_pam, g2_p_dtw_pam, ncol = 1)

# !!! No funciona
#table(Medoid = cl_class_ids(cl_medoid(clust5.p_dtw_pam_)),
# "True Classes"= rep(c(5L, 4L, 3L, 2L, 1L), each= 4L))

## 4.3 Hierarchical clustering with Euclidean distance and PAM centroids -------
# type = "hierarchical", distance = "L2", centroid = "pam"

cvi_p_L2_pam_rep = sapply(clust2.p_L2_pam, cvi, type = 'internal')
cvi_mean_p_L2_pam = apply(cvi_p_L2_pam_rep, MARGIN = 1, mean) # promedio métricas

cols_max = c('Sil', 'SF', 'CH', 'D')
# por votación se queda con la réplica 1 de mejor desempeño
sapply(row.names(cvi_p_L2_pam_rep), function(x) {
  if(x %in% cols_max){
    which.max(cvi_p_L2_pam_rep[row.names(cvi_p_L2_pam_rep)==x,])
  }else{
    which.min(cvi_p_L2_pam_rep[row.names(cvi_p_L2_pam_rep)==x,])
  }
})

workers<-makeCluster(8L)

invisible(clusterEvalQ(workers,library("dtwclust")))
# Registerthebackend;thisstepMUSTbedone
registerDoParallel(workers)

hclust_methods<-c("single","complete", "average", "mcquitty")
clust.h_L2_pam <- tsclust(mtrx.norm, 
                          type="hierarchical", 
                          k=5L, 
                          #control=hierarchical_control(method = "average"))
                          control = hierarchical_control(method =  hclust_methods))
                          
stopCluster(workers)
# Gobacktosequentialcomputation
registerDoSEQ()

cvi_h_L2_pam = sapply(clust.h_L2_pam, cvi, type = 'internal')

col_min = c('DB', 'DBstar', 'COP')
cvi_min = cvi_h_L2_pam[row.names(cvi_h_L2_pam) %in% col_min,] %>% 
  as.data.frame() %>% `colnames<-`(paste0('h_',hclust_methods)) %>% 
  mutate(cvi = row.names(.)) %>% 
  pivot_longer(cols = starts_with('h_'), names_to = "h_method", values_to = "index") %>% 
  pivot_wider(names_from = cvi, values_from = index)

# average | average | single
# k = 2 single | single | single
# k = 5 complete | complete | single
sapply(col_min, function(x) {
  cvi_min$h_method[cvi_min[,x] == max(cvi_min[,x])]
})

cvi_max = cvi_h_L2_pam[row.names(cvi_h_L2_pam) %in% cols_max,] %>% 
  as.data.frame() %>% `colnames<-`(paste0('h_', hclust_methods)) %>% 
  mutate(cvi = row.names(.)) %>% 
  pivot_longer(cols = starts_with('h_'), names_to = "h_method", values_to = "index") %>% 
  pivot_wider(names_from = cvi, values_from = index)

# complete | complete | average | single
# k = 2 average | mcquitty | complete | single
# k= 5 mcquitty | single | complete | single
sapply(cols_max, function(x) {
  cvi_max$h_method[cvi_max[,x] == max(cvi_max[,x])]
})

# fit hierarchical average con misma distancia L2
# por mayor votación (3) se elige "average"
clust.h_cmplt_L2_pam <- tsclust(mtrx.norm, 
                                type="hierarchical", 
                                k=5L, 
                                control = hierarchical_control(method = "complete"))

# el método "single" da resultados desproporcionados, casi todos a un cluster
plot(clust.h_cmplt_L2_pam) # gráfico dendograma
table(cutree(clust.h_cmplt_L2_pam, k =5))

# Error >> por lo que cvi es de la partición anterior
#cvi_h_avg_L2_pam = sapply(clust.h_avg_L2_pam, cvi, type = 'internal')

# compara con los valores de Simple partitional
df_L2_pam = cvi_h_L2_pam %>% 
  as.data.frame() %>% `colnames<-`(paste0('h_',hclust_methods)) %>% 
  dplyr::select('h_complete') %>% 
  bind_cols(as.data.frame(cvi_mean_p_L2_pam))

View(df_L2_pam)
# por votación se queda con el promedio de las réplicas de partitional - L2 - pam
sapply(row.names(df_L2_pam), function(x) {
  if(x %in% cols_max){
    which.max(df_L2_pam[row.names(df_L2_pam)==x,])
  }else{
    which.min(df_L2_pam[row.names(df_L2_pam)==x,])
  }
})

# crea la columna con las etiquetas del clustering
idx_medoid = cl_class_ids(cl_medoid(clust2.p_L2_pam))
# se elige la réplica 3 por tener las mejores cvi
idx_r1 = cl_class_ids(clust2.p_L2_pam$r_1)

g1_p_l2_pam <- plot(clust2.p_L2_pam$r_1, type = "sc")
g2_p_l2_pam <- plot(clust2.p_L2_pam$r_1, type = "centroids")

grid.arrange(g1_p_l2_pam, g2_p_l2_pam, ncol = 1)

df_pcp_monthly_normal5 = df_pcp_monthly_normal %>% ungroup() %>% 
  mutate(clust2 = factor(cl_class_ids(clust2.p_L2_pam$r_1)),
         clust5 = factor(cl_class_ids(clust.h_cmplt_L2_pam)) )

df_p_l2_pam_pivot_longer <- df_pcp_monthly_normal5 %>% 
  pivot_longer(cols = 2:13, names_to = 'mes', values_to = 'pcp_monthly_mean') %>% 
  mutate(mes = factor(mes, levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))

df_pcp_monthly_normal5 = st_as_sf(df_pcp_monthly_normal5, 
                                  coords = c('X', 'Y'),
                                  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

map2_p_l2_pam <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_pcp_monthly_normal5, 'Cluster'=clust2), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(17, 16))+
  scale_color_manual(values=c("#F8766D","#00BFC4"))+
  theme_minimal()+
  my_theme()+
  guides(color = guide_legend(override.aes = list(size = 4)))

bx2_p_l2_pam <- ggplot(rename(df_p_l2_pam_pivot_longer, 'Cluster'=clust2), 
                       aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal() +
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()

grid.arrange(map2_p_l2_pam, bx2_p_l2_pam, ncol = 2)

map5_h_l2_pam <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_pcp_monthly_normal5, 'Cluster'=clust5), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(4,17,18, 15, 16))+
  scale_color_manual(values=c("#C77CFF","#F8766D", "#A3A500", "#FB61D7", "#00BFC4"))+
  theme_minimal()+
  my_theme()+
  theme(legend.position = c(0.9, 0.85))+
  guides(color = guide_legend(override.aes = list(size = 4)))

bx5_h_l2_pam <- ggplot(rename(df_p_l2_pam_pivot_longer, 'Cluster'=clust5), 
                       aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#C77CFF","#F8766D", "#A3A500","#FB61D7", "#00BFC4"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal() +
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()+
  theme(legend.position = c(0.9, 0.85))

grid.arrange(map5_h_l2_pam, bx5_h_l2_pam, ncol = 2)

g1_h_l2_pam <- plot(clust.h_cmplt_L2_pam, type = "sc")
g2_h_l2_pam <- plot(clust.h_cmplt_L2_pam, type = "centroids")

grid.arrange(g1_h_l2_pam, g2_h_l2_pam, ncol = 1)

pct2_p_L2_pam <- df_pcp_monthly_normal5 %>% 
  count(clust2) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust2)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface ="bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

pct5_h_L2_pam <- df_pcp_monthly_normal5 %>% 
  count(clust5) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust5)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#C77CFF","#F8766D", "#A3A500","#FB61D7", "#00BFC4"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface ="bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

grid.arrange(pct2_p_L2_pam, pct5_h_L2_pam, ncol = 2,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

## 4.4 Hierarchical clustering with dtw distance -----------------------------
# type = "h", k = 6L, distance = "dtw"

cvi_p_dtw_pam_rep = sapply(clust4.p_dtw_pam, cvi, type = 'internal')
cvi_mean_p_dtw_pam = apply(cvi_p_dtw_pam_rep, MARGIN = 1, mean) # promedio métricas

# por votación se queda con la réplica 8 de mejor desempeño
sapply(row.names(cvi_p_dtw_pam_rep), function(x) {
  if(x %in% cols_max){
    which.max(cvi_p_dtw_pam_rep[row.names(cvi_p_dtw_pam_rep)==x,])
  }else{
    which.min(cvi_p_dtw_pam_rep[row.names(cvi_p_dtw_pam_rep)==x,])
  }
})

workers<-makeCluster(8L)

invisible(clusterEvalQ(workers,library("dtwclust")))
# Registerthebackend;thisstepMUSTbedone
registerDoParallel(workers)

clust.h_dtw_pam <- tsclust(mtrx.norm, 
                           type="hierarchical",
                           k = 4L, 
                           distance = "dtw",
                           control = hierarchical_control(method =  hclust_methods))

stopCluster(workers)
# Gobacktosequentialcomputation
registerDoSEQ()

cvi_h_dtw_pam = sapply(clust.h_dtw_pam, cvi, type = 'internal')

cvi_min = cvi_h_dtw_pam[row.names(cvi_h_dtw_pam) %in% col_min,] %>% 
  as.data.frame() %>% `colnames<-`(paste0('h_',hclust_methods)) %>% 
  mutate(cvi = row.names(.)) %>% 
  pivot_longer(cols = starts_with('h_'), names_to = "h_method", values_to = "index") %>% 
  pivot_wider(names_from = cvi, values_from = index)

# complete | complete | single
sapply(col_min, function(x) {
  cvi_min$h_method[cvi_min[,x] == max(cvi_min[,x])]
})

cvi_max = cvi_h_dtw_pam[row.names(cvi_h_dtw_pam) %in% cols_max,] %>% 
  as.data.frame() %>% `colnames<-`(paste0('h_', hclust_methods)) %>% 
  mutate(cvi = row.names(.)) %>% 
  pivot_longer(cols = starts_with('h_'), names_to = "h_method", values_to = "index") %>% 
  pivot_wider(names_from = cvi, values_from = index)

# average | single | average | single
sapply(cols_max, function(x) {
  cvi_max$h_method[cvi_max[,x] == max(cvi_max[,x])]
})

# fit hierarchical average con misma distancia dtw
# por mayor votación (3) se elige "complete"
clust.h_cmplt_dtw_pam <- tsclust(mtrx.norm, 
                                type="hierarchical", 
                                k=4L, 
                                distance = 'dtw',
                                control = hierarchical_control(method = "complete"))

# el método "single" da resultados desproporcionados, casi todos a un cluster
plot(clust.h_cmplt_dtw_pam) # gráfico dendograma
table(cutree(clust.h_cmplt_dtw_pam, k =4))

# Error >> por lo que cvi es de la partición anterior
#cvi_h_cmplt_L2_pam = sapply(clust.h_cmplt_L2_pam, cvi, type = 'internal')

# compara con los valores de Simple partitional
df_dtw_pam = cvi_h_dtw_pam %>% 
  as.data.frame() %>% `colnames<-`(paste0('h_',hclust_methods)) %>% 
  dplyr::select('h_complete') %>% 
  bind_cols(as.data.frame(cvi_mean_p_dtw_pam))

# por votación se queda con el promedio de las réplicas de partitional - dtw - pam
sapply(row.names(df_dtw_pam), function(x) {
  if(x %in% cols_max){
    which.max(df_dtw_pam[row.names(df_dtw_pam)==x,])
  }else{
    which.min(df_dtw_pam[row.names(df_dtw_pam)==x,])
  }
})

# crea la columna con las etiquetas del clustering
idx_medoid = cl_class_ids(cl_medoid(clust4.p_dtw_pam))
# se elige la réplica 3 por tener las mejores cvi
idx_r8 = cl_class_ids(clust4.p_dtw_pam$r_8)


plot(clust4.p_dtw_pam$r_8, type = "sc")
plot(clust4.p_dtw_pam$r_8, type = "centroids")

df_pcp_monthly_normal6 = df_pcp_monthly_normal %>% ungroup() %>% 
  mutate(clust4 = factor(cl_class_ids(clust4.p_dtw_pam$r_8)),
         clust4_h = factor(cl_class_ids(clust.h_cmplt_dtw_pam)))

df_p_dtw_pam_pivot_longer <- df_pcp_monthly_normal6 %>% 
  pivot_longer(cols = 2:13, names_to = 'mes', values_to = 'pcp_monthly_mean') %>% 
  mutate(mes = factor(mes, levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))

df_pcp_monthly_normal6 = st_as_sf(df_pcp_monthly_normal6, 
                                  coords = c('X', 'Y'),
                                  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

map4_p_dtw <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_pcp_monthly_normal6, 'Cluster'=clust4), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(4,17,16, 18))+
  scale_color_manual(values=c("#C77CFF","#F8766D",  "#00BFC4","#A3A500"))+
  theme_minimal()+
  my_theme()+
  theme(legend.position = c(0.9, 0.85))+
  guides(color = guide_legend(override.aes = list(size = 4)))

map4_h_dtw <- ggplot() + 
  geom_sf(data = shp_depto, fill = "white", color = "gray") +  # Muestra el mapa de los departamentos
  geom_sf(data = rename(df_pcp_monthly_normal6, 'Cluster'=clust4_h), 
          aes(fill = Cluster, pch = Cluster, color = Cluster), size = 1) +
  scale_shape_manual(values=c(16,17,4, 18))+
  scale_color_manual(values=c("#00BFC4","#F8766D", "#C77CFF","#A3A500"))+
  theme_minimal()+
  my_theme()+
  theme(legend.position = c(0.9, 0.85))+
  guides(color = guide_legend(override.aes = list(size = 4)))


bx4_p_dtw <- ggplot(rename(df_p_dtw_pam_pivot_longer, 'Cluster'=clust4), 
                    aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#C77CFF","#F8766D", "#00BFC4", "#A3A500"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal() +
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()+
  theme(legend.position = c(0.9, 0.85))

bx4_h_dtw <- ggplot(rename(df_p_dtw_pam_pivot_longer, 'Cluster'=clust4_h), 
                    aes(x = mes, y = pcp_monthly_mean, fill = Cluster)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#00BFC4","#F8766D","#C77CFF" , "#A3A500"))+
  scale_y_continuous(breaks = seq(0, 1200, 400)) +
  theme_minimal() +
  labs(y = "Normal precipitación mensual")+
  my_theme_bx()+
  theme(legend.position = c(0.9, 0.85))


grid.arrange(map4_p_dtw, bx4_p_dtw, ncol = 2)
grid.arrange(map4_h_dtw, bx4_h_dtw, ncol = 2)

g1_h_dtw_pam <- plot(clust.h_cmplt_dtw_pam, type = "sc")
g2_h_dtw_pam <- plot(clust.h_cmplt_dtw_pam, type = "centroids")

grid.arrange(g1_h_dtw_pam, g2_h_dtw_pam, ncol = 1)

pct4_p_dtw_pam <- df_pcp_monthly_normal6 %>% 
  count(clust4) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust4)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#C77CFF","#F8766D", "#00BFC4", "#A3A500"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface ="bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))

pct4_h_dtw_pam <- df_pcp_monthly_normal6 %>% 
  count(clust4_h) %>%
  mutate(pct = n / sum(n) * 100,
         label = ifelse(pct>1, paste0(round(pct, 1), "%"), "")) %>%
  ggplot(aes(x = "", y = n, fill = clust4_h)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#00BFC4","#F8766D","#C77CFF" , "#A3A500"))+
  geom_text(aes(label = label),
            position = position_stack(vjust =0.5),
            fontface ="bold", size = 5, color = "white") +
  labs(fill = "Cluster") +
  #scale_fill_viridis_d(option = "turbo") + 
  theme_void() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))


grid.arrange(pct4_p_dtw_pam, pct4_h_dtw_pam, ncol = 2,
             top = textGrob(
               "Porcentaje de estaciones metereológicas en cada cluster",
               gp = gpar(fontsize = 24, fontface = "bold")
             ))

save.image('clust_sttns_dtwclust.RData')
#load('clust_sttns_dtwclust.RData')

# 5. librería HiClimR ----

# 6. Hierachical clustering of spatially correlated fd -------------------------

data(CanadianWeather, package="fda")
str(CanadianWeather, max.level = 1)

## 6.1 Coordenadas -----------------------------------------------------------

#coords <- spTransform(coords,canada.CRS)
coords_df <- data.frame(lon = -CanadianWeather$coordinates[, "W.longitude"], 
                        lat = CanadianWeather$coordinates[, "N.latitude"])

convertlatlong2UTM <- function(area, units = 'm') {
  # temporary sf conversion
  area <- sf::st_as_sf(area)
  bounds <- sf::st_bbox(area)
  lat = mean(bounds[c(2, 4)]) # latitude
  long = mean(bounds[c(1, 3)]) # longitude
  # find UTM hemisphere (latitude)
  hemisphere <- ifelse(lat > 0, "north", "south")
  # find UTM zone
  zone <- (floor((long + 180) / 6) %% 60) + 1
  crs <- paste0("+proj=utm +zone=", zone, " +datum=WGS84 +ellps=WGS84 +", hemisphere,
                " +units=", units, " +no_defs")
  return(crs)
}

# 0. fda coordinates
coords0 <- st_as_sf(coords_df, coords = c("lon", "lat"), crs = 4326)  
XY_m0 <- st_coordinates(coords0)   
X_m0 <- XY_m0[,1]; Y_m0 <- XY_m0[,2]
X_km0 <- X_m0/1000; Y_km0 <- Y_m0/1000
XY_km0 <- matrix(c(X_km0, Y_km0), nc = 2)

# 1. github.com/mpbohorquezc/SpatFD-Functional-Geostatistics/man/sim_functional_process.Rd
# NO convierte a negativo la longitud
canada.CRS <- CRS("+init=epsg:4608")
coords1 <- SpatialPoints(CanadianWeather$coordinates,
                         proj4string = CRS("+init=epsg:4326"))

coords1 <- spTransform(coords1, canada.CRS) %>% st_as_sf()
XY_m1 <- st_coordinates(coords1)   
X_m1 <- XY_m1[,1]; Y_m1 <- XY_m1[,2]
X_km1 <- X_m1/1000; Y_km1 <- Y_m1/1000
XY_km1 <- matrix(c(X_km1, Y_km1), nc = 2)

# 2. chatGPT para Canada (epsg:3978)
coords_sf <- st_as_sf(coords_df, coords = c("lon", "lat"), crs = 4326)  # WGS84
coords2 <- st_transform(coords_sf, 3978) 
XY_m2 <- st_coordinates(coords2)   
X_m2 <- XY_m2[,1]; Y_m2 <- XY_m2[,2]
X_km2 <- X_m2/1000; Y_km2 <- Y_m2/1000
XY_km2 <- matrix(c(X_km2, Y_km2), nc = 2)

# 3. convertlatlong2UTM: coordenadas planas
CanadianWeather_planar = coords_df
CanadianWeather$coordinates[,2] = -CanadianWeather$coordinates[,2]

for(i in 1:nrow(CanadianWeather$coordinates)){
  CanadianWeather_planar[i,] = as.data.frame(CanadianWeather$coordinates)[i,] %>%
    st_as_sf(.,coords=c("W.longitude", "N.latitude"), crs = 4326) %>%
    st_transform(crs = convertlatlong2UTM(.)) %>%
    st_coordinates()}

#CanadianWeather_planar = as.data.frame(CanadianWeather_planar) %>% 
#  `colnames<-`(c("W.longitude", "N.latitude"))

XY_m3 <- CanadianWeather_planar
X_m3 <- XY_m3[,1]; Y_m3 <- XY_m3[,2]
X_km3 <- X_m3/1000; Y_km3 <- Y_m3/1000
XY_km3 <- matrix(c(X_km3, Y_km3), nc = 2)

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


#coord = as.data.frame(CanadianWeather$coordinates)
#fRegress(CanadianWeather$coordinates)
#coords = as.data.frame(coords)

## 6.2 Remover la tendencia espacial con el modelo de regresión funcional ------
## X_i(t) = \alpha(t) + \beta_1(t) Longitude_i + \beta_2(t) Latitude_i + e_i (t)

#TempRgn.f <- fRegress(Temp.fd ~ N.latitude + W.longitude , coord)
TempRgn.fm0 <- fRegress(Temp.fd ~ X_m0 + X_m0)
TempRgn.fkm0 <- fRegress(Temp.fd ~ X_km0 + X_km0)

TempRgn.fm1 <- fRegress(Temp.fd  ~ X_m1 + X_m1)
TempRgn.fkm1 <- fRegress(Temp.fd  ~ X_km1 + X_km1)

TempRgn.fm2 <- fRegress(Temp.fd ~ X_m2 + X_m2) # !!! In eigchk(Cmat) : Near singularity in coefficient matrix.
TempRgn.fkm2 <- fRegress(Temp.fd ~ X_km2 + X_km2)

TempRgn.fm3 <- fRegress(Temp.fd ~ X_m3 + X_m3) # !!! In eigchk(Cmat) : Near singularity in coefficient matrix.
TempRgn.fkm3 <- fRegress(Temp.fd ~ X_km3 + X_km3)
#TempRgn.f <- fRegress(Temp.fd ~ coords.x1 + coords.x2 , coords)
#TempRgn.f2 <- fRegress(Temp.fd ~ N.latitude + W.longitude , CanadianWeather_planar)

## 6.3 Obtiene los residuales del modelo funcional -----------------------------

#fdobj.res = TempRgn.f$yfdobj-TempRgn.f$yhatfdobj
#fdobj.res = Temp.fd-TempRgn.f$yhatfdobj
fdobj.res_m0 = Temp.fd-TempRgn.fm0$yhatfdobj
fdobj.res_km0 = Temp.fd-TempRgn.fkm0$yhatfdobj

fdobj.res_m1 = Temp.fd-TempRgn.fm1$yhatfdobj
fdobj.res_km1 = Temp.fd-TempRgn.fkm1$yhatfdobj

fdobj.res_m2 = Temp.fd-TempRgn.fm2$yhatfdobj
fdobj.res_km2 = Temp.fd-TempRgn.fkm2$yhatfdobj

fdobj.res_m3 = Temp.fd-TempRgn.fm3$yhatfdobj
fdobj.res_km3 = Temp.fd-TempRgn.fkm3$yhatfdobj

plot(fdobj.res_m0)
#summary(TempRgn.f)

# evalua los residuales
day_grid <- 1:365
res_m0 <- eval.fd(day_grid, fdobj.res_m0)
res_km0 <- eval.fd(day_grid, fdobj.res_km0)

res_m1 <- eval.fd(day_grid, fdobj.res_m1)
res_km1 <- eval.fd(day_grid, fdobj.res_km1)

res_m2 <- eval.fd(day_grid, fdobj.res_m2)
res_km2 <- eval.fd(day_grid, fdobj.res_km2)

res_m3 <- eval.fd(day_grid, fdobj.res_m3)
res_km3 <- eval.fd(day_grid, fdobj.res_km3)

## 6.4 Ajusta la función okfd -------------------------------------------------
# Fit a spherical model to the estimated trace-variogram by using the OLS technique
coords.cero <- data.frame(Lon = -64.06, Lat = 44.79)

okfd.res_m0 <- okfd(new.coords = coords.cero, coords=XY_m0,
                    data=res_m0, 
                    smooth.type='fourier', nbasis=65, argvals=day.5,
                    fix.nugget=FALSE, fix.kappa=FALSE)

okfd.res_m0$trace.vari.array[[1]] # spherical: 40751.21 y 25.44

okfd.res_km0 <- okfd(new.coords = coords.cero, coords=XY_km0,
                    data=res_km0, 
                    smooth.type='fourier', nbasis=65, argvals=day.5,
                    fix.nugget=FALSE, fix.kappa=FALSE)

okfd.res_km0$trace.vari.array[[1]] # spherical: 4.537189e+04 y 2.873431e-02

okfd.res_m1 <- okfd(new.coords = coords.cero, coords=XY_m1,
                    data=res_m1, 
                    smooth.type='fourier', nbasis=65, argvals=day.5,
                    fix.nugget=TRUE, nugget = 0, fix.kappa=FALSE)

okfd.res_m1$trace.vari.array[[1]] # spherical: 16533.4686  y 65.1114

okfd.res_km1 <- okfd(new.coords = coords.cero, coords=XY_km1,
                    data=res_km1, 
                    smooth.type='fourier', nbasis=65, argvals=day.5,
                    fix.nugget=TRUE, nugget = 0, fix.kappa=FALSE)

okfd.res_km1$trace.vari.array[[1]] # spherical: 1.860350e+04 y 6.161156e-02

okfd.res_m2 <- okfd(new.coords = coords.cero, coords=XY_m2,
                    data=res_m2, 
                    smooth.type='fourier', nbasis=65, argvals=day.5,
                    fix.nugget=FALSE, fix.kappa=FALSE)

okfd.res_m2$trace.vari.array[[1]] # spherical: 43814 y 2674732

okfd.res_km2 <- okfd(new.coords = coords.cero, coords=XY_km2,
                     data=res_km2, 
                     smooth.type='fourier', nbasis=65, argvals=day.5,
                     fix.nugget=TRUE, nugget = 0, fix.kappa=FALSE)

okfd.res_km2$trace.vari.array[[1]] # spherical: 1.860350e+04 y 6.161156e-02

okfd.res_m3 <- okfd(new.coords = coords.cero, coords=XY_m3,
                    data=res_m3, 
                    smooth.type='fourier', nbasis=65, argvals=day.5,
                    fix.nugget=FALSE, fix.kappa=FALSE)

okfd.res_m3$trace.vari.array[[1]] # spherical: 43814 y 2674732

okfd.res_km3 <- okfd(new.coords = coords.cero, coords=XY_km3,
                    data=res_km3, 
                    smooth.type='fourier', nbasis=65, argvals=day.5,
                    fix.nugget=FALSE, fix.kappa=FALSE)

okfd.res_km3$trace.vari.array[[1]] # spherical: 43814 y 2674732


okfd.res$trace.vari.array[[1]] # spherical: 8034.03 y 22.08
okfd.res$trace.vari.array[[4]] # 8034.03 y 22.08

plot(okfd.res)

## 6.5 Ajusta la función fit.tracevariog --------------------------------------
M0 <- fourierpen(fdobj.res_m0$basis,  Lfdobj=0)
res.fd_m0 <- smooth.basis(argvals = 1:365, y = res_m0, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_m0$coefs), res.fd_m0, M0)

new.emp.trace.vari_m0 <- trace.variog(coords=XY_m0,
                                      L2norm=L2norm, bin=FALSE)

fit_m0 = geofd::fit.tracevariog(new.emp.trace.vari_m0, models = "spherical",
                                sigma2.0 = 15000, phi.0 = 250,
                                fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                max.dist.variogram=NULL)

fit_m0$best$cov.pars # 40750.55, 25.44

kM0 <- fourierpen(fdobj.res_km0$basis,  Lfdobj=0)
res.fd_km0 <- smooth.basis(argvals = 1:365, y = res_km0, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_km0$coefs), res.fd_km0, kM0)

new.emp.trace.vari_km0 <- trace.variog(coords=XY_km0,
                                       L2norm=L2norm, bin=FALSE)

fit_km0 = geofd::fit.tracevariog(new.emp.trace.vari_km0, models = "spherical",
                                sigma2.0 = 15000, phi.0 = 250,
                                fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                max.dist.variogram=NULL)

fit_km0$best$cov.pars # 2.356310e+04 y 1.807513e-02

M1 <- fourierpen(fdobj.res_m1$basis,  Lfdobj=0)
res.fd_m1 <- smooth.basis(argvals = 1:365, y = res_m1, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_m1$coefs), res.fd_m1, M1)

new.emp.trace.vari_m1 <- trace.variog(coords=XY_m1,
                                      L2norm=L2norm, bin=FALSE)

fit_m1 = geofd::fit.tracevariog(new.emp.trace.vari_m1, models = "spherical",
                                sigma2.0 = 16500, phi.0 = 66,
                                fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                max.dist.variogram=NULL)

fit_m1$best$cov.pars # 8033, 22.08

kM1 <- fourierpen(fdobj.res_km1$basis,  Lfdobj=0)
res.fd_km1 <- smooth.basis(argvals = 1:365, y = res_km1, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_km1$coefs), res.fd_km1, kM1)

new.emp.trace.vari_km1 <- trace.variog(coords=XY_km1,
                                       L2norm=L2norm, bin=FALSE)

fit_km1 = geofd::fit.tracevariog(new.emp.trace.vari_km1, models = "spherical",
                                 sigma2.0 = 15000, phi.0 = 250,
                                 fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                 max.dist.variogram=NULL)

fit_km1$best$cov.pars # 1.835429e+04  y 6.079655e-02

M2 <- fourierpen(fdobj.res_m2$basis,  Lfdobj=0)
res.fd_m2 <- smooth.basis(argvals = 1:365, y = res_m2, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_m2$coefs), res.fd_m2, M2)

new.emp.trace.vari_m2 <- trace.variog(coords=XY_m2,
                                      L2norm=L2norm, bin=FALSE)

fit_m2 = geofd::fit.tracevariog(new.emp.trace.vari_m2, models = "spherical",
                                sigma2.0 = 16500, phi.0 = 66,
                                fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                max.dist.variogram=NULL)

fit_m2$best$cov.pars # 8033, 22.08

kM2 <- fourierpen(fdobj.res_km2$basis,  Lfdobj=0)
res.fd_km2 <- smooth.basis(argvals = 1:365, y = res_km2, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_km2$coefs), res.fd_km2, kM2)

new.emp.trace.vari_km2 <- trace.variog(coords=XY_km2,
                                       L2norm=L2norm, bin=FALSE)

fit_km2 = geofd::fit.tracevariog(new.emp.trace.vari_km2, models = "spherical",
                                 sigma2.0 = 15000, phi.0 = 250,
                                 fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                 max.dist.variogram=NULL)

fit_km2$best$cov.pars # 1.835429e+04  y 6.079655e-02

M3 <- fourierpen(fdobj.res_m3$basis,  Lfdobj=0)
res.fd_m3 <- smooth.basis(argvals = 1:365, y = res_m3, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_m3$coefs), res.fd_m3, M3)

new.emp.trace.vari_m3 <- trace.variog(coords=XY_m3,
                                      L2norm=L2norm, bin=FALSE)

fit_m3 = geofd::fit.tracevariog(new.emp.trace.vari_m3, models = "spherical",
                                sigma2.0 = 16500, phi.0 = 100,
                                fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                max.dist.variogram=NULL)

fit_m3$best$cov.pars # 8033, 22.08

kM3 <- fourierpen(fdobj.res_km3$basis,  Lfdobj=0)
res.fd_km3 <- smooth.basis(argvals = 1:365, y = res_km3, fdParobj = Temp.basis)$fd
L2norm = l2.norm(ncol(fdobj.res_km3$coefs), res.fd_km3, kM3)

new.emp.trace.vari_km3 <- trace.variog(coords=XY_km3,
                                       L2norm=L2norm, bin=FALSE)

fit_km3 = geofd::fit.tracevariog(new.emp.trace.vari_km3, models = "spherical",
                                 sigma2.0 = 15000, phi.0 = 250,
                                 fix.nugget=TRUE, nugget=0, fix.kappa=FALSE,
                                 max.dist.variogram=NULL)

fit_km3$best$cov.pars # 4371399.4  y 165403.4


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
