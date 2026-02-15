
library(tidyverse)
library(astsa)
library(forecast)
library(geobr)
library(ggspatial)
library(urca)
library(readxl)
library(sf)
library(ggplot2)
library(maps)
library(sp)
library(units) # convertir a km
library(forcats)
library(gridExtra) # grid.arrange
library(xtable)
library(fda)

setwd('~/TS_climatic/')
# Nonstationary seasonal DFM estimation of Nieto et al. (2016)

pol2_list <- readRDS('sttns_validation_chirps.RDS')
montlyAvg_pcp_mm2_ = readRDS("df_mean_pcp_sttns_clust_wtv.RDS")

cusum_g <- function(res){
  #res=residuals_f2_m1
  cum=cumsum(res)/sd(res)
  N=length(res)
  cumq=cumsum(res^2)/sum(res^2)
  Af=0.948 ###Cuantil del 95% para la estad?stica cusum
  co=0.14422####Valor del cuantil aproximado para cusumsq para n/2
  LS=Af*sqrt(N)+2*Af*c(1:length(res))/sqrt(N)
  LI=-LS
  LQS=co+(1:length(res))/N
  LQI=-co+(1:length(res))/N
  
  par(mfrow=c(2,1))
  plot(cum,type="l",ylim=c(min(LI),max(LS)), lwd = 2,
       xaxt = "n", yaxt = "n", 
       xlab="",ylab="",main="CUSUM")
  axis(1, cex.axis = 2)
  axis(2, cex.axis = 2)
  mtext("t", side = 1, line = 3, cex = 2)
  lines(LS,type="S",col="red", lwd = 2)
  lines(LI,type="S",col="red", lwd = 2)
  #CUSUM Square
  plot(cumq,type="l", lwd = 2,
       xaxt = "n", yaxt = "n",
       xlab="",ylab="",main="CUSUMSQ")
  axis(1, cex.axis = 2)
  axis(2, cex.axis = 2)
  mtext("t", side = 1, line = 3, cex = 2)
  lines(LQS,type="S",col="red", lwd = 2)
  lines(LQI,type="S",col="red", lwd = 2)
}

# 1. Data ----------------------------------------------------------------------
id_aero = c(16015010, 26155110, 27015330, 23195502, 21245040, 21201230)
filter(montlyAvg_pcp_mm2_, ID %in% id_aero) %>%  count(clust_wtv2)

Id_clust4 = filter(montlyAvg_pcp_mm2_, clust_wtv2 == 4) %>% distinct(ID)

list_pcp = list()
k = 1
for(i in 1:length(pol2_list)){
  list_pcp[[k]] = dplyr::select(mutate(pol2_list[[i]], Date = as.Date(Date)), ID, Date, sttns_fill, NombreEstacion, DEPARTAMENTO) 
  k = k+1
}

df_clust4 = bind_rows(list_pcp) %>% filter(ID %in% Id_clust4$ID)

df_clust4_2 = df_clust4 %>%  filter(Date < as.Date('2013-07-01'))

Data <- df_clust4_2 %>%
  dplyr::select(NombreEstacion, Date, sttns_fill) %>%
  spread(key = NombreEstacion, value = sttns_fill)

names(Data)
X <- as.matrix(Data[,2:ncol(Data)]) # elimina fila del mes

# Scaling data
sX <- scale(X) # escala los datos (x_{ij} - \bar{x}_{j})

Data2 <- df_clust4 %>%
  dplyr::select(NombreEstacion, Date, sttns_fill) %>%
  spread(key = NombreEstacion, value = sttns_fill)

X2 <- as.matrix(Data2[,2:ncol(Data2)]) # elimina fila del mes

sX2 <- scale(X2) # escala los datos (x_{ij} - \bar{x}_{j})

View(head(sX2))

# 2. Step 1: Finding the number and type of factors --------------------------

S <- 12
nObs <- nrow(sX) # 204 meses

EVMatrix <- t(sapply(0:36, FUN = function(h){
  
  # rezagos
  Xt <- sX[(h+1):nObs,]
  Xt_h <- sX[1:(nObs-h),]
  
  # Sample generalized autocovariance of lag h
  
  Ch <- ((S/nObs)^2) * (t(Xt) %*% Xt_h)
  
  
  # Five largest absolute eigenvalues
  
  abs(eigen(Ch)$values)[1:6]
}))

# Eigenvalues figure
EVData <- data.frame(Lag = 0:36)

# 222 filas (37 lags * 6 EV), 3 cols
# Valor C_h combinación lag-EV
EVData <- cbind(EVData, EVMatrix) %>%
  gather(key = "EV", value = "Value", 2:7)

my_theme_sil <- function() {
  theme(
    strip.text = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24), 
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "lightgrey", size = 0.25),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24)
  )
}

ggplot(filter(EVData, Lag >0), aes(x = Lag, y = Value, colour = EV)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 36, by = 3)) +
  scale_color_manual(values = c("darkblue", "darkgreen", "darkred",
                                "darkorange", "darkgrey", "darkcyan" ),
                     labels = c(expression(lambda[1]), expression(lambda[2]), 
                                expression(lambda[3]), expression(lambda[4]),
                                expression(lambda[5]), expression(lambda[6]))) +
  labs(y = "Valor absoluto", x = "Lag k", colour = expression("Valores propios \n de C(k)")) +
  theme_bw() + my_theme_sil()
  #theme(axis.title = element_text(size = 12),
  #      legend.title = element_text(size = 12),
  #      axis.text = element_text(size = 11),
  #      legend.text = element_text(size = 11))

length(Data$Date)/ncol(Data) #< 20

source("SLBDD_DFMPC_test_nfactors.R")

rmax <- 36 #round(ncol(X)/2,digits = 0)
k0 = 3 # default lagk = 0 Maximum number of lags considered in the combined matrix.

LamYaotest(sX, k0, rmax)$r1 # 4 # Lam y Yao (2012)
AhnHortest(sX, rmax)$r1 #2 # Ahn y Horenstein (2013)
CaroPenatest(sX, k0, rmax)$r1 #2 # Caro y Peña (2020)
#Chi2Test.eps(sX, k =rmax)$r # # Bolívar et al. (2021).


# 3. Step 2: Finding a model for the factors ----------------------------------

S <- 12
nObs <- nrow(sX) #204 meses
h <- 12

(nObs/S)

Xt <- sX[(h+1):nObs,]
Xt_h <- sX[1:(nObs-h),]

# matriz de autocovarianza muestral generalizada para rezago h
Ch <- ((S/nObs)^2) * (t(Xt) %*% Xt_h)

print(range(Ch)) # -0.1788231  0.2522219 # 0.2548381 0.6055674


Load <- Re(eigen(Ch)$vectors[,1:4])
dim(Load) # 352   4
# Factors
Fac <- sX %*% Load
# 390 meses
# num [1:390, 1:4] -15.38 -5.38 -5.14 29.48 33.92 ...

# Common components
Chi <- Fac %*% t(Load)
print(dim(Chi)) # 390 meses x 352 sttns
colnames(Chi) <- colnames(X)
View(head(Chi[,1:6]))

t(Load) %*% Load

Load_ = Load %*% solve(chol(t(Load) %*% Load))
t(Load_) %*%Load_

# the matrix C = PH, where H = [H(l, k)] is of dimension r × r^* and 
# its entries are given in the following way.
# For the ith row we have three cases: 

# (i) if 1 ≤ i ≤ r1, we set H(i, r_{11} + · · · + r{1,i-1} + 1) = 1 
# with the convention r{1,0} = 0. 

# (ii) if r1 + 1 ≤ i ≤ r1 + r2, we put 
# H(i, r_1^* + r_{21} + · · · + r_{2,i -r_1 - 1} + 1) = 1 with r{2,0} = 0. 

# (iii) if r1 + r2 + 1 ≤ i ≤ r, 
# we set H(i, r_1^* + r_2^* + r_31 + · · · + r_{3, i-r_1-r_2-1} + 1) = 1, 
# defining r_{3,0} = 0. The remaining entries are set equal to zero.

# 4. Modelos Factores ---------------------------------------------------------

## Factor 1

fac1 <- ts(Fac[,1], start = c(1981, 1), end = c(2013, 06), frequency = 12)
fac2 <- ts(Fac[,2], start = c(1981, 1), end = c(2013, 06), frequency = 12)
fac3 <- ts(Fac[,3], start = c(1981, 1), end = c(2013, 06), frequency = 12)
fac4 <- ts(Fac[,4], start = c(1981, 1), end = c(2013, 06), frequency = 12)

## 4.1 Gráfico de factores comunes --------------------------------------------
par(mfrow = c(2,2))
tsplot(fac1, col = "darkblue", lwd = 4,   
       xaxt = "n",   # remove x axis
       yaxt = "n",   # remove y axis
       xlab = "", mgpp = c(5,0.25,0), mar = c(1, 3, 1, 1),
       ylab = "")
axis(1, cex.axis = 2 , xlab = 'Fecha', cex.lab = 2)
axis(2, cex.axis = 2 , xlab = 'Fecha', cex.lab = 2)
mtext("Fecha", side = 1, line = 2, cex = 2)
mtext("Factor 1", side = 2, line = 2, cex = 2)

tsplot(fac2, col = "darkgreen", lwd = 4,   xaxt = "n",   # remove x axis
       yaxt = "n",   # remove y axis
       xlab = "", mgpp = c(5,0.25,0), mar = c(1, 3, 1, 1),
       ylab = "")
axis(1, cex.axis = 2 , cex.lab = 2)
axis(2, cex.axis = 2 , cex.lab = 2)
mtext("Fecha", side = 1, line = 2, cex = 2)
mtext("Factor 2", side = 2, line = 2, cex = 2)

tsplot(fac3, col = "darkred", lwd = 4,   xaxt = "n",   # remove x axis
       yaxt = "n",   # remove y axis
       xlab = "", mgpp = c(5,0.25,0), mar = c(1, 3, 1, 1),
       ylab = "")
axis(1, cex.axis = 2 , cex.lab = 2)
axis(2, cex.axis = 2 , cex.lab = 2)
mtext("Fecha", side = 1, line = 2, cex = 2)
mtext("Factor 3", side = 2, line = 2, cex = 2)

tsplot(fac4, col = "darkorange", lwd = 4,   xaxt = "n",   # remove x axis
       yaxt = "n",   # remove y axis
       xlab = "", mgpp = c(5,0.25,0), mar = c(1, 3, 1, 1),
       ylab = "")
axis(1, cex.axis = 2 , cex.lab = 2)
axis(2, cex.axis = 2 , cex.lab = 2)
mtext("Fecha", side = 1, line = 2, cex = 2)
mtext("Factor 4", side = 2, line = 2, cex = 2)

## 4.2 gráfico de subseries ---------------------------------------------------
p1 <- ggsubseriesplot(fac1, labels = 1:12, xlab = 'Mes') +
  labs(y = "Factor 1") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        #axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5),
        axis.text = element_text(size = 24))

p2 <- ggsubseriesplot(fac2, labels = 1:12, xlab = 'Mes') +
  labs(y = "Factor 2") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        #axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5),
        axis.text = element_text(size = 24))

p3 <- ggsubseriesplot(fac3, labels = 1:12, xlab = 'Mes') +
  labs(y = "Factor 3") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        #axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5),
        axis.text = element_text(size = 24))

p4 <- ggsubseriesplot(fac4, labels = 1:12, xlab = 'Mes') +
  labs(y = "Factor 4") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        #axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5),
        axis.text = element_text(size = 24))

grid.arrange(p1, p2, p3,p4, ncol = 2)

# primero chequea nsdiffs y luego prueba ADF
# H0: The time series is trend stationary.
# test_stat > crit_vals["5pct"] then "Reject H0: Series is NOT stationary"
#summary(ur.kpss(fac1)) # test-statistic = 0.2816 < Valor crítico No R. H0 (concluye estacionariedad)
#aTSA::adf.test(fac1, nlag = 12) # Type 1: no drift no trend
# p-value = 0.01 < 0.05 (R. H0, concluye que no hay raíz unitaria, por lo tanto es estacionaria)

## 4.3 Identificación modelos ARIMA/SARIMA ---------------------------
## 4.3.1 Factor 1 ----------------------------------------------------
#H0: No significant autocorrelation detected
#p.value < 0.05 then "Reject H0: The series has significant autocorrelation")
Box.test(fac1, lag = 24, type = "Ljung-Box") # p-value < 0.05 (R. H0: indep)

nsdiffs(fac1) # 1
# factor 1 estacional (no-estacionario) por tanto es de tipo r2
Dfac1 <- diff(fac1, lag = 12)
summary(ur.kpss(Dfac1)) # test-statistic = 0.0312 < Valor crítico No R. H0 (concluye estacionariedad)
aTSA::adf.test(Dfac1, nlag = 12)
# p-value = 0.01 < 0.05 (R. H0, concluye que no hay raíz unitaria, por lo tanto es estacionaria)

acf_f1 <- ggAcf(Dfac1, lag.max = 36) +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) # q = 1,2,3,4,5 y Q = 1

pacf_f1 <- ggPacf(Dfac1, lag.max = 36)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) 

grid.arrange(acf_f1, pacf_f1, ncol = 1)

# default: (5,2,5) x (2,1,2)
# ordenes sugeridos por ACF y PACF: (2,0,4) x (2,1,1)
# output: (1,0,2)(2,1,0)[12] 
auto.arima(fac1, d = 0, max.p = 2, max.q = 4, seasonal = TRUE)

# default (2,1,3) x (1,1,1)
# traza de ejecución:
# (2,0,4) x (2,1,1) -> Error
# (2,0,3) x (2,1,1) -> Error
# (2,0,3) x (2,1,1) -> Error
# (2,0,3) x (1,1,1) -> (1 0 1) x (1 0 1)
sarimaSpec(fac1, maxorder = c(2,0,2), maxsea = c(1,1,1), criterion = "aic",
           period = 12, output = TRUE,  method = "ML", include.mean = FALSE)
# para el Factor 1 se propone SARIMA(2,1,0)_12 y SARIMA(1,1,1)_12
## 4.3.2 Factor 2 ----------------------------------------------------
#H0: No significant autocorrelation detected
#p.value < 0.05 then "Reject H0: The series has significant autocorrelation")
Box.test(fac2, lag = 24, type = "Ljung-Box") # p-value < 0.05 (R. H0: indep)

nsdiffs(fac2) # 1
# factor 1 estacional (no-estacionario) por tanto es de tipo r2
Dfac2 <- diff(fac2, lag = 12)
summary(ur.kpss(Dfac2)) # test-statistic = 0.0337  < Valor crítico No R. H0 (concluye estacionariedad)
aTSA::adf.test(Dfac2, nlag = 12)
# p-value = 0.01 < 0.05 (R. H0, concluye que no hay raíz unitaria, por lo tanto es estacionaria)

acf_f2 <- ggAcf(Dfac2, lag.max = 36) +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) # q = 1,2,3,4,5 y Q = 1

pacf_f2 <- ggPacf(Dfac2, lag.max = 36)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) 

grid.arrange(acf_f2, pacf_f2, ncol = 1)

# default: (5,2,5) x (2,1,2)
# ordenes sugeridos por ACF y PACF: (1,0,1) x (3,1,1)
# output: (1,0,1)(2,1,0)[12]
auto.arima(fac2, d = 0, max.p = 1, max.q = 1, seasonal = TRUE)

# default (2,1,3) x (1,1,1)
# traza de ejecución:
# (1,0,1) x (3,1,1) -> Error
# (1,0,1) x (2,1,1) -> Error
# (1,0,1) x (2,1,1) -> (1 0 1) x (1 0 1)
sarimaSpec(fac2, maxorder = c(1,0,1), maxsea = c(1,1,1), criterion = "aic",
           period = 12, output = TRUE,  method = "ML", include.mean = FALSE)

# para el Factor 2 se propone SARIMA(2,1,0)_12 y SARIMA(1,1,1)_12
## 4.3.3 Factor 3 ----------------------------------------------------
#H0: No significant autocorrelation detected
#p.value < 0.05 then "Reject H0: The series has significant autocorrelation")
Box.test(fac3, lag = 24, type = "Ljung-Box") # p-value < 0.05 (R. H0: indep)

nsdiffs(fac3) # 1
# factor 1 estacional (no-estacionario) por tanto es de tipo r2
Dfac3 <- diff(fac3, lag = 12)
summary(ur.kpss(Dfac3)) # test-statistic = 0.0313   < Valor crítico No R. H0 (concluye estacionariedad)
aTSA::adf.test(Dfac3, nlag = 12)
# p-value = 0.01 < 0.05 (R. H0, concluye que no hay raíz unitaria, por lo tanto es estacionaria)

acf_f3 <- ggAcf(Dfac3, lag.max = 36) +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) # q = 1,2,3,4,5 y Q = 1

pacf_f3 <- ggPacf(Dfac3, lag.max = 36)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) 

grid.arrange(acf_f3, pacf_f3, ncol = 1)

# default: (5,2,5) x (2,1,2)
# ordenes sugeridos por ACF y PACF: (1,0,1) x (3,1,1)
# output: ARIMA(0,0,1)(1,1,0)[12]
auto.arima(fac3, d = 0, max.p = 1, max.q = 1, seasonal = TRUE)

# default (2,1,3) x (1,1,1)
# traza de ejecución:
# (1,0,1) x (3,1,1) -> (1 0 0) x (1 1 1)
sarimaSpec(fac3, maxorder = c(1,0,1), maxsea = c(3,1,1), criterion = "aic",
           period = 12, output = TRUE,  method = "ML", include.mean = FALSE)

# para el Factor 3 se propone SARIMA(1,1,0)_12 y SARIMA(1,1,1)_12

## 4.3.4 Factor 4 ----------------------------------------------------
#H0: No significant autocorrelation detected
#p.value < 0.05 then "Reject H0: The series has significant autocorrelation")
Box.test(fac4, lag = 24, type = "Ljung-Box") # p-value < 0.05 (R. H0: indep)

nsdiffs(fac4) # 0
# factor 1 estacional (no-estacionario) por tanto es de tipo r2
#Dfac3 <- diff(fac3, lag = 12)
summary(ur.kpss(fac4)) # test-statistic = 1.6639  > Valor crítico R. H0 (concluye no estacionariedad)
aTSA::adf.test(fac4, nlag = 12)
# p-value = 0.04 < 0.05 (R. H0, concluye que no hay raíz unitaria, por lo tanto es estacionaria)
dfac4 = diff(fac4, lag = 1)
summary(ur.kpss(dfac4)) # test-statistic = 0.0086 < Valor crítico No R. H0 (concluye estacionariedad)
aTSA::adf.test(dfac4, nlag = 12)
# p-value = 0.01 < 0.05 (R. H0, concluye que no hay raíz unitaria, por lo tanto es estacionaria)
nsdiffs(dfac4) # 0

acf_f4 <- ggAcf(dfac4, lag.max = 36) +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) # q = 1,2,3,4,5 y Q = 1

pacf_f4 <- ggPacf(dfac4, lag.max = 36)+
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        title = element_text(size = 24)) 

grid.arrange(acf_f4, pacf_f4, ncol = 1)

# default: (5,2,5)
# ordenes sugeridos por ACF y PACF: (11,1,3)
# output: ARIMA(3,1,1)
auto.arima(fac4, d = 1, max.p = 11, max.q = 3, seasonal = FALSE)

# default (2,1,3) x (1,1,1)
# traza de ejecución:
# (1,0,1) x (3,1,1) -> (1 0 0) x (1 1 1)
#sarimaSpec(fac4, maxorder = c(1,0,1), maxsea = c(3,1,1), criterion = "aic",
#           period = 12, output = TRUE,  method = "ML", include.mean = FALSE)

# default (5,1,4)
# traza de ejecución:
# (11,1,3) -> Error
# (7,1,3) -> Error
# (6,1,3) -> (4,0,3)
arimaSpec(fac4, maxorder = c(6,1,3), criterion="aic",
          output = FALSE, method = "CSS-ML", pv = 0.01)

# para el Factor 4 se propone ARIMA(3,1,1) y ARIMA(4,1,3)

# 4.4 Estimación de parámetros ----------------------------------------------
## 4.4.1 Factor 1 --------------------------------- 
# modelo1 factor1 SARIMA(2,1,0)_12
sarima(fac1, p = 0, d = 0, q = 0, P = 2, D = 1, Q = 0, S = 12,
       no.constant = T)
SARIMA_f1_m1 <- forecast::Arima(fac1, order = c(0, 0, 0), 
                                seasonal = list(order = c(2, 1, 0), period = 12))
summary(SARIMA_f1_m1) # coeficientes significativos
checkresiduals(SARIMA_f1_m1)

Box.test(SARIMA_f1_m1$residuals, lag = (length(SARIMA_f1_m1$residuals)/4), 
         type = "Ljung-Box", fitdf = 2)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f1_m1$residuals) # pvalue = 0.1439 > 0.05 (No R. H0, concluye normalidad)

cusum_g(SARIMA_f1_m1$residuals)

# modelo2 factor1 SARIMA(0,1,1)_12
sarima(fac1, p = 0, d = 0, q = 0, P = 0, D = 1, Q = 1, S = 12,
       no.constant = T)
SARIMA_f1_m2 <- forecast::Arima(fac1, order = c(0, 0, 0), 
                                seasonal = list(order = c(0, 1, 1), period = 12))
summary(SARIMA_f1_m2) # sar1 No significativo
checkresiduals(SARIMA_f1_m2)

Box.test(SARIMA_f1_m2$residuals, lag = (length(SARIMA_f1_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 1)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f1_m2$residuals) # pvalue = 2.034e-05 < 0.05 (R. H0, No concluye normalidad)
cusum_g(SARIMA_f1_m2$residuals)

## 4.4.2 Factor 2 --------------------------------- 
# modelo 1 factor2 SARIMA(2,1,0)_12 -> No cumple supuesto L-B, cumple supuesto JB al 1% p-value = 0.034
sarima(fac2, p = 0, d = 0, q = 0, P = 2, D = 1, Q = 0, S = 12,
       no.constant = T)
SARIMA_f2_m1 <- forecast::Arima(fac2, order = c(0, 0, 0), 
                                seasonal = list(order = c(2, 1, 0), period = 12))
summary(SARIMA_f2_m1) # coeficientes significativos
checkresiduals(SARIMA_f2_m1)

Box.test(SARIMA_f2_m1$residuals, lag = (length(SARIMA_f2_m1$residuals)/4), 
         type = "Ljung-Box", fitdf = 2)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f2_m1$residuals) # pvalue = 0.03468 < 0.05 (R. H0, No concluye normalidad)

cusum_g(SARIMA_f2_m1$residuals)

# modelo 2 factor2 SARIMA(0,1,1) -> No cumple supuesto JB
sarima(fac2, p = 0, d = 0, q = 0, P = 0, D = 1, Q = 1, S = 12,
       no.constant = T)
SARIMA_f2_m2 <- forecast::Arima(fac2, order = c(0, 0, 0), 
                                seasonal = list(order = c(0, 1, 1), period = 12))
summary(SARIMA_f2_m2) # sar1 No significativo
checkresiduals(SARIMA_f2_m2)

Box.test(SARIMA_f2_m2$residuals, lag = (length(SARIMA_f2_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 1)
# p-value > 0.05 (No R. H0, concluye no significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f2_m2$residuals) # pvalue = 2.034e-05 < 0.05 (R. H0, No concluye normalidad)
cusum_g(SARIMA_f2_m2$residuals)

## 4.4.3 Factor 3 --------------------------------- 
# modelo 1 factor 3 SARIMA(1,1,0) -> No cumple supuesto L-B
sarima(fac3, p = 0, d = 0, q = 0, P = 1, D = 1, Q = 0, S = 12,
       no.constant = T)
SARIMA_f3_m1 <- forecast::Arima(fac3, order = c(0, 0, 0), 
                                seasonal = list(order = c(1, 1, 0), period = 12))
summary(SARIMA_f3_m1) # coeficientes significativos
checkresiduals(SARIMA_f3_m1)

Box.test(SARIMA_f3_m1$residuals, lag = (length(SARIMA_f3_m1$residuals)/4), 
         type = "Ljung-Box", fitdf = 1)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f3_m1$residuals) # pvalue = 0.8861 > 0.05 (No R. H0, concluye normalidad)

cusum_g(SARIMA_f3_m1$residuals)

# modelo 2 factor 3 SARIMA(0,1,1) -> cumple supuesto JB al 1% p-value = 0.046
sarima(fac3, p = 0, d = 0, q = 0, P = 0, D = 1, Q = 1, S = 12,
       no.constant = T)
SARIMA_f3_m2 <- forecast::Arima(fac3, order = c(0, 0, 0), 
                                seasonal = list(order = c(0, 1, 1), period = 12))
summary(SARIMA_f3_m2) # sar1 No significativo
checkresiduals(SARIMA_f3_m2)

Box.test(SARIMA_f3_m2$residuals, lag = (length(SARIMA_f3_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 1)
# p-value > 0.05 (No R. H0, concluye no significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f3_m2$residuals) # pvalue = 0.04614 < 0.05 (R. H0, No concluye normalidad)
cusum_g(SARIMA_f3_m2$residuals)

## 4.4.4 Factor 4 --------------------------------- 
# modelo 1 factor 4 ARIMA(3,1,1) -> coef AR3 no significativo y no cumple supuestos L-B y JB
# modelo 1 factor 4 ARIMA(2,1,1) -> No cumple supuestos L-B y JB
# modelo 1 factor 4 ARIMA(3,1,0) -> No cumple supuesto L-B

sarima(fac4, p = 3, d = 1, q = 0, P = 0, D = 0, Q = 0, S = 0,
       no.constant = T)
ARIMA_f4_m1 <- forecast::Arima(fac4, order = c(3, 1, 0), 
                                seasonal = list(order = c(0, 0, 0), period = 0))
summary(ARIMA_f4_m1) # coeficientes significativos
checkresiduals(ARIMA_f4_m1)

Box.test(ARIMA_f4_m1$residuals, lag = (length(ARIMA_f4_m1$residuals)/4), 
         type = "Ljung-Box", fitdf = 3)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(ARIMA_f4_m1$residuals) # pvalue = 0.4983 > 0.05 (No R. H0, concluye normalidad)

cusum_g(ARIMA_f4_m1$residuals)

# modelo 2 factor 4 ARIMA(4,1,3) -> coef AR4 no significativo, No cumple supuestos L-B y JB
# modelo 2 factor 4 ARIMA(3,1,3) -> No cumple supuestos L-B y JB
# modelo 2 factor 4 ARIMA(3,1,2) -> coef AR2 y MA2 no significativos
# modelo 2 factor 4 ARIMA(2,1,2) -> coef MA2 no significativo y No cumple supuestos L-B y JB
# modelo 2 factor 4 ARIMA(2,1,1) -> No cumple supuestos L-B y JB
# modelo 2 factor 4 ARIMA(2,1,0) -> No cumple supuesto L-B

sarima(fac4, p = 2, d = 1, q = 0, P = 0, D = 0, Q = 0, S = 0,
       no.constant = T)

ARIMA_f4_m2 <- forecast::Arima(fac4, order = c(2, 1, 0), 
                                seasonal = list(order = c(0, 0, 0), period = 0))
summary(ARIMA_f4_m2) # sar1 No significativo
checkresiduals(ARIMA_f4_m2)

Box.test(ARIMA_f4_m2$residuals, lag = (length(ARIMA_f4_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 2)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(ARIMA_f4_m2$residuals) # pvalue = 0.7818 < 0.05 (R. H0, No concluye normalidad)
cusum_g(ARIMA_f4_m2$residuals)

# 4.5 Tabla resumen modelos y chequeo supuestos ------------------------------
pcoef_ = function(x){paste(round(as.vector(x$coef),2), collapse = ',')}

p_vale_box_test = function(x){
  Box.test(x$residuals, lag = (length(x$residuals)/4), type = "Ljung-Box", 
           fitdf = length(x$coef))$p.value}

p_vale_Jarque_Bera = function(x){
  tseries::jarque.bera.test(x$residuals)$p.value}

tbl_models_factors = data.frame(
  Factor = c('f_1', 'f_1', 'f_2', 'f_2', 'f_3', 'f_3', 'f_4'),
  Modelo = c('SARIMA(2,1,0)_12', 'SARIMA(0,1,1)_12',
             'SARIMA(2,1,0)_12', 'SARIMA(0,1,1)_12',
             'SARIMA(1,1,0)_12', 'SARIMA(0,1,1)_12',
             'ARIMA(3,1,0)'),
  Coeficientes = c(pcoef_(SARIMA_f1_m1), pcoef_(SARIMA_f1_m2),
                   pcoef_(SARIMA_f2_m1), pcoef_(SARIMA_f2_m2),
                   pcoef_(SARIMA_f3_m1), pcoef_(SARIMA_f3_m2),
                   pcoef_(ARIMA_f4_m1)),
  p_valor_Ljung_Box = c(p_vale_box_test(SARIMA_f1_m1), p_vale_box_test(SARIMA_f1_m2),
                        p_vale_box_test(SARIMA_f2_m1), p_vale_box_test(SARIMA_f2_m2),
                        p_vale_box_test(SARIMA_f3_m1), p_vale_box_test(SARIMA_f3_m2),
                        p_vale_box_test(ARIMA_f4_m1)),
  p_valor_Jarque_Bera = c(p_vale_Jarque_Bera(SARIMA_f1_m1), p_vale_Jarque_Bera(SARIMA_f1_m2),
                          p_vale_Jarque_Bera(SARIMA_f2_m1), p_vale_Jarque_Bera(SARIMA_f2_m2),
                          p_vale_Jarque_Bera(SARIMA_f3_m1), p_vale_Jarque_Bera(SARIMA_f3_m2),
                          p_vale_Jarque_Bera(ARIMA_f4_m1))
)

xtable(
  tbl_models_factors,
  caption = "",
  label = "tab:km4",
  align = c("l", "c", "c", "c","c", "c") # column alignment: left, right, center
)

# 5. rolling forecast ----------------------------------------------------

## 5.1 scale(X2): datos transformados ----------------------------------------
h = 1
ntrain = length(fac1)-1 #390 is the same for length(fac2)
ntest <- nrow(sX2) - ntrain # 90
time(fac1)[ntrain] # 2013.417 is the same for time(fac2)[ntrain]

eh_ <- matrix(0, nrow = ntest, ncol = ncol(sX2))
seh <- list(eh1 = eh_, eh12 = eh_, eh21 = eh_, eh2 = eh_)

Chi_hat_ <- matrix(0, nrow = ntest, ncol = ncol(sX2))
sChi_hat <- list(Chi_hat1 = Chi_hat_, Chi_hat12 = Chi_hat_, 
                Chi_hat21 = Chi_hat_, Chi_hat2 = Chi_hat_)

for(i in 1:(ntest)){
  step = ntrain + i
  Fac_ = sX2[1:step,] %*% Load
  fac1_ <- ts(Fac_[,1], start = c(1981, 1), end = time(fac1)[ntrain]+(i-1)/12, frequency = 12)
  fac2_ <- ts(Fac_[,2], start = c(1981, 1), end = time(fac2)[ntrain]+(i-1)/12, frequency = 12)
  fac3_ <- ts(Fac_[,3], start = c(1981, 1), end = time(fac3)[ntrain]+(i-1)/12, frequency = 12)
  fac4_ <- ts(Fac_[,4], start = c(1981, 1), end = time(fac4)[ntrain]+(i-1)/12, frequency = 12)
  
  fac1_=window(fac1_,end=time(fac1)[ntrain]+(i-1)/12)
  fac2_=window(fac2_,end=time(fac2)[ntrain]+(i-1)/12)
  fac3_=window(fac3_,end=time(fac3)[ntrain]+(i-1)/12)
  fac4_=window(fac4_,end=time(fac4)[ntrain]+(i-1)/12)
  # refit SARIMA models with trainning data
  #re_Sarima_f1_m1 = Arima(y = fac1_, c(0, 1, 0), seasonal = list(order = c(0, 1, 1), period = 12)) # model paper SARIMA(0,1,1)
  #re_Sarima_f2_m1 = Arima(y = fac2_, c(0, 0, 0), seasonal = list(order = c(0, 1, 1), period = 12)) # model paper SARIMA(1,1,0)
  
  # add new data to model with fixed params
  
  Sarima_f1_m1 = Arima(y = fac1_, c(0, 0, 0), fixed=SARIMA_f1_m1$coef,
                       seasonal = list(order = c(2, 1, 0), period = 12)) # f1-m1 SARIMA(2,1,0)
  Sarima_f2_m1 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m1$coef,
                       seasonal = list(order = c(2, 1, 0), period = 12)) # f2-m1 SARIMA(2,1,0)
  Sarima_f2_m2 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m2$coef,
                       seasonal = list(order = c(0, 1, 1), period = 12)) # f2-m2 SARIMA(0,1,1)
  Sarima_f3_m1 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m1$coef,
                       seasonal = list(order = c(1, 1, 0), period = 12)) # f3-m1 SARIMA(1,1,0)
  Sarima_f3_m2 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m2$coef,
                       seasonal = list(order = c(0, 1, 1), period = 12)) # f3-m1 SARIMA(0,1,1)
  arima_f4_m1 = Arima(y = fac4_, c(3, 1, 0), fixed=ARIMA_f4_m1$coef,
                       seasonal = list(order = c(0, 0, 0), period = 0)) # f4-m1 ARIMA(3,1,0)
  
  # forecast with factors values 1-step ahead
  f1_hat_m1 = forecast::forecast(Sarima_f1_m1, h = h)$mean
  f2_hat_m1 = forecast::forecast(Sarima_f2_m1, h = h)$mean
  f2_hat_m2 = forecast::forecast(Sarima_f2_m2, h = h)$mean
  f3_hat_m1 = forecast::forecast(Sarima_f3_m1, h = h)$mean
  f3_hat_m2 = forecast::forecast(Sarima_f3_m2, h = h)$mean
  f4_hat_m1 = forecast::forecast(arima_f4_m1, h = h)$mean
  # get forecast for the real values of each time series in scaled data object sX2
  sChi_hat$Chi_hat1[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m1, f4_hat_m1) %*% t(Load)
  sChi_hat$Chi_hat12[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m2, f4_hat_m1) %*% t(Load)
  sChi_hat$Chi_hat21[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m1, f4_hat_m1) %*% t(Load)
  sChi_hat$Chi_hat2[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m2, f4_hat_m1) %*% t(Load)
  # get real data
  obs <- sX2[ntrain + i, ]
  # save error
  seh$eh1[i, ] <- (obs - sChi_hat$Chi_hat1[i,])^2
  seh$eh12[i, ] <- (obs - sChi_hat$Chi_hat12[i,])^2
  seh$eh21[i, ] <- (obs - sChi_hat$Chi_hat21[i,])^2
  seh$eh2[i, ] <- (obs - sChi_hat$Chi_hat2[i,])^2
}

lapply(seh, mean)

## 5.1.1 rolling h-pasos adelante --------------------------------------------
rolling_sX2 <- list()
for(h in 1:12){
  ntrain = length(fac1) - 1  #390 is the same for length(fac2)
  ntest <- nrow(sX2) - ntrain - h  # 90
  time(fac1)[ntrain] # 2013.417 is the same for time(fac2)[ntrain]
  
  eh_ <- matrix(0, nrow = ntest, ncol = ncol(sX2))
  seh2 <- list(eh1 = eh_, eh12 = eh_, eh21 = eh_, eh2 = eh_)
  
  Chi_hat_ <- matrix(0, nrow = ntest, ncol = ncol(sX2))
  sChi_hat2 <- list(Chi_hat1 = Chi_hat_, Chi_hat12 = Chi_hat_, 
                   Chi_hat21 = Chi_hat_, Chi_hat2 = Chi_hat_)
  
  for(i in 1:ntest){
    step = ntrain + i
    Fac_ = sX2[1:step,] %*% Load
    fac1_ <- ts(Fac_[,1], start = c(1981, 1), end = time(fac1)[ntrain]+(i)/12, frequency = 12)
    fac2_ <- ts(Fac_[,2], start = c(1981, 1), end = time(fac2)[ntrain]+(i)/12, frequency = 12)
    fac3_ <- ts(Fac_[,3], start = c(1981, 1), end = time(fac3)[ntrain]+(i)/12, frequency = 12)
    fac4_ <- ts(Fac_[,4], start = c(1981, 1), end = time(fac4)[ntrain]+(i)/12, frequency = 12)
    
    fac1_=window(fac1_,end=time(fac1)[ntrain]+(i)/12)
    fac2_=window(fac2_,end=time(fac2)[ntrain]+(i)/12)
    fac3_=window(fac3_,end=time(fac3)[ntrain]+(i)/12)
    fac4_=window(fac4_,end=time(fac4)[ntrain]+(i)/12)
    
    # add new data to model with fixed params
    
    Sarima_f1_m1 = Arima(y = fac1_, c(0, 0, 0), fixed=SARIMA_f1_m1$coef,
                         seasonal = list(order = c(2, 1, 0), period = 12)) # f1-m1 SARIMA(2,1,0)
    Sarima_f2_m1 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m1$coef,
                         seasonal = list(order = c(2, 1, 0), period = 12)) # f2-m1 SARIMA(2,1,0)
    Sarima_f2_m2 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m2$coef,
                         seasonal = list(order = c(0, 1, 1), period = 12)) # f2-m2 SARIMA(0,1,1)
    Sarima_f3_m1 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m1$coef,
                         seasonal = list(order = c(1, 1, 0), period = 12)) # f3-m1 SARIMA(1,1,0)
    Sarima_f3_m2 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m2$coef,
                         seasonal = list(order = c(0, 1, 1), period = 12)) # f3-m1 SARIMA(0,1,1)
    arima_f4_m1 = Arima(y = fac4_, c(3, 1, 0), fixed=ARIMA_f4_m1$coef,
                        seasonal = list(order = c(0, 0, 0), period = 0)) # f4-m1 ARIMA(3,1,0)
    
    # forecast with factors values 1-step ahead
    f1_hat_m1 = forecast::forecast(Sarima_f1_m1, h = h)$mean[h]
    f2_hat_m1 = forecast::forecast(Sarima_f2_m1, h = h)$mean[h]
    f2_hat_m2 = forecast::forecast(Sarima_f2_m2, h = h)$mean[h]
    f3_hat_m1 = forecast::forecast(Sarima_f3_m1, h = h)$mean[h]
    f3_hat_m2 = forecast::forecast(Sarima_f3_m2, h = h)$mean[h]
    f4_hat_m1 = forecast::forecast(arima_f4_m1, h = h)$mean[h]
    # get forecast for the real values of each time series in scaled data object sX2
    sChi_hat2$Chi_hat1[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m1, f4_hat_m1) %*% t(Load)
    sChi_hat2$Chi_hat12[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m2, f4_hat_m1) %*% t(Load)
    sChi_hat2$Chi_hat21[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m1, f4_hat_m1) %*% t(Load)
    sChi_hat2$Chi_hat2[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m2, f4_hat_m1) %*% t(Load)
    # get real data
    obs <- sX2[ntrain + i + h, ]
    # save error
    seh2$eh1[i, ] <- (obs - sChi_hat2$Chi_hat1[i,])^2
    seh2$eh12[i, ] <- (obs - sChi_hat2$Chi_hat12[i,])^2
    seh2$eh21[i, ] <- (obs - sChi_hat2$Chi_hat21[i,])^2
    seh2$eh2[i, ] <- (obs - sChi_hat2$Chi_hat2[i,])^2
  }
  
  rolling_sX2[[h]] <- lapply(seh2, mean) %>% bind_cols()
}

df_MSE_rolling_sX2 = bind_rows(rolling_sX2)

xtable(
  df_MSE_rolling_sX2[c(1,3,6,12),],
  caption = "", digits = 3,
  label = "tab:km4",
  align = c("l", "c", "c", "c", "c") # column alignment: left, right, center
)

## 5.2 X2: datos originales  --------------------------------------------

## 5.2.1 rolling 1-paso adelante -----------------------------------------
# misma conclusión que con scale(X2), la combinación m1-m2-m2-m1 la mejor
h = 1
ntrain = length(fac1)-1 #390 is the same for length(fac2)
ntest <- nrow(X2) - ntrain # 90
time(fac1)[ntrain] # 2013.417 is the same for time(fac2)[ntrain]

eh_ <- matrix(0, nrow = ntest, ncol = ncol(X2))
eh <- list(eh1 = eh_, eh12 = eh_, eh21 = eh_, eh2 = eh_)

Chi_hat_ <- matrix(0, nrow = ntest, ncol = ncol(X2))
Chi_hat <- list(Chi_hat1 = Chi_hat_, Chi_hat12 = Chi_hat_, 
                Chi_hat21 = Chi_hat_, Chi_hat2 = Chi_hat_)

for(i in 1:(ntest)){
  step = ntrain + i
  Fac_ = X2[1:step,] %*% Load
  fac1_ <- ts(Fac_[,1], start = c(1981, 1), end = time(fac1)[ntrain]+(i-1)/12, frequency = 12)
  fac2_ <- ts(Fac_[,2], start = c(1981, 1), end = time(fac2)[ntrain]+(i-1)/12, frequency = 12)
  fac3_ <- ts(Fac_[,3], start = c(1981, 1), end = time(fac3)[ntrain]+(i-1)/12, frequency = 12)
  fac4_ <- ts(Fac_[,4], start = c(1981, 1), end = time(fac4)[ntrain]+(i-1)/12, frequency = 12)
  
  fac1_=window(fac1_,end=time(fac1)[ntrain]+(i-1)/12)
  fac2_=window(fac2_,end=time(fac2)[ntrain]+(i-1)/12)
  fac3_=window(fac3_,end=time(fac3)[ntrain]+(i-1)/12)
  fac4_=window(fac4_,end=time(fac4)[ntrain]+(i-1)/12)
  
  # add new data to model with fixed params
  
  Sarima_f1_m1 = Arima(y = fac1_, c(0, 0, 0), fixed=SARIMA_f1_m1$coef,
                       seasonal = list(order = c(2, 1, 0), period = 12)) # f1-m1 SARIMA(2,1,0)
  Sarima_f2_m1 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m1$coef,
                       seasonal = list(order = c(2, 1, 0), period = 12)) # f2-m1 SARIMA(2,1,0)
  Sarima_f2_m2 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m2$coef,
                       seasonal = list(order = c(0, 1, 1), period = 12)) # f2-m2 SARIMA(0,1,1)
  Sarima_f3_m1 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m1$coef,
                       seasonal = list(order = c(1, 1, 0), period = 12)) # f3-m1 SARIMA(1,1,0)
  Sarima_f3_m2 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m2$coef,
                       seasonal = list(order = c(0, 1, 1), period = 12)) # f3-m1 SARIMA(0,1,1)
  arima_f4_m1 = Arima(y = fac4_, c(3, 1, 0), fixed=ARIMA_f4_m1$coef,
                      seasonal = list(order = c(0, 0, 0), period = 0)) # f4-m1 ARIMA(3,1,0)
  
  # forecast with factors values 1-step ahead
  f1_hat_m1 = forecast::forecast(Sarima_f1_m1, h = h)$mean
  f2_hat_m1 = forecast::forecast(Sarima_f2_m1, h = h)$mean
  f2_hat_m2 = forecast::forecast(Sarima_f2_m2, h = h)$mean
  f3_hat_m1 = forecast::forecast(Sarima_f3_m1, h = h)$mean
  f3_hat_m2 = forecast::forecast(Sarima_f3_m2, h = h)$mean
  f4_hat_m1 = forecast::forecast(arima_f4_m1, h = h)$mean
  # get forecast for the real values of each time series in scaled data object sX2
  Chi_hat$Chi_hat1[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m1, f4_hat_m1) %*% t(Load)
  Chi_hat$Chi_hat12[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m2, f4_hat_m1) %*% t(Load)
  Chi_hat$Chi_hat21[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m1, f4_hat_m1) %*% t(Load)
  Chi_hat$Chi_hat2[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m2, f4_hat_m1) %*% t(Load)
  # get real data
  obs <- X2[ntrain + i, ]
  # save error
  eh$eh1[i, ] <- (obs - Chi_hat$Chi_hat1[i,])^2
  eh$eh12[i, ] <- (obs - Chi_hat$Chi_hat12[i,])^2
  eh$eh21[i, ] <- (obs - Chi_hat$Chi_hat21[i,])^2
  eh$eh2[i, ] <- (obs - Chi_hat$Chi_hat2[i,])^2
}

lapply(eh, mean)

## 5.2.2 rolling h-paso adelante -----------------------------------------
# misma conclusión que con scale(X2), la combinación m1-m2-m2-m1 la mejor
rolling_X2 <- list()
Chi_hat2_ <- list()
for(h in 1:12){
  ntrain = length(fac1)-1 #390 is the same for length(fac2)
  #ntrain = length(fac1)#-1 #390 is the same for length(fac2)
  ntest <- nrow(X2) - ntrain - h # 90
  #ntest <- nrow(X2) - ntrain - h # 90
  time(fac1)[ntrain] # 2013.417 is the same for time(fac2)[ntrain]
  
  eh_ <- matrix(0, nrow = ntest, ncol = ncol(X2))
  eh2 <- list(eh1 = eh_, eh12 = eh_, eh21 = eh_, eh2 = eh_)
  
  Chi_hat_ <- matrix(0, nrow = ntest, ncol = ncol(X2))
  Chi_hat2_[[h]] <- list(Chi_hat1 = Chi_hat_, Chi_hat12 = Chi_hat_, 
                         Chi_hat21 = Chi_hat_, Chi_hat2 = Chi_hat_)
  
  for(i in 1:ntest){
    step = ntrain + i
    Fac_ = X2[1:step,] %*% Load
    fac1_ <- ts(Fac_[,1], start = c(1981, 1), end = time(fac1)[ntrain]+(i)/12, frequency = 12)
    fac2_ <- ts(Fac_[,2], start = c(1981, 1), end = time(fac2)[ntrain]+(i)/12, frequency = 12)
    fac3_ <- ts(Fac_[,3], start = c(1981, 1), end = time(fac3)[ntrain]+(i)/12, frequency = 12)
    fac4_ <- ts(Fac_[,4], start = c(1981, 1), end = time(fac4)[ntrain]+(i)/12, frequency = 12)
    
    fac1_=window(fac1_,end=time(fac1)[ntrain]+(i)/12)
    fac2_=window(fac2_,end=time(fac2)[ntrain]+(i)/12)
    fac3_=window(fac3_,end=time(fac3)[ntrain]+(i)/12)
    fac4_=window(fac4_,end=time(fac4)[ntrain]+(i)/12)
    
    # add new data to model with fixed params
    
    Sarima_f1_m1 = Arima(y = fac1_, c(0, 0, 0), fixed=SARIMA_f1_m1$coef,
                         seasonal = list(order = c(2, 1, 0), period = 12)) # f1-m1 SARIMA(2,1,0)
    Sarima_f2_m1 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m1$coef,
                         seasonal = list(order = c(2, 1, 0), period = 12)) # f2-m1 SARIMA(2,1,0)
    Sarima_f2_m2 = Arima(y = fac2_, c(0, 0, 0), fixed=SARIMA_f2_m2$coef,
                         seasonal = list(order = c(0, 1, 1), period = 12)) # f2-m2 SARIMA(0,1,1)
    Sarima_f3_m1 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m1$coef,
                         seasonal = list(order = c(1, 1, 0), period = 12)) # f3-m1 SARIMA(1,1,0)
    Sarima_f3_m2 = Arima(y = fac3_, c(0, 0, 0), fixed=SARIMA_f3_m2$coef,
                         seasonal = list(order = c(0, 1, 1), period = 12)) # f3-m1 SARIMA(0,1,1)
    arima_f4_m1 = Arima(y = fac4_, c(3, 1, 0), fixed=ARIMA_f4_m1$coef,
                        seasonal = list(order = c(0, 0, 0), period = 0)) # f4-m1 ARIMA(3,1,0)
    
    # forecast with factors values 1-step ahead
    f1_hat_m1 = forecast::forecast(Sarima_f1_m1, h = h)$mean[h]
    f2_hat_m1 = forecast::forecast(Sarima_f2_m1, h = h)$mean[h]
    f2_hat_m2 = forecast::forecast(Sarima_f2_m2, h = h)$mean[h]
    f3_hat_m1 = forecast::forecast(Sarima_f3_m1, h = h)$mean[h]
    f3_hat_m2 = forecast::forecast(Sarima_f3_m2, h = h)$mean[h]
    f4_hat_m1 = forecast::forecast(arima_f4_m1, h = h)$mean[h]
    # get forecast for the real values of each time series in scaled data object sX2
    Chi_hat2_[[h]]$Chi_hat1[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m1, f4_hat_m1) %*% t(Load)
    Chi_hat2_[[h]]$Chi_hat12[i,] <- cbind(f1_hat_m1, f2_hat_m1, f3_hat_m2, f4_hat_m1) %*% t(Load)
    Chi_hat2_[[h]]$Chi_hat21[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m1, f4_hat_m1) %*% t(Load)
    Chi_hat2_[[h]]$Chi_hat2[i,] <- cbind(f1_hat_m1, f2_hat_m2, f3_hat_m2, f4_hat_m1) %*% t(Load)
    # get real data
    obs <- X2[ntrain + i + h, ]
    # save error
    eh2$eh1[i, ] <- (obs - Chi_hat2_[[h]]$Chi_hat1[i,])^2
    eh2$eh12[i, ] <- (obs - Chi_hat2_[[h]]$Chi_hat12[i,])^2
    eh2$eh21[i, ] <- (obs - Chi_hat2_[[h]]$Chi_hat21[i,])^2
    eh2$eh2[i, ] <- (obs - Chi_hat2_[[h]]$Chi_hat2[i,])^2
  }
  
  rolling_X2[[h]] <- lapply(eh2, mean) %>% bind_cols()
}

df_MSE_rolling = bind_rows(rolling_X2)

df_MSE_rolling[c(1,3,6,12),]

# 6. Visualización Factores y Cargas ------------------------------------------

## 6.1 Identifica la curva más profunda y las menos profundas del cluster 4 -----
df_clust4_wide = filter(montlyAvg_pcp_mm2_, clust_wtv2 == 4)%>% 
  pivot_wider(id_cols = 'mes', names_from = 'ID', values_from = 'pcp') %>%
  tibble::column_to_rownames("mes") %>% as.matrix()

dim(df_clust4_wide) # 12 fil x 352 col

monthbasis11 = monthbasis11 <- create.fourier.basis(rangeval=c(0, 12), nbasis=11, period = 12,
                                                    axes=list('axesIntervals'))
pcp.fd = smooth.basisPar(argvals = 1:12, df_clust4_wide, monthbasis11)$fd

eval_pcp.fd <- eval.fd(1:12, pcp.fd)

pal5 = scales::hue_pal()(5)

bx_clust4 = boxplot(pcp.fd)
quantile(bx_clust4$depth) 
# Q2: 0.4931152  medcurve
depths_ = as.vector(bx_clust4$depth)
which(abs(depths_ - 0.4931152)<0.001)
which(bx_clust4$depth == depths_[309]) # ID: 26100070 
min(bx_clust4$depth)

# Q_: 0.1079837 (menos profunda)
which(abs(depths_ - 0.1079837  )<0.05)
which(bx_clust4$depth == depths_[87]) # ID: 21206060
which(bx_clust4$depth == depths_[159]) # ID: 21080070   

boxplot(pcp.fd)
lines(eval_pcp.fd[,309], col = 'green', lwd = 2)
lines(eval_pcp.fd[,87], col = 'darkorange', lwd = 2)
lines(eval_pcp.fd[,159], col = 'darkorange', lwd = 2)

## 6.1 Plot Rolling forecast "LUCERNA HACIENDA [26100070]" (curva mediana) cluster 4 -----
dim(Chi_hat2_[[1]]$Chi_hat2) # 90 meses x 352 sttns
names(Data)[str_detect(string = names(Data), pattern = '26100070') ]
n_id = which(names(Data) == "LUCERNA HACIENDA [26100070]")

df_clust4_med <- tibble(
  Time = Data2$Date[391:480],
  pcp = X2[391:480, n_id],
  DFM_h1 = Chi_hat2_[[1]]$Chi_hat2[,n_id],
  DFM_h12 = c(rep(NA,11),Chi_hat2_[[12]]$Chi_hat2[,n_id])) %>%
  pivot_longer(cols = -Time, names_to = "Series", values_to = "Value")

df_MSE_clust4_med = df_clust4_med %>%  
  pivot_wider(id_cols = Time, values_from = 'Value', names_from = 'Series') %>% 
  mutate(MSE_1 = (pcp-DFM_h1)^2, MSE_12 = (pcp-DFM_h12)^2) 

RMSE_clust4_fbx = data.frame(h = c(1,12), RMSE_med = c(0,0), RMSE_low = c(0,0))

RMSE_clust4_fbx$RMSE_med[1] <- sqrt(mean(df_MSE_clust4_med$MSE_1))
RMSE_clust4_fbx$RMSE_med[2] <- sqrt(mean(df_MSE_clust4_med$MSE_12, na.rm = TRUE))

labels_ = c('Observado','Pronóstico 1-paso adelante', 'Pronóstico 12-pasos adelante')
rol_clust4_med <- ggplot(df_clust4_med, aes(x = Time, y = Value, color = Series, linetype = Series)) +
  geom_line(linewidth = 0.9) +
  scale_linetype_manual(values = c(pcp = "solid", DFM_h1 = "dashed", DFM_h12 = "dashed"),
                        name ='',
                        breaks = c("pcp","DFM_h1","DFM_h12"),
                        labels = labels_)+
  scale_color_manual(values = c(pcp = "black", DFM_h1= "red", DFM_h12= "purple"),
                     name ='',
                     breaks = c("pcp","DFM_h1","DFM_h12"),
                     labels = labels_) +
  labs(x = "Fecha", y = "Precipitación (mm)", title = str_to_title(names(Data2)[n_id])) +
  coord_cartesian(ylim =c(0,400))+
  theme_bw() +
  theme(
    title = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 24),
    legend.position = c(0.2, 0.9),     # top-left inside
    legend.text = element_text(size = 24),
    legend.background = element_blank())

## 6.1.2 Plot Rolling forecast "LUCERNA HACIENDA [26100070]" (curva mediana) cluster 4 -----
names(Data2)[str_detect(string = names(Data2), pattern = '21206060') ]
n_id2 = which(names(Data2) == "CASABLANCA [21206060]")

#names(Data2)[str_detect(string = names(Data2), pattern = '21080070') ]
#which(names(Data2) == "SANTA ROSA HACIENDA [21080070]")

df_clust4_low <- tibble(
  Time = Data2$Date[391:480],
  pcp = X2[391:480, n_id2],
  DFM_h1 = Chi_hat2_[[1]]$Chi_hat2[,n_id2],
  DFM_h12 = c(rep(NA,11),Chi_hat2_[[12]]$Chi_hat2[,n_id2])) %>%
  pivot_longer(cols = -Time, names_to = "Series", values_to = "Value")

df_MSE_clust4_low = df_clust4_low %>%  
  pivot_wider(id_cols = Time, values_from = 'Value', names_from = 'Series') %>% 
  mutate(MSE_1 = (pcp-DFM_h1)^2, MSE_12 = (pcp-DFM_h12)^2) 

RMSE_clust4_fbx$RMSE_low[1] <- sqrt(mean(df_MSE_clust4_low$MSE_1))
RMSE_clust4_fbx$RMSE_low[2] <- sqrt(mean(df_MSE_clust4_low$MSE_12, na.rm = TRUE))

rol_clust4_low <-ggplot(df_clust4_low, aes(x = Time, y = Value, color = Series, linetype = Series)) +
  geom_line(linewidth = 0.9) +
  scale_linetype_manual(values = c(pcp = "solid", DFM_h1 = "dashed", DFM_h12 = "dashed"),
                        name ='',
                        breaks = c("pcp","DFM_h1","DFM_h12"),
                        labels = labels_)+
  scale_color_manual(values = c(pcp = "black", DFM_h1= "red", DFM_h12= "purple"),
                     name ='',
                     breaks = c("pcp","DFM_h1","DFM_h12"),
                     labels = labels_) +
  labs(x = "Fecha", y = "Precipitación (mm)", title = str_to_title(names(Data2)[n_id2])) +
  coord_cartesian(ylim =c(0,400))+
  theme_bw() +
  theme(
    title = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 24),
    legend.position = c(0.2, 0.9),     # top-left inside
    legend.text = element_text(size = 24),
    legend.background = element_blank())

grid.arrange(rol_clust4_med, rol_clust4_low, ncol = 1)

names(RMSE_clust4_fbx) = c('h',str_to_title(names(Data2)[c(n_id, n_id2)]))
RMSE_clust4_fbx

xtable(
  RMSE_clust4_fbx,
  caption = "", 
  label = "tab:km4",
  align = c("l", "c", "c", "c") # column alignment: left, right, center
)

## 6.2 Grafico de las cargas -------------------------------------------------

DataLoad <- data.frame(NombreEstacion = colnames(Data2)[-1],
                       Load1 = Load[,1], Load2 = Load[,2],
                       Load3 = Load[,3], Load4 = Load[,4])

range(Load[,1]) # -0.008611787  0.084528239
range(Load[,2]) # -0.08636567  0.12061910
range(Load[,3]) # -0.1642059  0.1418977
range(Load[,4]) # -0.1498870  0.1432364
min(Load)

top_L1 =  slice_max(DataLoad, abs(Load1), n = 35, with_ties = FALSE)

probs_ = rep(0.5/nrow(DataLoad),nrow(DataLoad))
sttns_fbx = c("LUCERNA HACIENDA [26100070]", "CASABLANCA [21206060]")
ids_fbx = which(DataLoad$NombreEstacion %in% sttns_fbx)
probs_[ids_fbx] = 0.5/nrow(DataLoad) + 0.5/2

set.seed(3520)
sample_L1 = filter(DataLoad, NombreEstacion %in% 
                     sample(NombreEstacion, size = 35, replace = FALSE, prob = probs_))

sample_L1$flag <- ifelse(sample_L1$NombreEstacion %in% sttns_fbx,"highlight","normal")
table(sample_L1$flag)

Load_f1 <- ggplot(sample_L1, aes(NombreEstacion, Load1, fill = flag)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c(highlight = pal5[4], normal = "lightgrey"), guide = "none")+
  coord_flip(ylim = c(min(Load), max(Load))) +
  labs(x = "Nombre Estación", y = expression("Cargas del"~f[1])) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))

Load_f2 <- ggplot(sample_L1, aes(NombreEstacion, Load2, fill = flag)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c(highlight = pal5[4], normal = "lightgrey"), guide = "none")+
  coord_flip(ylim = c(min(Load), max(Load))) +
  labs(x = "Nombre Estación", y = expression("Cargas del"~f[2])) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))

Load_f3 <- ggplot(sample_L1, aes(NombreEstacion, Load3, fill = flag)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c(highlight = pal5[4], normal = "lightgrey"), guide = "none")+
  coord_flip(ylim = c(min(Load), max(Load))) +
  labs(x = "Nombre Estación", y = expression("Cargas del"~f[3])) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))

Load_f4 <- ggplot(sample_L1, aes(NombreEstacion, Load4, fill = flag)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c(highlight = pal5[4], normal = "lightgrey"), guide = "none")+
  coord_flip(ylim = c(min(Load), max(Load))) +
  labs(x = "Nombre Estación", y = expression("Cargas del"~f[4])) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14))
