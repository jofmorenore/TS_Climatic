
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

sarima(fac1, p = 0, d = 0, q = 0, P = 0, D = 1, Q = 1, S = 12,
       no.constant = T)
SARIMA_f1_m2 <- forecast::Arima(fac1, order = c(0, 0, 0), 
                                seasonal = list(order = c(0, 1, 1), period = 12))
summary(SARIMA_f1_m2) # sar1 No significativo
checkresiduals(SARIMA_f1_m2)

Box.test(SARIMA_f1_m2$residuals, lag = (length(SARIMA_f1_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 2)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f1_m2$residuals) # pvalue = 2.034e-05 < 0.05 (R. H0, No concluye normalidad)
cusum_g(SARIMA_f1_m2$residuals)

## 4.4.2 Factor 2 --------------------------------- 
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

sarima(fac2, p = 0, d = 0, q = 0, P = 0, D = 1, Q = 1, S = 12,
       no.constant = T)
SARIMA_f2_m2 <- forecast::Arima(fac2, order = c(0, 0, 0), 
                                seasonal = list(order = c(0, 1, 1), period = 12))
summary(SARIMA_f2_m2) # sar1 No significativo
checkresiduals(SARIMA_f2_m2)

Box.test(SARIMA_f2_m2$residuals, lag = (length(SARIMA_f2_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 2)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f2_m2$residuals) # pvalue = 2.034e-05 < 0.05 (R. H0, No concluye normalidad)
cusum_g(SARIMA_f2_m2$residuals)

## 4.4.3 Factor 3 --------------------------------- 
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

sarima(fac3, p = 0, d = 0, q = 0, P = 0, D = 1, Q = 1, S = 12,
       no.constant = T)
SARIMA_f3_m2 <- forecast::Arima(fac3, order = c(0, 0, 0), 
                                seasonal = list(order = c(0, 1, 1), period = 12))
summary(SARIMA_f3_m2) # sar1 No significativo
checkresiduals(SARIMA_f3_m2)

Box.test(SARIMA_f3_m2$residuals, lag = (length(SARIMA_f3_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 2)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(SARIMA_f3_m2$residuals) # pvalue = 2.034e-05 < 0.05 (R. H0, No concluye normalidad)
cusum_g(SARIMA_f3_m2$residuals)

## 4.4.4 Factor 4 --------------------------------- 
sarima(fac4, p = 2, d = 1, q = 1, P = 0, D = 0, Q = 0, S = 0,
       no.constant = T)
ARIMA_f4_m1 <- forecast::Arima(fac4, order = c(2, 1, 1), 
                                seasonal = list(order = c(0, 0, 0), period = 0))
summary(ARIMA_f4_m1) # coeficientes significativos
checkresiduals(ARIMA_f4_m1)

Box.test(ARIMA_f4_m1$residuals, lag = (length(ARIMA_f4_m1$residuals)/4), 
         type = "Ljung-Box", fitdf = 4)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(ARIMA_f4_m1$residuals) # pvalue = 0.8861 > 0.05 (No R. H0, concluye normalidad)

cusum_g(ARIMA_f4_m1$residuals)

sarima(fac4, p = 3, d = 1, q = 3, P = 0, D = 0, Q = 0, S = 0,
       no.constant = T)

ARIMA_f4_m2 <- forecast::Arima(fac4, order = c(3, 1, 3), 
                                seasonal = list(order = c(0, 0, 0), period = 0))
summary(ARIMA_f4_m2) # sar1 No significativo
checkresiduals(ARIMA_f4_m2)

Box.test(ARIMA_f4_m2$residuals, lag = (length(ARIMA_f4_m2$residuals)/4), 
         type = "Ljung-Box", fitdf = 2)
# p-value < 0.05 (R. H0, concluye significant autocorrelation detected)

tseries::jarque.bera.test(ARIMA_f4_m2$residuals) # pvalue = 0.007278 < 0.05 (R. H0, No concluye normalidad)
cusum_g(ARIMA_f4_m2$residuals)

SARIMA_f1_m2$arma
