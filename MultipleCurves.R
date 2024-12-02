####################################################################################################################
####################################################################################################################
#------------------------------------- Packages --------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots
library(Quandl) ### Data
library(dplyr) ### Data manipulations
library(stats) ### ACF plots
library(matrixcalc) ### Matrix calculations
library("RColorBrewer") ### Colors
library(latex2exp) ### Text for plots
library(matrixStats) ### ColSds


####################################################################################################################
####################################################################################################################
#----------------------------------- Working Directory -------------------------------------------------------------
####################################################################################################################
####################################################################################################################


setwd("~/Documents/KU/Seminar - Asset Prices and Financial Markets/RCode")


####################################################################################################################
####################################################################################################################
#------------------------------- Parameters & Time Discretization --------------------------------------------------
####################################################################################################################
####################################################################################################################
### Vasicek
kappa <- 0.25       ### Short Rate Speed of mean reversion
theta <- 0.03       ### Short Rate Long-term mean
sigma <- 0.002      ### Short Rate Volatility
r0 <- 0.025         ### initial short rate
### Time discretization
tau <- 0.25         ### 3M tenor
T_maturities <- seq(0.25, 30, by=0.25) ### Maturities from 3M to 30Y
### Basis Factor
eta <- 0.007  ### Volatility for the chi process
### Spread
alpha <- -0.5  ### Impact of forward rate movements on spread ~ Correlation
beta <- 0.1   ### Impact of chi process on spread ~ Volatility
S0 <- 0.0025   ### Initial spread value
chi0 <- 1



####################################################################################################################
####################################################################################################################
#----------------------------------- Functions ---------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Function to simulate Vasicek short rate path
### Simulate short rate using Euler-Maruyama method for the Vasicek model
VasicekShortRate <- function(kappa, theta, sigma, r0, t_grid) {
  dt <- diff(c(0, t_grid))  ### Time steps
  r <- numeric(length(t_grid))
  r[1] <- r0
  
  for (i in 1:length(t_grid)) {
    r[i + 1] <- r[i] + kappa * (theta - r[i]) * dt[i] + sigma * sqrt(dt[i]) * rnorm(1)
  }
  return(r)
}

### Function to compute A(t,T) for the Vasicek model
A_Vasicek <- function(t, T, kappa, theta, sigma) {
  B_tT <- (1 - exp(-kappa * (T - t))) / kappa
  term1 <- (theta - (sigma^2) / (2 * kappa^2)) * (B_tT - (T - t))
  term2 <- (sigma^2) / (4 * kappa) * B_tT^2
  return(term1 - term2)
}

### Function to compute B(t,T) for the Vasicek model
B_Vasicek <- function(t, T, kappa) {
  (1 - exp(-kappa * (T - t))) / kappa
}

### Function to compute zero-coupon bond price p_D(t,T) under Vasicek
p_D <- function(t, T, r_t, kappa, theta, sigma) {
  A_tT <- A_Vasicek(t, T, kappa, theta, sigma)
  B_tT <- B_Vasicek(t, T, kappa)
  return(exp(A_tT - B_tT * r_t))
}

### Compute OIS forward rates
OISForward <- function(t, T_k, r_t, kappa, theta, sigma, tau) {
  P_T_k_minus_1 <- p_D(t, T_k - tau, r_t, kappa, theta, sigma)
  P_T_k <- p_D(t, T_k, r_t, kappa, theta, sigma)
  return((1 / tau) * (P_T_k_minus_1 / P_T_k - 1))
}


####################################################################################################################
####################################################################################################################
#------------------------------------ OIS + LIBOR ------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Seed
set.seed(13)
### Time discretization
t_grid <- seq(0, 30, by = tau)  ### Discretize time up to 30 years

### Initialize matrices to store OIS and LIBOR rates
OIS <- matrix(0, nrow = length(t_grid), ncol = length(T_maturities))
LIBOR <- matrix(0, nrow = length(t_grid), ncol = length(T_maturities))

### Simulate chi process (Brownian motion)
Z <- rnorm(length(t_grid))  ### Brownian increments
chi <- exp(-0.5 * (eta^2) * t_grid + eta * Z * sqrt(t_grid))

# Loop over the grid of t-values
for (t_index in 1:length(t_grid)) {
  t <- t_grid[t_index]
  
  # Simulate short rate using Vasicek model up to the current t
  VasiSim <- VasicekShortRate(kappa, theta, sigma, r0, t_grid[1:t_index])
  
  # Loop over maturities and compute forward OIS and LIBOR rates
  for (k in 1:length(T_maturities)) {
    T_k <- T_maturities[k]
    
    # Forward OIS rate at time t
    F_0 <- OISForward(0, T_k, VasiSim[1], kappa, theta, sigma, tau)
    F_t <- OISForward(t, T_k, VasiSim[t_index], kappa, theta, sigma, tau)
    
    # Store OIS rate
    OIS[t_index, k] <- F_t
    
    # Compute spread
    S_t <- S0 + alpha * (F_t - F_0) + beta * (chi[t_index] - chi0)
    
    # Forward LIBOR rate
    LIBOR[t_index, k] <- F_t + S_t
  }
}


t_fixed <- 3.0  # Choose a specific time
plot(T_maturities, OIS[which(t_grid == t_fixed),], type = "l", col = "#901a1E", lwd = 2,
     ylab = "Forward Rates", xlab = "Maturities (years)",
     main = paste("Forward OIS and LIBOR Curves at t = ", t_fixed),
     ylim = c(min(OIS[which(t_grid == t_fixed),]), max(LIBOR[which(t_grid == t_fixed),])))
lines(T_maturities, LIBOR[which(t_grid == t_fixed),], col = "#666666", lwd = 2)
legend("bottomright", legend=c("Forward OIS", "Forward LIBOR"), col=c("#901a1E", "#666666"), lty=1, lwd=2)


max(LIBOR[which(t_grid == t_fixed),]-OIS[which(t_grid == t_fixed),])
### alpha / beta  0.001               0.01              0.1
### -0.1       0.00285087        0.002702928      0.001223512     
### -0.5       0.004320099        0.004172158      0.002692742
min(LIBOR[which(t_grid == t_fixed),]-OIS[which(t_grid == t_fixed),])
### alpha / beta  0.001               0.01              0.1
### -0.1       0.002483772        0.002335831      0.0008564146     
### -0.5       0.002484614        0.002336672      0.0008572562
mean(LIBOR[which(t_grid == t_fixed),]-OIS[which(t_grid == t_fixed),])
### alpha / beta  0.001               0.01              0.1
### -0.1       0.002533394        0.002385452      0.0009060359     
### -0.5       0.002732721        0.002584779      0.001105363




t_fixed <- 3.0  # Choose a specific time
plot(T_maturities, OIS[which(t_grid == t_fixed),], type = "l", col = "#901a1E", lwd = 2,
     ylab = "Forward Rates", xlab = "Maturities (years)",
     main = paste("Forward OIS and LIBOR Curves"),
     ylim = c(min(OIS[which(t_grid == t_fixed),]), max(LIBOR[which(t_grid == t_fixed),])))
lines(T_maturities, LIBOR[which(t_grid == t_fixed),], col = "#666666", lwd = 2)
lines(T_maturities, OIS[which(t_grid == 0),], col = "#901a1E", lwd = 2, lty = 2)
lines(T_maturities, LIBOR[which(t_grid == 0),], col = "#666666", lwd = 2, lty = 2)
legend("bottomright", legend=c("Forward OIS at t = 3", "Forward LIBOR at t = 3", "Forward OIS at t = 0", "Forward LIBOR at t = 0"), col=c("#901a1E", "#666666","#901a1E", "#666666"), lty=c(1,1,2,2), lwd=2)



T_k_fixed <- 30.0  # Choose a specific maturity
plot(t_grid, OIS[,which(T_maturities == T_k_fixed)], type = "l", 
     ylab = expression(F[30 * Y](t) ~ "&" ~ L[30 * Y](t)), xlab = "t", 
     main = "OIS Forward Rate F_30Y(t) & LIBOR Forward Rate L_30Y(t)")
lines(t_grid, LIBOR[,which(T_maturities == T_k_fixed)], col = "blue")
legend("bottomleft", legend=c("Forward OIS", "Forward LIBOR"), col=c("black", "blue"), lty=1)







### Plot OIS and LIBOR curves for each t on the grid
par(mfrow=c(1,1))  ### 2 plots vertically

### Plot OIS rates for different t-values
matplot(T_maturities, t(OIS), type = "l", col = rainbow(length(t_grid)),
        ylab = "OIS Forward Rate", xlab = "Maturity (years)",
        main = "OIS Forward Rates as a Function of Maturity for Different t-values")
legend("topright", legend = paste0("t = ", round(t_grid, 2)), col = rainbow(length(t_grid)), lty = 1)

### Plot LIBOR rates for different t-values
matplot(T_maturities, t(LIBOR), type = "l", col = rainbow(length(t_grid)),
        ylab = "LIBOR Forward Rate", xlab = "Maturity (years)",
        main = "LIBOR Forward Rates as a Function of Maturity for Different t-values")
legend("topright", legend = paste0("t = ", round(t_grid, 2)), col = rainbow(length(t_grid)), lty = 1)


####################################################################################################################
####################################################################################################################
#--------------------------------- OIS Forward Curve ---------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Time discretization
dt <- 0.01
t_grid_Fine <- seq(0, 30, by = dt)  ### Discretize time up to 30 years

### Simulate the short rate over time
VasiSimFine <- VasicekShortRate(kappa, theta, sigma, r0, t_grid_Fine)

### Compute forward rates for different maturities
OISFine <- sapply(T_maturities, function(T_k) {
  sapply(t_grid_Fine, function(t) {
    OISForward(t, T_k, VasiSimFine[floor(t / dt) + 1], kappa, theta, sigma, tau)
  })
})

### Plot forward rates as a function of t for a fixed maturity T_k (say T_k = 2.0)
plot(t_grid_Fine, OISFine[,which(T_maturities == 30.0)], type = "l", 
     ylab = expression(F[30*Y](t)), xlab = "t", 
     main = "OIS Forward Rate F_30Y(t)")

### Plot forward rates as a function of maturity T_k for a fixed t (say t = 10.0)
plot(T_maturities, OISFine[which(t_grid_Fine == 10.0),], type = "l", 
     ylab = "rate", xlab = "T", 
     main = "Forward Rates vs Maturity (t = 10.0)")



####################################################################################################################
####################################################################################################################
#---------------------------- Spread + LIBOR as function of t ------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Simulate basis factor process for t_grid using geometric Brownian motion
ChiFine <- exp(-0.5 * (eta^2) * t_grid_Fine + eta * cumsum(rnorm(length(t_grid_Fine))) * sqrt(dt))

### Function to compute the spread S_k^Delta(t)
spread_S <- function(t, F_t, F_0, chi_t, chi_0, alpha, beta, S_0) {
  S_0 + alpha * (F_t - F_0) + beta * (chi_t - chi_0)
}

### Compute forward LIBOR rates for different maturities
LIBORFine <- sapply(T_maturities, function(T_k) {
  F_0 <- OISForward(0, T_k, VasiSimFine[1], kappa, theta, sigma, tau)  # F_k^Delta(0)
  sapply(t_grid_Fine, function(t) {
    F_t <- OISForward(t, T_k, VasiSimFine[floor(t / dt) + 1], kappa, theta, sigma, tau)
    S_t <- spread_S(t, F_t, F_0, ChiFine[floor(t / dt) + 1], chi0, alpha, beta, S0)
    return(F_t + S_t)
  })
})

### Plot the forward OIS curve and LIBOR curve for a fixed maturity T_k (say T_k = 1.0)
T_k_fixed <- 30.0  # Choose a specific maturity
plot(t_grid_Fine, OISFine[,which(T_maturities == T_k_fixed)], type = "l", col = "black", 
     ylab = "Forward Rates", xlab = "t", 
     main = paste("Forward Rates for Maturity T_k =", T_k_fixed))
lines(t_grid_Fine, LIBORFine[,which(T_maturities == T_k_fixed)], col = "blue")
legend("topright", legend=c("Forward OIS", "Forward LIBOR"), col=c("black", "blue"), lty=1)



####################################################################################################################
####################################################################################################################
#---------------------------- Spread + LIBOR as function of T ------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Set a fixed time t (e.g., t = 10)
t_fix <- 10

### Initialize empty vectors to store the forward OIS and LIBOR rates for each maturity
OIS_t <- numeric(length(T_maturities))
LIBOR_t <- numeric(length(T_maturities))

### Simulate chi process at t_fixed (assuming chi_Delta(0) = 1)
chi_t <- exp(-0.5 * (eta^2) * t_fix + eta * rnorm(1) * sqrt(t_fix))

### Loop through maturities and compute forward OIS and LIBOR rates
for (i in 1:length(T_maturities)) {
  T_k <- T_maturities[i]
  
  ### Compute forward OIS rate F_k^Delta(t_fixed)
  F_0 <- OISForward(0, T_k, VasiSimFine[1], kappa, theta, sigma, tau)  ### Forward rate at t=0
  F_t <- OISForward(t_fix, T_k, VasiSimFine[floor(t_fix / dt) + 1], kappa, theta, sigma, tau)  ### At t_fixed
  
  ### Store forward OIS rate
  OIS_t[i] <- F_t
  
  ### Compute spread S_k^Delta(t_fixed)
  S_t <- S0 + alpha * (F_t - F_0) + beta * (chi_t - chi0)
  
  ### Compute forward LIBOR rate
  LIBOR_t[i] <- F_t + S_t
}

### Plot the forward OIS and LIBOR curves as a function of maturities
plot(T_maturities, OIS_t, type = "l", col = "black", lwd = 2,
     ylab = "Forward Rates", xlab = "Maturities (years)",
     main = paste("Forward OIS and LIBOR Curves at t =", t_fix))
lines(T_maturities, LIBOR_t, col = "blue", lwd = 2)
legend("bottomright", legend=c("Forward OIS", "Forward LIBOR"), col=c("black", "blue"), lty=1, lwd=2)





