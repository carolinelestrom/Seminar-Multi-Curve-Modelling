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
### Basis Factor
eta <- 0.007  ### Volatility for the chi process
### Spread
alpha <- -0.05  ### Impact of forward rate movements on spread ~ Correlation
beta <- 0.07   ### Impact of chi process on spread ~ Volatility
S0 <- 0.0025   ### Initial spread value
chi0 <- 1


####################################################################################################################
####################################################################################################################
#------------------------------------ Parameters -------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Simulate the trajectories and compute exposure
n_trajectories <- 100   ### Number of simulation paths
swap_horizon <- 3       ### Swap horizon in years (3Y payer swap)
tau_FLOAT <- 0.25       ### 3M tenor for floating leg
tau_FIX <- 0.5          ### 6M tenor for fixed leg
N <- 1e6                ### Notional
t_grid <- seq(0, swap_horizon, by = 1/360)  ### Time discretization
T_FLOAT <- seq(0.0, 3, by=tau_FLOAT)         ### Floating maturities from 3M to 3Y
T_FIX <- seq(0.0, 3, by=tau_FIX)         ### Fixed maturities from 3M to 3Y




####################################################################################################################
####################################################################################################################
#---------------------------------- Par Swap Rate ------------------------------------------------------------------
####################################################################################################################
####################################################################################################################


### Function to compute the par swap rate (K_par)
ParSwapRate <- function(t, T_FLOAT, T_FIX, tau_FLOAT, tau_FIX, VasiSim, kappa, theta, sigma, S0, alpha, beta, chi_t) {
  numerator <- 0
  denominator <- 0
  
  ### Loop over maturities to compute numerator (floating leg)
  for (k in 1:length(T_FLOAT)) {
    T_k <- T_FLOAT[k]
    
    # Ensure that the current time t is less than the maturity T_k
    if (t >= T_k) {
      next  # Skip this iteration if the condition is not met
    }
    
    ### Compute discount factor for floating leg
    P_T_k <- p_D(t, T_k, VasiSim, kappa, theta, sigma)
    
    ### Compute forward OIS and LIBOR rates
    F_0 <- OISForward(0, T_k, VasiSim, kappa, theta, sigma, tau_FLOAT)
    F_t <- OISForward(t, T_k, VasiSim, kappa, theta, sigma, tau_FLOAT)
    
    S_t <- S0 + alpha * (F_t - F_0) + beta * (chi_t - chi0)
    LIBOR_t <- F_t + S_t
    
    ### Numerator for the floating leg (forward LIBOR)
    numerator <- numerator + tau_FLOAT * P_T_k * LIBOR_t
  }
  ### Loop over maturities to compute denominator (fixed leg)
  for (j in 1:length(T_FIX)) {
    T_j <- T_FIX[j]
    
    # Ensure that the current time t is less than the maturity T_k
    if (t >= T_j) {
      next  # Skip this iteration if the condition is not met
    }
    
    ### Compute discount factor for floating leg
    P_T_j <- p_D(t, T_j, VasiSim, kappa, theta, sigma)
    
    ### Denominator for the fixed leg (discount factor only)
    denominator <- denominator + tau_FIX * P_T_j
  }
  
  ### Par swap rate
  return(numerator / denominator)
}

SwpRate <- numeric(length(t_grid))
### Loop over time steps to compute exposures
for (t_index in 1:length(t_grid)) {
  t <- t_grid[t_index]
  
  ### Compute the par swap rate K_par at each time t
  SwpRate[t_index] <- ParSwapRate(t_grid[t_index], T_FLOAT, T_FIX, tau_FLOAT, tau_FIX, VasiSim[t_index], kappa, theta, sigma, S0, alpha, beta, chi[t_index])
}

plot(SwpRate)

####################################################################################################################
####################################################################################################################
#--------------------------------- PV of Plain Vanilla Swap --------------------------------------------------------
####################################################################################################################
####################################################################################################################


### Initialize matrices to store exposures
PV_FLOAT <- matrix(0, nrow = n_trajectories, ncol = length(t_grid))
PV_FIX <- matrix(0, nrow = n_trajectories, ncol = length(t_grid))
K_par <- numeric(length(t_grid))

### Run 1000 simulations
for (n in 1:n_trajectories) {
  ### Simulate chi process (Brownian motion)
  Z <- rnorm(length(t_grid))  ### Brownian increments
  chi <- exp(-0.5 * (eta^2) * t_grid + eta * Z * sqrt(t_grid))
  
  ### Simulate short rate using Vasicek model
  VasiSim <- VasicekShortRate(kappa, theta, sigma, r0, t_grid)
  
  ### Loop over time steps to compute exposures
  for (t_index in 1:length(t_grid)) {
    t <- t_grid[t_index]
    
    ### Compute the par swap rate K_par at each time t
    K_par[t_index] <- ParSwapRate(t_grid[t_index], T_FLOAT, T_FIX, tau_FLOAT, tau_FIX, VasiSim[t_index], kappa, theta, sigma, S0, alpha, beta, chi[t_index])
    
    
    for (k in 1:length(T_FLOAT)) {
      T_k <- T_FLOAT[k]
      
      # Ensure that the current time t is less than the maturity T_k
      if (t >= T_k) {
        next  # Skip this iteration if the condition is not met
      }
      
      ### Compute OIS and LIBOR forward rates as before
      F_0 <- OISForward(0, T_k, VasiSim[1], kappa, theta, sigma, tau_FLOAT)
      F_t <- OISForward(t, T_k, VasiSim[t_index], kappa, theta, sigma, tau_FLOAT)
      
      S_t <- S0 + alpha * (F_t - F_0) + beta * (chi[t_index] - chi0)
      
      ### Compute LIBOR forward rate
      LIBOR_t <- F_t + S_t
      
      ### Store present values for floating leg (LIBOR leg)
      PV_FLOAT[n, t_index] <- PV_FLOAT[n, t_index] +
        tau_FLOAT * p_D(t, T_k, VasiSim[t_index], kappa, theta, sigma) * LIBOR_t * N
    }
    
    for (j in 1:length(T_FIX)) {
      T_j <- T_FIX[j]
      
      # Ensure that the current time t is less than the maturity T_j
      if (t >= T_j) {
        next  # Skip this iteration if the condition is not met
      }
      
      ### Store present value of the fixed leg using K_par
      PV_FIX[n, t_index] <- PV_FIX[n, t_index] +
        tau_FIX * p_D(t, T_j, VasiSim[t_index], kappa, theta, sigma) * K_par[t_index] * N
    }
  }
}

### Compute net exposure for each path and each time step
PV_Payer <- PV_FLOAT - PV_FIX

mean <- q5 <- q95 <- numeric(length(t_grid))
for (i in 1:length(t_grid)){
  mean[i] <- mean(PV_Payer[,i])
  q5[i] <- quantile(PV_Payer[,i], probs = 0.05)
  q95[i] <- quantile(PV_Payer[,i], probs = 0.95)
}

### Plot PV trajectories over time
plot(t_grid, PV_Payer[1,], type = "l", col = 1, lwd = 2,
     xlab = "Time (years)", ylab = "Present Value (USD)",
     main = "1000 Simulated Present Value Paths of a 3Y Payer Plain Vanilla Swap",
     #ylim = c(-200000, 200000))
     ylim = c(min(PV_Payer), max(PV_Payer)))
for (i in 2:n_trajectories){
  lines(t_grid, PV_Payer[i,], t = 'l', ylim = c(min(PV_Payer), max(PV_Payer)), col = i)
}


### Plot exposure profile over time
plot(t_grid, mean, type = "l", col = "#901a1E", lwd = 2,
     xlab = "Time (years)", ylab = "Present Value (USD)",
     main = "Exposure Profile of 3Y Payer Plain Vanilla Swap",
     #ylim = c(-200000, 200000))
     ylim = c(min(q5), max(q95)))
lines(t_grid, q5, t = 'l', ylim = c(min(q5), max(q95)), col = 1, lwd = 3)
lines(t_grid, q95, t = 'l', ylim = c(min(q5), max(q95)), col = 1, lwd = 3)


