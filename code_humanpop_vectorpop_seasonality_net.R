library(deSolve)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Function to calculate the seasonal variation in transmission
seasonal_modulation <- function(time, amplitude, phase_shift) {
  return(1 + amplitude * sin(2 * pi *(time - phase_shift) / 365))
}

# Function to calculate the effective reduction in transmission due to mosquito nets
calculate_net_effect <- function(N_net, u_net, k_net, Nh) {
  net_coverage <- (N_net / Nh) * u_net
  effective_reduction <- net_coverage * k_net
  return(effective_reduction)
}

# Define the system of differential equations
disease_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    if (time >= net_start_time && use_net_intervention) {
      effective_reduction <- calculate_net_effect(N_net, u_net, k_net, N_h)
    } else {
      effective_reduction <- 0
    }
    
    seas <- seasonal_modulation(time, amplitude, phase_shift)
    FOIh <- (beta_hv * I_v / N_h) * seas * (1 - effective_reduction)
    
    # Human population equations
    dS_h <- -FOIh * S_h + delta * R
    dE_h <- FOIh * S_h - psi * (alpha_s + alpha_c * (1 - alpha_s) + alpha_a * (1 - alpha_c) * (1 - alpha_s) + (1 - alpha_a) * (1 - alpha_c) * (1 - alpha_s)) * E_h
    dI_s <- psi * alpha_s * E_h + sigma_s * R + mu_c * I_c - (gamma_s + theta_s) * I_s
    dI_c <- alpha_c * (1 - alpha_s) * psi * E_h + sigma_c * R + gamma_s * I_s + mu_a * I_a - (theta_c + gamma_c + mu_c) * I_c
    dI_a <- alpha_a * (1 - alpha_c) * (1 - alpha_s) * psi * E_h + gamma_c * I_c + mu_n * I_n + sigma_a * R - (mu_a + gamma_a) * I_a
    dI_n <- psi * (1 - alpha_a) * (1 - alpha_c) * (1 - alpha_s) * E_h + gamma_a * I_a + sigma_n * R - (mu_n + kappa_n) * I_n
    dR <- kappa_s * T_s + kappa_c * T_c + kappa_n * I_n - (delta + sigma_s + sigma_c + sigma_a + sigma_n) * R
    dT_s <- theta_s * I_s - kappa_s * T_s
    dT_c <- theta_c * I_c - kappa_c * T_c
    
    # Vector population equations
    FOIv <- ((I_c + I_a + I_s + I_n) / N_v) * beta_vh * seas * (1 - effective_reduction)
    dS_v <- -FOIv * S_v  
    dE_v <- FOIv * S_v - omega * E_v
    dI_v <- omega * E_v
    
    # Return the rates of change
    list(c(dS_h, dE_h, dI_s, dI_c, dI_a, dI_n, dR, dT_s, dT_c, dS_v, dE_v, dI_v))
  })
}

# Initial conditions
initial_state <- c(
  S_h = 9990, E_h = 6, I_s = 1, I_c = 1, I_a = 1, I_n = 1, R = 0, T_s = 0, T_c = 0,
  S_v = 500, E_v = 10, I_v = 5
)

# Parameters
parameters <- c(
  beta_hv = 0.9,
  beta_vh = 0.5,
  alpha_s = 0.5, alpha_c = 0.5, alpha_a = 0.5, 
  sigma_s = 0.01, sigma_c = 0.01, sigma_a = 0.01, sigma_n = 0.01,
  gamma_s = 0.1, gamma_c = 0.1, gamma_a = 0.1,
  theta_s = 0.05, theta_c = 0.05,
  mu_n = 0.01, mu_c = 0.01, mu_a = 0.01, kappa_n = 0.01,
  kappa_s = 0.5, kappa_c = 0.5,
  delta = 0.05,
  omega = 0.1,
  N_h = 1000,
  N_v = 1000,
  psi = 0.2,
  amplitude = 0.3,
  phase_shift = 180,
  N_net = 6000,
  u_net = 0.5,
  k_net = 0.2,
  net_start_time =1
)

# Time sequence
time <- seq(0, 300, by = 1)

# Solve the model with nets
parameters["use_net_intervention"] <- TRUE
output_with_nets <- ode(y = initial_state, times = time, func = disease_model, parms = parameters)

# Solve the model without nets
parameters["use_net_intervention"] <- FALSE
output_without_nets <- ode(y = initial_state, times = time, func = disease_model, parms = parameters)

# Convert outputs to data frames for plotting
output_with_nets_df <- as.data.frame(output_with_nets)
output_with_nets_df$scenario <- "With Nets"
output_without_nets_df <- as.data.frame(output_without_nets)
output_without_nets_df$scenario <- "Without Nets"

# Combine the data frames
combined_output_df <- rbind(output_with_nets_df, output_without_nets_df)

# Melt the data for ggplot2
output_long <- melt(combined_output_df, id.vars = c("time", "scenario"))

# Define compartments in the same order as in the model
compartments <- c("S_h", "E_h", "I_s", "I_c", "I_a", "I_n", "R", "T_s", "T_c", "S_v", "E_v", "I_v")

# Create a list to store individual plots
plot_list <- list()

# Generate plots for each compartment
for (compartment in compartments) {
  p <- ggplot(subset(output_long, variable == compartment), aes(x = time, y = value, color = scenario)) +
    geom_line() +
    labs(title = compartment, x = "Time", y = "Value") +
    theme_minimal()
  plot_list[[compartment]] <- p
}

# Arrange plots in a 4x3 grid
grid.arrange(grobs = plot_list, ncol = 3, nrow = 4)

rowSums(output_with_nets[,2:11])