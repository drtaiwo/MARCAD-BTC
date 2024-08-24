# Load necessary libraries
library(deSolve)  # For solving differential equations
library(ggplot2)  # For plotting
library(gridExtra)  # For arranging multiple plots
library(reshape2)  # For reshaping data

# Function to calculate the seasonal variation in transmission
# time: Time variable
# amplitude: Amplitude of the seasonal variation
# phase_shift: Shift in the phase of the seasonal cycle
seasonal_modulation <- function(time, amplitude, phase_shift) {
  return(1 + amplitude * cos(2 * pi * (time - phase_shift) / 365))
}

# Function to calculate the effective reduction in transmission due to mosquito nets
# N_net: Number of mosquito nets distributed
# u_net: Usage rate of distributed nets (0-1)
# k_net: Effectiveness of nets (0-1, where 1 is 100% effective)
# Nh: Total human population
calculate_net_effect <- function(N_net, u_net, k_net, Nh) {
  net_coverage = (N_net / Nh) * u_net
  effective_reduction = net_coverage * k_net
  return(effective_reduction)
}

# Define the system of differential equations (SEIRV model)
# time: Current time in the simulation
# state: Current state of the system (populations in each compartment)
# parameters: Parameters of the model
disease_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Calculate the reduction in transmission due to nets if intervention is active
    if (time >= net_start_time && use_net_intervention) {
      effective_reduction <- calculate_net_effect(N_net, u_net, k_net, N_h)
    } else {
      effective_reduction <- 0
    }
    
    # Calculate seasonal modulation of transmission
    seas <- seasonal_modulation(time, amplitude, phase_shift)
    
    # Force of infection for humans
    FOIh <- (beta_hv * I_v / N_h) * seas * (1 - effective_reduction)
    
    # Human population equations
    dS_h <- -FOIh * S_h + delta * R  # Susceptible humans
    dE_h <- FOIh * S_h - psi * (alpha_s + alpha_c * (1 - alpha_s) + alpha_a * (1 - alpha_c) * (1 - alpha_s) + (1 - alpha_a) * (1 - alpha_c) * (1 - alpha_s)) * E_h  # Exposed humans
    dI_s <- psi * alpha_s * E_h + sigma_s * R + mu_c * I_c - (gamma_s + theta_s) * I_s  # Symptomatic infected humans
    dI_c <- alpha_c * (1 - alpha_s) * psi * E_h + sigma_c * R + gamma_s * I_s + mu_a * I_a - (theta_c + gamma_c + mu_c) * I_c  # Clinically infected humans
    dI_a <- alpha_a * (1 - alpha_c) * (1 - alpha_s) * psi * E_h + gamma_c * I_c + mu_n * I_n + sigma_a * R - (mu_a + gamma_a) * I_a  # Asymptomatic infected humans
    dI_n <- psi * (1 - alpha_a) * (1 - alpha_c) * (1 - alpha_s) * E_h + gamma_a * I_a + sigma_n * R - (mu_n + kappa_n) * I_n  # Non-specific infected humans
    dR <- kappa_s * T_s + kappa_c * T_c + kappa_n * I_n - (delta + sigma_s + sigma_c + sigma_a + sigma_n) * R  # Recovered humans
    dT_s <- theta_s * I_s - kappa_s * T_s  # Treated symptomatic humans
    dT_c <- theta_c * I_c - kappa_c * T_c  # Treated clinically infected humans
    
    # Vector population equations
    FOIv <- ((I_c + I_a + I_s + I_n) / N_v) * beta_vh * seas * (1 - effective_reduction)
    dS_v <- -FOIv * S_v  # Susceptible vectors
    dE_v <- FOIv * S_v - omega * E_v  # Exposed vectors
    dI_v <- omega * E_v  # Infected vectors
    
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
  beta_hv = 0.9,        # Transmission rate from vectors to humans: Probability of transmission per contact from a vector to a human.
  beta_vh = 0.5,        # Transmission rate from humans to vectors: Probability of transmission per contact from a human to a vector.
  
  alpha_s = 0.5,        # Proportion of exposed individuals progressing to symptomatic infection (I_s).
  alpha_c = 0.5,        # Proportion of exposed individuals progressing to clinical infection (I_c).
  alpha_a = 0.5,        # Proportion of exposed individuals progressing to asymptomatic infection (I_a).
  
  sigma_s = 0.01,       # Rate at which symptomatic infections progress to recovery or death.
  sigma_c = 0.01,       # Rate at which clinical infections progress to recovery or death.
  sigma_a = 0.01,       # Rate at which asymptomatic infections progress to recovery or death.
  sigma_n = 0.01,       # Rate at which infection in vectors (I_n) progresses.
  
  gamma_s = 0.1,        # Rate of recovery or death from symptomatic infections (I_s).
  gamma_c = 0.1,        # Rate of recovery or death from clinical infections (I_c).
  gamma_a = 0.1,        # Rate of recovery or death from asymptomatic infections (I_a).
  
  theta_s = 0.05,       # Rate of treatment or intervention for symptomatic infections (I_s).
  theta_c = 0.05,       # Rate of treatment or intervention for clinical infections (I_c).
  
  mu_n = 0.01,          # Mortality rate for vectors (I_n).
  mu_c = 0.01,          # Mortality rate for clinical infections (I_c).
  mu_a = 0.01,          # Mortality rate for asymptomatic infections (I_a).
  kappa_n = 0.01,       # Rate at which vectors (I_n) are removed from the population.
  
  kappa_s = 0.5,        # Clearance rate for symptomatic infections (T_s).
  kappa_c = 0.5,        # Clearance rate for clinical infections (T_c).
  
  delta = 0.05,         # Natural death rate for humans.
  
  omega = 0.1,          # Transmission rate from vector to vector.
  
  N_h = 1000,           # Total number of humans in the population.
  N_v = 1000,           # Total number of vectors in the population.
  
  psi = 0.2,            # Probability of transmission per contact between humans and vectors.
  
  amplitude = 0.8,      # Amplitude of the seasonal variation in transmission.
  
  phase_shift = 180,    # Phase shift of the seasonal variation in transmission, in days.
  
  N_net = 6000,         # Total number of mosquito nets available.
  u_net = 0.5,          # Usage rate of mosquito nets (fraction of nets used).
  k_net = 0.2,          # Effectiveness of mosquito nets in reducing transmission.
  
  net_start_time = 1    # Time (in days) when mosquito nets are introduced into the model.
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

rowSums(output_with_nets[,2:10])
rowSums(output_with_nets[,11:13])
