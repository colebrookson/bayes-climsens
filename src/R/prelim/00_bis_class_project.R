#' set up ======================================================================
library(ggplot2)
set.seed(123)

T <- 100 # Number of time points
trend <- seq(10, 20, length.out = T) # Increasing mean
seasonality <- 2 * sin(2 * pi * (1:T) / 20) # Seasonal component
noise <- rnorm(T, mean = 1, sd = 0.2) # Random noise
temperature <- trend * seasonality * noise # Multiplicative model

# Convert data into a time series object
# Assume seasonality repeats every 20 points
temperature_ts <- ts(temperature, frequency = 20)

# Perform multiplicative decomposition
decomp <- stats::decompose(temperature_ts, type = "multiplicative")

# Plot the decomposed components
plot(decomp)

# Extract components for individual plotting
trend_component <- decomp$trend
seasonal_component <- decomp$seasonal
random_component <- decomp$random

# Plot original data and components
p1 <- ggplot(data.frame(Time = 1:T, Temperature = temperature), aes(x = Time, y = Temperature)) +
    geom_line(color = "blue") +
    labs(title = "Original Temperature Data", y = "Temperature", x = "Time") +
    theme_minimal()

p2 <- ggplot(data.frame(Time = 1:T, Trend = trend_component), aes(x = Time, y = Trend)) +
    geom_line(color = "green") +
    labs(title = "Trend Component", y = "Trend", x = "Time") +
    theme_minimal()

p3 <- ggplot(data.frame(Time = 1:T, Seasonal = seasonal_component), aes(x = Time, y = Seasonal)) +
    geom_line(color = "orange") +
    labs(title = "Seasonal Component", y = "Seasonal Effect", x = "Time") +
    theme_minimal()

p4 <- ggplot(data.frame(Time = 1:T, Random = random_component), aes(x = Time, y = Random)) +
    geom_line(color = "purple") +
    labs(title = "Random Component", y = "Remainder", x = "Time") +
    theme_minimal()

# Print the plots
print(p1)
print(p2)
print(p3)
print(p4)
#' data simulation =============================================================

# Parameters
T <- 100 # Number of time points
beta_0 <- 1.0
beta_1 <- 0.5
beta_2 <- -0.3
rho <- 0.8
tau <- 0.5
phi <- 2 # Dispersion parameter

# Generate X_sigma: Sinusoidal pattern with autoregressive noise
X_sigma <- numeric(T)
X_sigma[1] <- rnorm(1, mean = 0, sd = 0.5)
for (t in 2:T) {
    X_sigma[t] <- 0.8 * X_sigma[t - 1] + sin(2 * pi * t / 20) + rnorm(1, mean = 0, sd = 0.2)
}

# Generate X_mu: Linear trend with autoregressive noise
X_mu <- numeric(T)
X_mu[1] <- 20
for (t in 2:T) {
    X_mu[t] <- 0.8 * X_mu[t - 1] + 0.05 * t + rnorm(1, mean = 0, sd = 0.2)
}

# Generate latent process U_t
U <- numeric(T)
U[1] <- rnorm(1, mean = 0, sd = tau) # Initial latent state
for (t in 2:T) {
    U[t] <- rnorm(1, mean = rho * U[t - 1], sd = tau)
}

# Calculate lambda_t (rate parameter for Negative Binomial)
lambda <- numeric(T)
lambda[1] <- exp(beta_0 + beta_1 * X_mu[1] + beta_2 * X_sigma[1] + U[1])
Y <- numeric(T)

# Simulate Y_{t+1} from Negative Binomial
Y[1] <- 100
for (t in 1:(T - 1)) {
    lambda[t + 1] <- exp(beta_0 + beta_1 * X_mu[t] + beta_2 * X_sigma[t] + U[t])
    size <- 1 / phi # Convert dispersion phi to size parameter
    prob <- size / (size + lambda[t + 1])
    Y[t + 1] <- rnbinom(1, size = size, prob = prob)
}

# Create data frame for plotting
sim_data <- data.frame(
    Time = 1:T,
    X_mu = X_mu,
    X_sigma = X_sigma,
    U = U,
    Y = Y
)

# Plot covariates and latent process
p1 <- ggplot(sim_data, aes(x = Time, y = X_mu)) +
    geom_line(color = "blue") +
    labs(title = "Covariate X_mu (Linear Trend)", y = "X_mu", x = "Time") +
    theme_minimal()

p2 <- ggplot(sim_data, aes(x = Time, y = X_sigma)) +
    geom_line(color = "green") +
    labs(title = "Covariate X_sigma (Sinusoidal Pattern)", y = "X_sigma", x = "Time") +
    theme_minimal()

p3 <- ggplot(sim_data, aes(x = Time, y = U)) +
    geom_line(color = "purple") +
    labs(title = "Latent Process U", y = "U", x = "Time") +
    theme_minimal()

p4 <- ggplot(sim_data, aes(x = Time, y = Y)) +
    geom_line(color = "red") +
    labs(title = "Observed Counts Y", y = "Y", x = "Time") +
    theme_minimal()

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)
