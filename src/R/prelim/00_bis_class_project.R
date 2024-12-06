##' File Description
##' AUTHOR: Cole B. Brookson
##' DATE OF CREATION: 2024-11-19
#'
#' This file contains the simple analysis that will be displayed for the final
#' project in BIS567
#'
#'

#' set up ======================================================================
library(ggplot2)
source(here::here("./src/R/functions/00_global_functions.R"))
set.seed(123)

# Load necessary libraries
if (!requireNamespace("rstan", quietly = TRUE)) {
    install.packages("rstan")
}
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# temperature generation =======================================================

years <- 50
T <- years * 365 # Number of time points
trend <- seq(20, 23, length.out = T) # Increasing mean
seasonality <- 5 * sin(2 * pi * (1:T) / 180) # Seasonal component
noise <- rnorm(T, mean = 0, sd = 3) # Random noise
temperature <- trend + seasonality + noise # Multiplicative model
plot(temperature)

# Convert data into a time series object
# Assume seasonality repeats every 20 points
temperature_ts <- stats::ts(temperature, frequency = 180)

# Perform multiplicative decomposition
decomp <- stats::decompose(temperature_ts, type = "additive")

# Plot the decomposed components
plot(decomp)

# Extract components for individual plotting
trend_component <- decomp$trend
seasonal_component <- decomp$seasonal
random_component <- decomp$random

# Plot original data and components
p1 <- ggplot(
    data.frame(Time = 1:T, Temperature = temperature),
    aes(x = Time, y = Temperature)
) +
    geom_line(color = "blue") +
    labs(title = "Original Temperature Data", y = "Temperature", x = "Time") +
    theme_base()

p2 <- ggplot(
    data.frame(Time = 1:T, Trend = trend_component),
    aes(x = Time, y = Trend)
) +
    geom_line(color = "green") +
    labs(title = "Trend Component", y = "Trend", x = "Time") +
    theme_base()

p3 <- ggplot(
    data.frame(Time = 1:T, Seasonal = seasonal_component),
    aes(x = Time, y = Seasonal)
) +
    geom_line(color = "orange") +
    labs(title = "Seasonal Component", y = "Seasonal Effect", x = "Time") +
    theme_base()

p4 <- ggplot(
    data.frame(Time = 1:T, Random = random_component),
    aes(x = Time, y = Random)
) +
    geom_line(color = "purple") +
    labs(title = "Random Component", y = "Remainder", x = "Time") +
    theme_base()

# Print the plots
print(p1)
print(p2)
print(p3)
print(p4)

# Ensure patchwork is installed
if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
}

# Combine the plots using patchwork
combined_plot <- p1 + p2 + p3 + p4 +
    patchwork::plot_layout(ncol = 2) + # Arrange in a 2x2 grid
    patchwork::plot_annotation(
        # title = "Multiplicative Decomposition of Temperature Data",
        # subtitle = "Original Data, Trend, Seasonal, and Random Components",
        theme = theme_base()
    )

# Display the combined plot
print(combined_plot)

#' data simulation =============================================================

# Parameters
T <- years * 365 # Number of time points
beta_0 <- 1.0
beta_1 <- 0.9
beta_2 <- 0.3
beta_3 <- 0.5
rho <- 0.8
tau <- 0.5
phi <- 2 # Dispersion parameter

# Generate X_sigma: Sinusoidal pattern with autoregressive noise

# X_sigma <- numeric(T)
# X_sigma[1] <- rnorm(1, mean = 0, sd = 0.5)
# for (t in 2:T) {
#     X_sigma[t] <- 0.8 * X_sigma[t - 1] + sin(2 * pi * t / 20) +
#      rnorm(1, mean = 0, sd = 0.2)
# }

#' NOTE: instead of doing this with noise, I'm just going to use the data that
#' got decomposed from the timeseries
#' # taking the bit I can get the mean from
X_sigma <- seasonal_component[which(!is.na(trend_component))]

# Generate X_mu: Linear trend with autoregressive noise
# X_mu <- numeric(T)
# X_mu[1] <- 20
# for (t in 2:T) {
#     X_mu[t] <- 0.8 * X_mu[t - 1] + 0.05 * t + rnorm(1, mean = 0, sd = 0.2)
# }

#' NOTE: instead of doing this with noise, I'm just going to use the data that
#' got decomposed from the timeseries
#' # taking the bit I can get the mean from
X_mu <- trend_component[which(!is.na(trend_component))]
t_usable <- length(X_mu)
if (length(X_mu) != length(X_sigma)) {
    print(length(X_mu))
    print(length(X_sigma))
    stop("X vectors not the same length")
}
# Generate latent process U_t based on temperature T
generate_U <- function(X_mu, rho, tau, t_usable) {
    U <- numeric(t_usable)
    U[1] <- rnorm(1, mean = rho * X_mu[1], sd = tau) # Initial latent state
    for (t in 2:t_usable) {
        epsilon_t <- rnorm(1, mean = 0, sd = tau) # Stochastic component
        deterministic <- rho * X_mu[t] # Correlation with T
        stochastic <- (1 - rho) * U[t - 1] # Autoregressive component
        U[t] <- deterministic + sqrt(1 - rho^2) * epsilon_t + stochastic
    }
    return(U)
}

tau <- 1
rho <- 0.2 # Correlation between U and T
U <- generate_U(X_mu, rho, tau, t_usable)

# Plot the results
library(ggplot2)
df <- data.frame(Time = 1:t_usable, Temperature = X_mu, Latent_Process = U)
df_long <- df %>%
    tidyr::pivot_longer(
        cols = c(Temperature, Latent_Process),
        names_to = "Variable",
        values_to = "Value"
    )

# Create the plot with faceting
ggplot(df_long, aes(x = Time, y = Value)) +
    geom_line(aes(colour = Variable)) +
    facet_wrap(~Variable, scales = "free_y", ncol = 1) +
    labs(y = "Value", x = "Time") +
    theme_base()
# Calculate lambda_t (rate parameter for Negative Binomial)
lambda <- numeric(t_usable)
lambda[1] <- exp(beta_0 + beta_1 * X_mu[1] + beta_2 * X_sigma[1] + U[1])
Y <- numeric(t_usable)

# Simulate Y_{t+1} from Negative Binomial
Y[1] <- 10
# FOR NOW, we're going to scale:
X_mu <- scale(X_mu)
X_sigma <- scale(X_sigma)
for (t in 1:(t_usable - 1)) {
    lambda[t + 1] <- exp(beta_0 + beta_1 * X_mu[t] +
        beta_2 * X_sigma[t] + beta_3 * U[t])
    if (is.na(lambda[t + 1])) {
        print(t)
        print(lambda[t + 1])
        stop("NA")
    }
    size <- 1 / phi # Convert dispersion phi to size parameter
    prob <- size / (size + lambda[t + 1])
    Y[t + 1] <- rnbinom(1, size = size, prob = prob)
}

# Create data frame for plotting
sim_data <- data.frame(
    Time = 1:t_usable,
    X_mu = X_mu,
    X_sigma = X_sigma,
    U = U,
    Y = Y
)

# Plot covariates and latent process
x_mu_p <- ggplot(sim_data, aes(x = Time, y = X_mu)) +
    geom_line(color = "blue") +
    labs(title = "Covariate X_mu (Linear Trend)", y = "X_mu", x = "Time") +
    theme_base()

x_sigma_p <- ggplot(sim_data, aes(x = Time, y = X_sigma)) +
    geom_line(color = "green") +
    labs(
        title = "Covariate X_sigma (Sinusoidal Pattern)",
        y = "X_sigma", x = "Time"
    ) +
    theme_base()

u_p <- ggplot(sim_data, aes(x = Time, y = U)) +
    geom_line(color = "purple") +
    labs(title = "Latent Process U", y = "U", x = "Time") +
    theme_base()

y_p <- ggplot(sim_data, aes(x = Time, y = Y)) +
    geom_line(color = "red") +
    labs(title = "Observed Counts Y", y = "Y", x = "Time") +
    theme_base()

# Combine the plots using patchwork
response_plot <- x_mu_p + x_sigma_p + u_p + y_p +
    patchwork::plot_layout(ncol = 2) + # Arrange in a 2x2 grid
    patchwork::plot_annotation(
        # title = "Multiplicative Decomposition of Temperature Data",
        # subtitle = "Original Data, Trend, Seasonal, and Random Components",
        theme = theme_base()
    )

# Display the combined plot
print(response_plot)

#' model =======================================================================

# Simulated Data
set.seed(123)
N <- t_usable

# Prepare data for Stan
Y_stan <- as.integer(as.character(Y))
stan_data <- list(
    N = N,
    Y = Y,
    X_mu = as.vector(X_mu), # Convert to vector
    X_sigma = as.vector(X_sigma) # Convert to vector
)

# Compile Stan model
stan_model <- rstan::stan_model(
    file = here::here("./src/stan/00_negative_binomial_model.stan")
)

# Fit the model
fit <- rstan::sampling(stan_model,
    data = stan_data,
    iter = 2000,
    chains = 4,
    seed = 123
)

# Print summary of the results
print(fit, pars = c("beta_0", "beta_1", "beta_2", "tau", "phi"))

# Plot the results
rstan::traceplot(fit, pars = c("beta_0", "beta_1", "beta_2", "tau", "phi"))

# Extract posterior samples
posterior_samples <- rstan::extract(fit)

# Visualize posterior distributions
library(ggplot2)
posterior_df <- data.frame(
    beta_0 = posterior_samples$beta_0,
    beta_1 = posterior_samples$beta_1,
    beta_2 = posterior_samples$beta_2,
    tau = posterior_samples$tau,
    phi = posterior_samples$phi
)

ggplot(posterior_df, aes(x = beta_1)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = "Posterior Distribution of beta_1", x = "beta_1", y = "Density") +
    theme_minimal()

#' sensitivity analysis ========================================================

# function to run simulation and estimation for given rho
simulate_and_fit <- function(
    rho_vals, beta_3, stan_model_path,
    iterations = 2000, chains = 4) {
    results <- data.frame(rho = numeric(), beta_1_mean = numeric(), beta_1_sd = numeric())
    fit_list <- list()
    for (i in seq_along(rho_vals)) {
        rho <- rho_vals[i]
        cat(sprintf(
            "Starting simulation for rho = %.2f (%d/%d)\n",
            rho, i, length(rho_vals)
        ))

        # Simulate data
        U <- generate_U(X_mu = X_mu, rho = rho, tau = tau, t_usable = t_usable)

        lambda <- numeric(t_usable)
        lambda[1] <- exp(beta_0 + beta_1 * X_mu[1] + beta_2 * X_sigma[1] + U[1])
        Y <- numeric(t_usable)
        Y[1] <- 10

        for (t in 1:(t_usable - 1)) {
            lambda[t + 1] <- exp(beta_0 + beta_1 * X_mu[t] +
                beta_2 * X_sigma[t] + beta_3 * U[t])
            size <- 1 / phi # Convert dispersion phi to size parameter
            prob <- size / (size + lambda[t + 1])
            Y[t + 1] <- rnbinom(1, size = size, prob = prob)
        }

        # Prepare data for Stan
        Y_stan <- as.integer(as.character(Y))
        stan_data <- list(
            N = t_usable,
            Y = Y_stan,
            X_mu = as.vector(X_mu), # Convert to vector
            X_sigma = as.vector(X_sigma) # Convert to vector
        )

        # Fit the model
        cat("Compiling and sampling the Stan model...\n")
        fit <- rstan::sampling(
            stan_model,
            data = stan_data,
            iter = iterations,
            chains = chains,
            seed = 123,
            verbose = FALSE
        )
        fit_list[[paste("rho_", rho, sep = "")]] <- fit

        # Extract posterior samples
        posterior_samples <- rstan::extract(fit)
        beta_1_samples <- posterior_samples$beta_1

        # Summarize results
        beta_1_mean <- mean(beta_1_samples)
        beta_1_sd <- sd(beta_1_samples)

        cat(sprintf(
            "Finished simulation for rho = %.2f: beta_1_mean = %.4f, beta_1_sd = %.4f\n\n",
            rho, beta_1_mean, beta_1_sd
        ))

        # Append results
        results <- rbind(results, data.frame(rho = rho, beta_1_mean = beta_1_mean, beta_1_sd = beta_1_sd))
    }

    return(list(results, fit_list))
}

# explore impact of rho on beta_1 estimation
rho_vals <- seq(0, 1, by = 0.2)
beta_1_estimates <- numeric(length(rho_values))

results <- simulate_and_fit(
    rho_vals = rho_vals,
    beta_3 = beta_3,
    stan_model_path = "./src/stan/00_negative_binomial_model.stan",
    iterations = 2000,
    chains = 4
)

# compare true vs estimated beta_1
results <- data.frame(
    rho = rho_values,
    estimated_beta_1 = beta_1_estimates,
    true_beta_1 = beta_1
)

# plot results
ggplot(results, aes(x = rho)) +
    geom_line(aes(y = estimated_beta_1, color = "Estimated beta_1")) +
    geom_hline(aes(yintercept = true_beta_1, color = "True beta_1"), linetype = "dashed") +
    labs(
        title = "Impact of Correlation (rho) on beta_1 Estimation",
        x = "rho",
        y = "beta_1"
    ) +
    scale_color_manual(values = c("Estimated beta_1" = "blue", "True beta_1" = "red")) +
    theme_minimal()


ggplot(results, aes(x = rho, y = beta_1_mean)) +
    geom_line(color = "blue") +
    geom_point(size = 2) +
    labs(
        title = "Mean of Posterior Samples for beta_1 Across rho",
        x = "rho",
        y = "Mean of beta_1"
    ) +
    theme_minimal()

# Plot beta_1 standard deviation vs. rho
ggplot(results, aes(x = rho, y = beta_1_sd)) +
    geom_line(color = "red") +
    geom_point(size = 2) +
    labs(
        title = "Standard Deviation of Posterior Samples for beta_1 Across rho",
        x = "rho",
        y = "Standard Deviation of beta_1"
    ) +
    theme_minimal()
