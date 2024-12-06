data {
  int<lower=1> N; // Number of observations
  array[N] int<lower=0> Y; // Observed counts
  vector[N] X_mu; // Covariate 1
  vector[N] X_sigma; // Covariate 2
}
parameters {
  real beta_0; // Intercept
  real beta_1; // Coefficient for X_mu
  real beta_2; // Coefficient for X_sigma
  real<lower=0> tau; // Standard deviation
  real<lower=0> phi; // Dispersion parameter
}
model {
  vector[N] lambda; // Rate parameter for NegBinom
  
  // Priors
  beta_0 ~ normal(0, 100); // Less informative prior
  beta_1 ~ normal(0, 100); // Less informative prior
  beta_2 ~ normal(0, 100); // Less informative prior
  tau ~ normal(0, 10); // Weakly informative half-normal
  phi ~ normal(0, 10); // Weakly informative half-normal
  
  // Likelihood
  for (n in 1 : N) {
    lambda[n] = exp(beta_0 + beta_1 * X_mu[n] + beta_2 * X_sigma[n]);
    Y[n] ~ neg_binomial_2(lambda[n], phi);
  }
}
