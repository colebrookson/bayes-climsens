data {
  int<lower=1> N;              // Number of observations
  int<lower=0> Y[N];           // Observed counts
  vector[N] X_mu;              // Covariate 1
  vector[N] X_sigma;           // Covariate 2
}
parameters {
  real beta_0;                 // Intercept
  real beta_1;                 // Coefficient for X_mu
  real beta_2;                 // Coefficient for X_sigma
  vector[N] U;                 // Latent process
  real<lower=0> tau;           // Standard deviation of U
  real<lower=0> phi;           // Dispersion parameter
}
model {
  vector[N] lambda;            // Rate parameter for NegBinom
  
  // Priors
  beta_0 ~ normal(0, 10);
  beta_1 ~ normal(0, 10);
  beta_2 ~ normal(0, 10);
  tau ~ normal(0, 2);
  phi ~ normal(0, 2);
  
  // Latent process prior (AR(1))
  U[1] ~ normal(0, tau);
  for (n in 2:N) {
    U[n] ~ normal(U[n - 1], tau);
  }
  
  // Likelihood
  for (n in 1:N) {
    lambda[n] = exp(beta_0 + beta_1 * X_mu[n] + beta_2 * X_sigma[n] + U[n]);
    Y[n] ~ neg_binomial_2(lambda[n], phi);
  }
}
