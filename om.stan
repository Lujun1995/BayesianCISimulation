data {
  int<lower=0> N1;         // Number of observations
  vector[N1] Y;            // Response variable
  vector[N1] X1;           // Predictor variable X1
  vector[N1] X2;           // Predictor variable X2
  vector[N1] X3;           // Predictor variable X3
  vector[N1] A;            // Additional predictor (A)
}

parameters {
  real alpha;             // Intercept
  real beta1;            // Coefficient for X1
  real beta2;            // Coefficient for X2
  real beta3;            // Coefficient for X3
  real gamma;            // Coefficient for A
  real delta1;           // Coefficient for the interaction A*X1
  real delta2;           // Coefficient for the interaction A*X2
  real<lower=0> sigma;   // Standard deviation of the residuals
}

model {
  // Prior distributions
  alpha ~ normal(0, 10);   // Prior for the intercept
  beta1 ~ normal(0, 10);   // Prior for the coefficient for X1
  beta2 ~ normal(0, 10);   // Prior for the coefficient for X2
  beta3 ~ normal(0, 10);   // Prior for the coefficient for X3
  gamma ~ normal(0, 10);   // Prior for the coefficient for A
  delta1 ~ normal(0, 10);  // Prior for the coefficient for A*X1 interaction
  delta2 ~ normal(0, 10);  // Prior for the coefficient for A*X2 interaction
  sigma ~ cauchy(0, 5);    // Prior for the standard deviation
  // Likelihood
  Y ~ normal(alpha + beta1 * X1 + beta2 * X2 + beta3 * X3 + gamma * A  + delta2 * A .* X2, sigma);
}



