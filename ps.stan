data {
  int<lower=0> N; // Number of observations
  int<lower=0,upper=1> S[N];  // Binary outcome variable
  vector[N] X1;           // Predictor variable X1
  vector[N] X2;           // Predictor variable X2
  vector[N] X3;           // Predictor variable X3
}


parameters {
  real theta0;            // Intercept
  real theta1;            // Coefficient for X1
  real theta2;            // Coefficient for X2
  real theta3;            // Coefficient for X3
}


model {
  // Likelihood
  for (i in 1:N) {
    S[i] ~ bernoulli_logit(theta0 + theta1 * X1[i] + theta2 * X2[i] + theta3*X3[i]);
  }
}

