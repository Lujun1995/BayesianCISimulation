
# 1. Generate data under linear case
gen_data = function(N, b0, b1, b2, b3) {
  # Baseline Co variate
  X1 = rnorm(N, 1, 1)
  X2 = rpois(N, 2)
  X3 = runif(N, -1, 1)
  # Sample Selection 
  ps = plogis(b0 + b1*X1 + b2*X2 + b3*X3)
  S = rbinom(N, 1, ps)
  # Treatment assignment
  A = rbinom(N, 1, 0.5)
  # Outcome
  error = rnorm(N, 0, 2)
  Y = X1 + X2 + X3 + A + A*(X1-X2) + error
  # Data
  data = data = cbind(X1, X2, X3, A, S, Y) %>% as_tibble()
  return(data)
}

gen_data2 = function(N, b0, b1, b2, b3) {
  # Baseline Co variate
  X1 = rnorm(N, 1, 1)
  X2 = rpois(N, 2)
  X3 = runif(N, -1, 1)
  # Sample Selection 
  ps = plogis(b0 + b1*X1 + b2*X2 + b3*X3)
  S = rbinom(N, 1, ps)
  # Treatment assignment
  A = rbinom(N, 1, 0.5)
  # Outcome
  error = rnorm(N, 0, 2)
  Y = X1 + X2 + X3 + A + A*(X1-X2+X3) + error
  # Data
  data = data = cbind(X1, X2, X3, A, S, Y) %>% as_tibble()
  return(data)
}

gen_data3 = function(N, b0, b1, b2, b3) {
  # Baseline Co variate
  X1 = rnorm(N, 1, 1)
  X2 = rpois(N, 2)
  X3 = runif(N, -1, 1)
  # Sample Selection 
  ps = plogis(b0 + b1*X1 + b2*X2 + b3*X3)
  S = rbinom(N, 1, ps)
  # Treatment assignment
  A = rbinom(N, 1, 0.5)
  # Outcome
  error = rnorm(N, 0, 2)
  Y = X1 + X2 + X3 + A + A*(10*X1-X2+X3) + error
  # Data
  data = data = cbind(X1, X2, X3, A, S, Y) %>% as_tibble()
  return(data)
}

# 2. Get posterior mean
norm_post = function(prior_mean = 0, prior_sd = 100, d) {
  likelihood_mean <- mean(d)
  likelihood_sd <- sd(d)
  tau <- 1 / (likelihood_sd ^ 2)
  prior_tau <- 1 / (prior_sd ^ 2)
  posterior_tau <- prior_tau + length(d) * tau
  posterior_sd <- 1 / sqrt(posterior_tau)
  posterior_mean <- (prior_mean * prior_tau + sum(d) * tau) / posterior_tau
  return(
    c(posterior_sd, 
      posterior_mean))
}

# Get unadjusted ATE
unadjusted_ATE = function(prior_mean = 0, prior_sd = 100, data) {
  selected_data = data[data$S == 1, ]
  d1 = selected_data$Y[selected_data$A == 1]
  d0 = selected_data$Y[selected_data$A == 0]
  post1 = norm_post(prior_mean = prior_mean, prior_sd = prior_sd, d1)
  post0 = norm_post(prior_mean = prior_mean, prior_sd = prior_sd, d0)
  ate_mean = post1[2] - post0[2]
  ate_lci = qnorm(0.025, ate_mean, sqrt(post1[1]^2 + post0[1]^2))
  ate_uci = qnorm(0.975, ate_mean, sqrt(post1[1]^2 + post0[1]^2))
  truth = 1 + mean(data$X1) - mean(data$X2)
  rb = ate_mean- truth
  if (truth >= ate_lci & truth<= ate_uci) {
    cover = 1
  } else {
    cover = 0
  }
  return(c(ate_mean, ate_lci, ate_uci, cover, rb, rb^2))
}

bayesian_om = function(data) {
  loaded_model <- readRDS("compiled_model.rds")
  selected_data = data[data$S == 1, ]
  input <- list(
    N1 = nrow(selected_data),
    Y = selected_data$Y,
    X1 = selected_data$X1,
    X2 = selected_data$X2,
    X3 = selected_data$X3,
    A = selected_data$A
  )
  fit = sampling(loaded_model, data = input,iter = 1500, chains = 4, warmup=500)
  post = extract(fit)
  ate_post = post$gamma + post$delta1*mean(data$X1) + post$delta2*mean(data$X2)
  ate_mean = mean(ate_post)
  ate_lci = quantile(ate_post, 0.025)
  ate_uci = quantile(ate_post, 0.975)
  truth = 1 + mean(data$X1) - mean(data$X2)
  rb = ate_mean- truth
  if (truth >= ate_lci & truth<= ate_uci) {
    cover = 1
  } else {
    cover = 0
  }
  return(c(ate_mean, ate_lci, ate_uci, cover, rb, rb^2))
}


optim_loss1 = function(data, weight) {
  data = as.matrix(data)
  # use glmnet with no penalization to find the optimal
  fit_S = glmnet::glmnet(data[ ,1:3], 
                         data[, 5], 
                         family = "binomial",
                         lambda = 0,
                         weights = weight)
  
  # calculate the opt weight for the outcome regression
  hat_ps = predict(fit_S, data[data[, 5]==1, 1:3], type="response")
  w =(1/hat_ps)*weight[data[, 5]==1]
  
  data_sample = as.data.frame(data[data[,5] ==1, ])
  model_S = lm(Y ~ A,
               data = data_sample,
               weights = w
  )
  
  return(coef(model_S)[2])
}

optim_loss2 = function(data, weight) {
  data = as.matrix(data)
  # use glmnet with no penalization to find the optimal
  fit_S = glmnet::glmnet(data[ ,1:3], 
                         data[, 5], 
                         family = "binomial",
                         lambda = 0,
                         weights = weight)
  
  # calculate the opt weight for the outcome regression
  hat_ps = predict(fit_S, data[data[, 5]==1, 1:3], type="response")
  w = 1/hat_ps*weight[data[, 5]==1]
  
  # weighted regression
  data_sample = as.data.frame(data[data[,5] ==1, ])
  model_S = lm(Y ~ X1 + X2 + X3 + A + A*X1 +  A*X2,
               data = data_sample,
               weights = w
  )
  return(est = coef(model_S)[5] + coef(model_S)[6]*mean(data[,1]) + coef(model_S)[7]*mean(data[,2]))
}

inverse_weighting = function(data) {
  ate_post = rep(0, 1000)
  N = nrow(data)
  for (j in 1:1000) {
    # dir.weig = as.numeric(rdirichlet(1, rep(1.0, N))*N)
    dir.weig = as.numeric(rmultinom(1, N, rep(1.0/N, N)))
    ate_post[j] = optim_loss1(data, dir.weig)
  }
  ate_mean = mean(ate_post)
  ate_lci = quantile(ate_post, 0.025)
  ate_uci = quantile(ate_post, 0.975)
  truth = 1 + mean(data$X1) - mean(data$X2)
  rb = ate_mean- truth
  if (truth >= ate_lci & truth<= ate_uci) {
    cover = 1
  } else {
    cover = 0
  }
  return(c(ate_mean, ate_lci, ate_uci, cover, rb, rb^2))
}
  

inverse_weighting2 = function(data) {
  ate_post = rep(0, 1000)
  N = nrow(data)
  for (j in 1:1000) {
    # dir.weig = as.numeric(rdirichlet(1, rep(1.0, N))*N)
    dir.weig = as.numeric(rmultinom(1, N, rep(1.0/N, N)))
    ate_post[j] = optim_loss2(data, dir.weig)
  }
  ate_mean = mean(ate_post)
  ate_lci = quantile(ate_post, 0.025)
  ate_uci = quantile(ate_post, 0.975)
  truth = 1 + mean(data$X1) - mean(data$X2)
  rb = ate_mean- truth
  if (truth >= ate_lci & truth<= ate_uci) {
    cover = 1
  } else {
    cover = 0
  }
  return(c(ate_mean, ate_lci, ate_uci, cover, rb, rb^2))
}


