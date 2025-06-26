gen_irt_responses <- function(N,
                              J,
                              status,
                              p_status,
                              mean_theta_vec, sd_theta_vec,
                              a_param, b_param,
                              var_names){
  if(length(status)!=length(mean_theta_vec)) stop()
  if(length(status)!=length(sd_theta_vec)) stop()
  mean_theta_df <- as.data.frame(t(mean_theta_vec))
  sd_theta_df <- as.data.frame(t(sd_theta_vec))
  colnames(mean_theta_df) <- colnames(sd_theta_df) <-  status

  the_sample <- sample(status, N, TRUE, p_status)

  mean_theta_sample <- as.numeric(mean_theta_df[the_sample])
  sd_theta_sample <- as.numeric(sd_theta_df[the_sample])

  theta_sample <- rnorm(N, mean_theta_sample, sd_theta_sample)

  ## Generate data from a 2PLM with the model
  ## parameters 'b' and 'a' specified above
  irtmodel_prob <- function(logit) exp(logit) / (1 + exp(logit))
  init.matrix <- matrix(rep(theta_sample, J), ncol = J) #theta values for each item
  step1logit <- init.matrix - matrix(b, nrow = N, ncol = J, byrow = TRUE)
  step2logit <- step1logit * matrix(a,  nrow = N, ncol = J, byrow = TRUE)
  Pxj_th <- t(apply(step1logit, 1, irtmodel_prob)) #probabilities
  responses <- matrix(rbinom(length(Pxj_th), size = 1, prob = Pxj_th),
                      nrow = N)
  sumscores <- apply(responses, 1 ,sum)
  df <- data.frame(responses)
  colnames(df) <- paste0('Item', 1:J)
  df$Status <- the_sample
  df$sumscore <- sumscores
  return(df)
}
