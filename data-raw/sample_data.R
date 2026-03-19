# Number of individuals
n_id <- 200

# Draw person-specific parameters
id_params <- data.frame(
  id = 1:n_id,
  mean = runif(n_id, min = 1, max = 4),
  sd   = runif(n_id, min = 1, max = 2),
  n_obs = sample(3:7, n_id, replace = TRUE)
)

# Simulate observations for each individual
sim_data <- do.call(rbind, lapply(1:n_id, function(i) {
  data.frame(
    id = id_params$id[i],
    value = rnorm(
      n = id_params$n_obs[i],
      mean = id_params$mean[i],
      sd = id_params$sd[i]
    )
  )
}))

sample_data <- sim_data
save(sample_data, file = "data/sample_data.rda")
