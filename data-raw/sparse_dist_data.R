set.seed(123)

n_ids_per_group <- 200
group_levels <- c(0, 1)

id_params <- do.call(
  rbind,
  lapply(group_levels, function(group_value) {
    data.frame(
      group = group_value,
      id_within_group = seq_len(n_ids_per_group),
      n_obs = sample(3:7, n_ids_per_group, replace = TRUE)
    )
  })
)

id_params$id <- seq_len(nrow(id_params))

simulate_group_values <- function(group_value, n_obs) {
  if (group_value == 0) {
    return(rnorm(n_obs, mean = 0, sd = 1))
  }

  component <- stats::rbinom(n_obs, size = 1, prob = 0.05)
  stats::rnorm(
    n_obs,
    mean = ifelse(component == 1, 2, 0),
    sd = ifelse(component == 1, 0.25, 1)
  )
}

sparse_dist_data <- do.call(
  rbind,
  lapply(seq_len(nrow(id_params)), function(i) {
    group_value <- id_params$group[i]
    n_obs <- id_params$n_obs[i]

    data.frame(
      group = group_value,
      id = id_params$id[i],
      value = simulate_group_values(group_value, n_obs)
    )
  })
)

save(sparse_dist_data, file = "data/sparse_dist_data.rda")
