set.seed(123)

n_ids_per_group <- 200
group_levels <- c(1, 2)

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

group1_idx <- id_params$group == 1
group2_idx <- id_params$group == 2

id_params$mean_bulk <- NA_real_
id_params$sd_bulk <- NA_real_
id_params$tail_prob <- NA_real_
id_params$tail_mean <- NA_real_
id_params$tail_sd <- NA_real_

# Group 1 remains centered near N(0, 1), but each id has its own
# location/scale perturbation.
id_params$mean_bulk[group1_idx] <- stats::rnorm(sum(group1_idx), mean = 0, sd = 0.35)
id_params$sd_bulk[group1_idx] <- exp(stats::rnorm(sum(group1_idx), mean = 0, sd = 0.15))

# Group 2 keeps the same overall bimodal structure, with id-specific
# variation in the bulk component and the tail component.
id_params$mean_bulk[group2_idx] <- stats::rnorm(sum(group2_idx), mean = 0, sd = 0.35)
id_params$sd_bulk[group2_idx] <- exp(stats::rnorm(sum(group2_idx), mean = 0, sd = 0.15))
id_params$tail_prob[group2_idx] <- stats::rbeta(sum(group2_idx), shape1 = 2, shape2 = 38)
id_params$tail_mean[group2_idx] <- stats::rnorm(sum(group2_idx), mean = 2, sd = 0.2)
id_params$tail_sd[group2_idx] <- exp(stats::rnorm(
  sum(group2_idx),
  mean = log(0.25),
  sd = 0.12
))

simulate_group_values <- function(group_value,
                                  n_obs,
                                  mean_bulk,
                                  sd_bulk,
                                  tail_prob = NA_real_,
                                  tail_mean = NA_real_,
                                  tail_sd = NA_real_) {
  if (group_value == 1) {
    return(stats::rnorm(n_obs, mean = mean_bulk, sd = sd_bulk))
  }

  component <- stats::rbinom(n_obs, size = 1, prob = tail_prob)
  stats::rnorm(
    n_obs,
    mean = ifelse(component == 1, tail_mean, mean_bulk),
    sd = ifelse(component == 1, tail_sd, sd_bulk)
  )
}

sparse_dist_data <- do.call(
  rbind,
  lapply(seq_len(nrow(id_params)), function(i) {
    row_i <- id_params[i, ]

    data.frame(
      group = row_i$group,
      id = row_i$id,
      value = simulate_group_values(
        group_value = row_i$group,
        n_obs = row_i$n_obs,
        mean_bulk = row_i$mean_bulk,
        sd_bulk = row_i$sd_bulk,
        tail_prob = row_i$tail_prob,
        tail_mean = row_i$tail_mean,
        tail_sd = row_i$tail_sd
      )
    )
  })
)

# Hand-tune the first three ids so the plotting examples have stable,
# easy-to-interpret shapes.
example_ids <- data.frame(
  group = c(
    rep(1, 3),
    rep(1, 5),
    rep(1, 6)
  ),
  id = c(
    rep(1, 3),
    rep(2, 5),
    rep(3, 6)
  ),
  value = c(
    rep(0.2, 3),
    -0.35, -0.12, 0.08, 0.21, 3.25,
    -3.1, -1.8, -0.5, 0.6, 1.9, 3.4
  )
)

sparse_dist_data <- rbind(
  subset(sparse_dist_data, !(id %in% 1:3)),
  example_ids
)

sparse_dist_data <- sparse_dist_data[order(sparse_dist_data$id), ]
row.names(sparse_dist_data) <- NULL

save(sparse_dist_data, file = "data/sparse_dist_data.rda")
