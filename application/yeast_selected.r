library(tidyverse)

data_name <- "yeast_data_selected"

source("R/bcd_additive.R")

yeast_data <- readRDS(paste0(
  "application/data/yeast_cell_cycle/", data_name, ".rds"
))

TRANSFORMER <- "sigmoid"
transformer <- pracma::sigmoid

family <- "gaussian"


res <- bcd.ode.additive(
  list(yeast_data$ydata),
  seq_along(yeast_data$times),
  family = family,
  dist_args = NULL,
  lam_th_min = 1e0, lam_th_max = 1e0,
  delta_update_method = "decrease",
  lam_gam_list = 10^seq(-5, 2, 0.5),
  lam_gam_ratio = "auto",
  knots.theta.dt = 1,
  dt.in = 0.1,
  scaling_factor = 1,
  norder_f = 4, nknots_f = 8,
  knots_position = "quantile",
  init_method = "baseline",
  lam_gam_init = 10^seq(-4, 1, 0.1),
  lam_gam_init2 = 10^seq(-4, 1, 0.2),
  penalty = "grpLasso",
  penalty.weighted = TRUE,
  smooth_package = "stats",
  theta_method = "1dH",
  initialize_penalty_weights = TRUE,
  MAXITER = 4,
  NUM_REWEIGHT = 2,
  tuning_tolerant = 1,
  random_order = TRUE,
  BICobj = "ode",
  eps.conv = 1e-4,
  ada.gam = 1,
  know_truth = FALSE,
  fit.trajectory = FALSE,
  mode.test = FALSE
)

saveRDS(
  res,
  paste0("application/result/yeast_cell_cycle/res_", data_name, ".rds")
)