library(tidyverse)

source("R/bcd_additive.R")

stock_data <- readRDS("application/data/stock/old_data/stock_data_old.rds")

ydata <- stock_data$price_bin3

TRANSFORMER <- "sigmoid"
transformer <- pracma::sigmoid

family <- "binomial"

res2 <- bcd.ode.additive(
  ydata, seq_len(nrow(ydata[[1]])),
  family = family,
  dist_args = NULL,
  lam_th_min = 1e-0, lam_th_max = 1e-0,
  delta_update_method = "decrease",
  lam_gam_list = 10^seq(-4, 1, 0.5),
  lam_gam_ratio = "auto",
  delta_lam_gam = 10^0.2,
  knots.theta.dt = 1,
  dt.in = 0.1,
  scaling_factor = 1,
  norder_f = 4, nknots_f = 6,
  knots_position = "quantile",
  init_method = "baseline",
  lam_gam_init = 10^seq(-4, 1, 0.1),
  lam_gam_init2 = 10^seq(-4, 2, 0.2),
  penalty = "grpLasso",
  penalty.weighted = TRUE,
  smooth_package = "stats",
  theta_method = "1dH",
  initialize_penalty_weights = FALSE,
  MAXITER = 4,
  NUM_REWEIGHT = 2,
  final_result_only = FALSE,
  random_order = TRUE,
  BICobj = "ode",
  tuning_tolerant = 1,
  eps.conv = 1e-4,
  ada.gam = 1,
  L_lower_bound_factor = 1e-8,
  know_truth = FALSE,
  fit.trajectory = FALSE,
  mode.test = FALSE
)

saveRDS(
  res,
  "application/result/stock/res_stock_old_bin3.rds"
)