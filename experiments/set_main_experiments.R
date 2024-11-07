d_t <- 20 / N
# consider N=20,40,100,200

if (family == "g") {
  # R <- 1
  num_rep <- 1
  init_method <- "baseline"
  initw <- TRUE
  MAXITER <- 4
  lam_gam_list <- 10^(seq(-6, 2, 0.25))
  lam_gam_list.grade <- 10^seq(-6, 2, 0.1)
  lam_gam_list.saode <- 10^seq(-6, 2, 0.1)
  lam_gam_list2.saode <- 10^seq(-6, 2, 0.2)
  lam_gam_init <- 10^seq(-6, 1, 0.1)
  lam_gam_init2 <- 10^seq(-6, 2, 0.2)
  all_grpreg <- TRUE
  # if (SNR > 10) {
  #   all_grpreg <- TRUE
  # } else {
  #   all_grpreg <- FALSE
  # }
} else if (family == "p") {
  num_rep <- 1
  init_method <- "baseline"
  initw <- FALSE
  MAXITER <- 4
  lam_gam_list <- 10^(seq(-6, 2, 0.25))
  lam_gam_list.grade <- 10^seq(-6, 2, 0.1)
  lam_gam_list.saode <- 10^seq(-6, 2, 0.1)
  lam_gam_list2.saode <- 10^seq(-6, 2, 0.2)
  lam_gam_init <- 10^seq(-6, 1, 0.1)
  lam_gam_init2 <- 10^seq(-6, 2, 0.2)
  all_grpreg <- TRUE
} else if (family == "b") {
  num_rep <- 1
  init_method <- "baseline"
  initw <- FALSE
  MAXITER <- 4
  lam_gam_list <- 10^(seq(-6, 2, 0.25))
  lam_gam_list.grade <- 10^seq(-6, 2, 0.1)
  lam_gam_list.saode <- 10^seq(-6, 2, 0.1)
  lam_gam_list2.saode <- 10^seq(-6, 2, 0.2)
  lam_gam_init <- 10^seq(-6, 1, 0.1)
  lam_gam_init2 <- 10^seq(-6, 2, 0.2)
  all_grpreg <- TRUE
}