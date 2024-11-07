if (N==20) {
  d_t <- 1
} else if (N==40) {
  d_t <- 0.5
} else {
  stop("Haven't consider N other than 20, 40")
}
if (family=="p") {
  if (R != 5) stop("Is R wrong?")
}
if (family=="b") {
  if (!(R %in% c(10,20,40))) stop("Is R wrong?")
}

lam_gam_list = 10^seq(-5,3,0.25)
lam_gam_list0 = 10^seq(-2,2,0.2)  # no use
init_method <- "GRADE"
initw_method <- "jade"
if (family=="b" || family=="p") {
  lam_gam_init = 10^seq(-3,2,0.1)
  lam_gam_init2 = 10^seq(-3,4,0.2)
} else if (family=="g") {
  lam_gam_init = 10^seq(-2,2,0.1)
  lam_gam_init2 = 10^seq(-3,4,0.2)
}
lam_gam_list.grade = 10^seq(-6,3,0.1)
lam_gam_list.saode = 10^seq(-6,2,0.1)
lam_gam_list2.saode = 10^seq(-6,4,0.2)
MAXITER = 5
MAXITER0 = 0