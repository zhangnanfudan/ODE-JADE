# aggregate estimation result


library(tidyverse)


## helpers =============

# numerical tools
estAUC <- function(tpr, fpr) {
  # tpr, fpr of an decreasng order of lambdas
  ord <- order(fpr)
  tpr <- tpr[ord]
  fpr <- fpr[ord]
  if (min(tpr) != 0) {
    tpr <- c(0, tpr)
    fpr <- c(0, fpr)
  }
  if (max(tpr != 1)) {
    # worst case assumption
    tpr <- c(tpr, max(tpr))
    fpr <- c(fpr, 1)
  }
  daf <- data.frame(tpr = tpr, fpr = fpr) %>%
    group_by(fpr) %>%
    summarise(meantpr = mean(tpr))
  return(sum(0.5 * diff(daf$fpr) * (daf$meantpr[-1] + daf$meantpr[-nrow(daf)])))
}

# format tools
sciNote <- function(x) {
  if (x == 0) {
    return("0")
  } else {
    expo <- as.integer(floor(log10(x)))
    coeff <- round(x / (10^expo), 3)
    if (is.na(coeff)) browser()
    if (expo == 0) {
      return(as.character(coeff))
    } else {
      return(paste0(coeff, "\\times ", "10^{", expo, "}"))
    }
  }
}

div10Expo <- function(x, expo = 0, digits = 3) {
  return(round(x / 10^expo, digits = digits))
}

format_num <- function(x, expo = 0, digits = 3, num_float = 3) {
  y <- div10Expo(x, expo = expo, digits = digits)
  is_zero <- (y == 0)
  is_float <- (y %% 1 != 0)
  if (is_zero) {
    return(stringi::stri_pad_right(paste0("0."), max(num_float + 2, 1), 0))
  }
  num_int <- max(floor(log10(y)) + 1, 1)
  ndec <- max(num_int - 1, 0)
  if (is_float) {
    y <- as.character(y)
  } else {
    y <- paste0(as.character(y), ".")
  }
  pad_width <- num_int + 1 + num_float - ndec
  return(stringi::stri_pad_right(y, pad_width, 0))
}

makecell <- function(x, y, pt_reduce = 0) {
  # add alignment
  if (pt_reduce == 0) {
    shrink <- ""
  } else {
    shrink <- paste0("[-", pt_reduce, "pt]")
  }
  return(paste0("\\bigcell{r}{ ", x, " \\\\", shrink, " ", " (", y, ") }"))
}

add.ampersand <- function(x) {
  nx <- length(x)
  x[-nx] <- paste0(x[-nx], " &")
  return(x)
}

add.phantom <- function(x, loc = NULL) {
  nx <- length(x)
  if (is.null(loc)) {
    loc <- rep(1, nx - 1)
  }
  # loc=c(2,2,3,3,3)
  loc <- cumsum(loc + 1)
  y <- character(nx + length(loc))
  y[loc] <- "\\phantom{abc}"
  y[-loc] <- x
  return(y)
}

# combine.mean_and_sd <- function(x) {
#   x <- matrix(as.character(x), nrow = 2)
#   x[2,] <- paste0("(",x[2,],")")
#   x <- apply(x, 2, paste, collapse = " ")
#   return(x)
# }

combine.mean_and_sd <- function(x, pt_reduce = 0) {
  x <- matrix(as.character(x), nrow = 2)
  y <- makecell(x[1, ], x[2, ], pt_reduce = pt_reduce)
  return(y)
}


generate.cmidrule <- function(width, first_loc = 1, phantom_columns = TRUE) {
  # \cmidrule{8-9}
  nc <- length(width)
  if (phantom_columns) {
    start_loc <- c(0, cumsum(width[-nc] + 1)) + first_loc
  } else {
    start_loc <- c(0, cumsum(width[-nc])) + first_loc
  }
  end_loc <- start_loc + width - 1
  out <- paste0(r"{\cmidrule{}", start_loc, "-", end_loc, "}")
  return(out)
}


## aggregate Gaussian results =====================

filenames <- paste0("result/", c(
  "poisson_a_N200R10",
  "poisson_a_N100R10",
  "poisson_a_N40R10",
  "binomial_a_N200R40",
  "binomial_a_N100R40",
  "binomial_a_N40R40"
), ".txt")

nN <- 3
nF <- 2
nset <- nN * nF

N <- rep(c(200, 100, 40), nF)
Fam <- rep(c("Poisson", "Bernoulli"), each = nN)


dat <- lapply(
  seq_along(filenames),
  function(i) {
    dat <- read.delim(filenames[i], header = TRUE, sep = ",")
    dat["Fam"] <- Fam[i]
    dat["N"] <- N[i]
    return(dat)
  }
) %>% do.call(rbind, .)

selected_methods <- c("jointest", "GRADE", "SAODE")
nmethod <- length(selected_methods)

tab <- dat %>%
  mutate(
    method = factor(method, levels = c("jointest", "GRADE", "SAODE")),
    Fam = factor(Fam, levels = c("Poisson", "Bernoulli"))
  ) %>%
  dplyr::filter(method %in% selected_methods) %>%
  group_by(Fam, N, method) %>%
  summarise(
    theta.mse.mean = mean(theta.mse),
    theta.mse.sd = sd(theta.mse),
    dtheta.mse.mean = mean(dtheta.mse),
    dtheta.mse.sd = sd(dtheta.mse),
    f.mse.mean = mean(f.mse.t),
    f.mse.sd = sd(f.mse.t),
    f.mse.t.mean = mean(f.mse.t),
    f.mse.t.sd = sd(f.mse.t),
    f.mse.f.mean = mean(f.mse.f),
    f.mse.f.sd = sd(f.mse.f),
    ftheta.mse.mean = mean(ftheta.mse.t),
    ftheta.mse.sd = sd(ftheta.mse.t),
    ftheta.mse.t.mean = mean(ftheta.mse.t),
    ftheta.mse.t.sd = sd(ftheta.mse.t),
    ftheta.mse.f.mean = mean(ftheta.mse.f),
    ftheta.mse.f.sd = sd(ftheta.mse.f),
    tpr.mean = mean(tpr),
    tpr.sd = sd(tpr),
    fpr.mean = mean(fpr),
    fpr.sd = sd(fpr)
  ) %>%
  arrange(Fam, desc(N)) %>%
  ungroup()


# calculate aucs
# filenames2 <- paste0("result/roc/", c(
#   "poisson_a_N40R5",
#   "poisson_a_N20R5",
#   "binomial_a_N40R40",
#   "binomial_a_N20R40"
# ),"roc.txt")


# dat2 <- lapply(seq_along(filenames2),
#                function(i) {
#                  res <- read.delim(filenames2[i],header = TRUE,sep=",")
#                  ids <- cumsum(c(1,!(res$method[-nrow(res)]==res$method[-1])*1))
#                  nexp <- max(unique(ids))/3
#                  if (nexp%%1!=0) stop("something wrong")
#                  tmp <- data.frame(
#                    method = rep(c("jointest","GRADE","SAODE"), nexp),
#                    Fam = rep(Fam[i],nexp*3), N = rep(N[i],nexp*3)
#                  )
#                  aucs <- sapply(
#                    unique(ids),
#                    function(id) {
#                      estAUC(res$tpr[ids==id], res$fpr[ids==id])
#                    }
#                  )
#                  tmp$auc <- aucs
#                  return(tmp)
#                }) %>% do.call(rbind, .)
#
#
# tab2 <- dat2 %>%
#   mutate(method = factor(method, levels = c("jointest","GRADE","SAODE"))) %>%
#   dplyr::filter(method %in% selected_methods) %>%
#   group_by(Fam, N, method) %>%
#   summarise(auc.mean = mean(auc),
#             auc.sd = sd(auc)) %>%
#   arrange(Fam, desc(N)) %>%
#   ungroup() %>%
#   mutate(
#     auc.mean = round(auc.mean, digits = 3),
#     auc.sd = round(auc.sd, digits = 3)
#   )



# crit_names <- c("theta.mse", "dtheta.mse", "f.mse.t", "f.mse.f", "tpr", "fpr")
crit_names <- c("theta.mse", "dtheta.mse", "ftheta.mse.t", "ftheta.mse.f", "tpr", "fpr")
ncrit <- length(crit_names)
nsetcol <- 3
var_names <- c(
  rbind(
    paste0(crit_names, ".mean"),
    paste0(crit_names, ".sd")
  )
)
set_names <- paste0("Fam", Fam, "N", N)
expos <- c(-2, 0, 0, 0, 0, 0)
ndigits <- numfloats <- c(3, 3, 3, 3, 3, 3)


# tab <- left_join(tab, tab2, by = c("Fam", "N", "method"))
tab <- tab %>%
  select(all_of(c("Fam", "N", "method", var_names))) %>%
  relocate(method, .after = N) %>%
  mutate(
    method = case_when(
      method == "jointest" ~ "JADE",
      method == "GRADE" ~ "GRADE",
      method == "SAODE" ~ "SA-ODE"
    )
  )
tab <- tab %>%
  relocate(Fam, N, method, .before = Fam)
for (i in seq_along(crit_names)) {
  for (ii in 1:2) {
    tab[, nsetcol + 2 * i - 2 + ii] <- sapply(tab[[nsetcol + 2 * i - 2 + ii]],
      format_num,
      expo = expos[i],
      digits = ndigits[i],
      num_float = numfloats[i]
    )
  }
}
tab <- tab %>% mutate(Fam = as.character(Fam))

latex_simple <- c()
for (i in seq_len(nrow(tab))) {
  rowval <- as.character(tab[i, ])
  rowlatex <- c(
    as.character(rowval[1]), rowval[2:3],
    combine.mean_and_sd(rowval[4:length(rowval)],
      pt_reduce = 7
    )
  )
  if (rowval[3] == "JADE") {
    rowlatex[1:2] <- paste0(
      r"{\multirow{3}{*}{ }", rowlatex[1:2], r"{ }}"
    )
  }
  if (rowval[3] == "GRADE") {
    rowlatex[1:2] <- " "
    rowlatex[4:5] <- paste0(
      r"{\multirow{2}{*}{ }", rowlatex[4:5], r"{ }}"
    )
  }
  if (rowval[3] == "SA-ODE") {
    rowlatex[1:2] <- " "
    rowlatex[4:5] <- " "
  }
  # rowlatex <- add.phantom(rowlatex, loc = c(4,1,1,1))
  rowlatex <- add.ampersand(rowlatex)
  rowlatex <- paste(paste(rowlatex, collapse = " "), r"{\\}", "\n")
  latex_simple <- c(latex_simple, rowlatex)
  if (i %% 3 == 0) {
    latex_simple <- c(latex_simple, paste0(r"{\cmidrule{1-}", ncrit + nsetcol, r"{}}", " \n"))
    latex_simple <- c(latex_simple, "\n")
  }
}

cat(latex_simple)
