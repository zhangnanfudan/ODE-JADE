# aggregate estimation result


library(tidyverse)


## helpers =============

# numerical tools
estAUC <- function(tpr, fpr) {
  # tpr, fpr of an decreasng order of lambdas
  ord <- order(fpr)
  tpr <- tpr[ord]; fpr <- fpr[ord]
  if (min(tpr) !=0 ) {
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
  if (x==0) {
    return("0")
  } else {
    expo <- as.integer(floor(log10(x)))
    coeff <- round(x/(10^expo),3)
    if (is.na(coeff)) browser()
    if (expo==0) {
      return(as.character(coeff))
    } else {
      return(paste0(coeff, "\\times ", "10^{",expo,"}"))
    }
  }
}

div10Expo <- function(x, expo=0, digits=3) {
  return(round(x/10^expo, digits = digits))
}

format_num <- function(x, expo=0, digits=3, num_float=3) {
  y <- div10Expo(x, expo = expo, digits = digits)
  is_zero <- (y == 0)
  is_float <- (y%%1 != 0)
  if (is_zero) {
    return(stringi::stri_pad_right(paste0("0."), max(num_float+2,1), 0))
  }
  num_int <- max(floor(log10(y))+1, 1)
  ndec <- max(num_int - 1, 0)
  if (is_float) {
    y <- as.character(y)
  } else {
    y <- paste0(as.character(y),".")
  }
  pad_width <- num_int + 1 + num_float - ndec
  return(stringi::stri_pad_right(y, pad_width, 0))
}

makecell <- function(x, y, pt_reduce = 0) {
  # add alignment
  if (pt_reduce==0) {
    shrink <- ""
  } else {
    shrink <- paste0("[-",pt_reduce,"pt]")
  }
  return(paste0("\\bigcell{r}{ ", x," \\\\",shrink," "," (", y,") }"))
}

add.ampersand <- function(x) {
  nx <- length(x)
  x[-nx] <- paste0(x[-nx], " &")
  return(x)
}

add.phantom <- function(x, loc = NULL) {
  nx <- length(x)
  if (is.null(loc)) {
    loc <- rep(1,nx-1)
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
  y <- makecell(x[1,], x[2,], pt_reduce = pt_reduce)
  return(y)
}


generate.cmidrule <- function(width, first_loc = 1, phantom_columns = TRUE) {
  # \cmidrule{8-9}
  nc <- length(width)
  if (phantom_columns) {
    start_loc <- c(0, cumsum(width[-nc]+1)) + first_loc
  } else {
    start_loc <- c(0, cumsum(width[-nc])) + first_loc
  }
  end_loc <- start_loc + width - 1
  out <- paste0(r"{\cmidrule{}",start_loc,"-",end_loc,"}")
  return(out)
}


## aggregate Gaussian results =====================

filenames <- paste0("result/", c(
  "gaussian_a_N40R1SNR4",
  "gaussian_a_N40R1SNR10",
  "gaussian_a_N40R1SNR25",
  "gaussian_a_N100R1SNR4",
  "gaussian_a_N100R1SNR10",
  "gaussian_a_N100R1SNR25",
  "gaussian_a_N200R1SNR4",
  "gaussian_a_N200R1SNR10",
  "gaussian_a_N200R1SNR25"
),".txt")

nN <- 3
nSNR <- 3
nset <- nN * nSNR

N <- rep(c(40,100,200), each=nSNR)
R <- rep(1,nN*nSNR)
SNR <- rep(c(4,10,25),nN)


dat <- lapply(seq_along(filenames),
       function(i) {
         dat <- read.delim(filenames[i],header = TRUE,sep=",")
         dat['N'] <- N[i]
         dat['R'] <- R[i]
         dat['SNR'] <- SNR[i]
         return(dat)
       }) %>% do.call(rbind, .)

selected_methods <- c("jointest","GRADE","SAODE")
nmethod <- length(selected_methods)

tab <- dat %>%
  mutate(method = factor(method, levels = c("jointest","GRADE","SAODE"))) %>% 
  dplyr::filter(method %in% selected_methods) %>% 
  group_by(method, N, R, SNR) %>% 
  summarise(theta.mse.mean = mean(theta.mse),
            theta.mse.sd = sd(theta.mse),
            dtheta.mse.mean = mean(dtheta.mse),
            dtheta.mse.sd = sd(dtheta.mse),
            f.mse.mean = mean(f.mse),
            f.mse.sd = sd(f.mse),
            f.mse.t.mean = mean(f.mse.t),
            f.mse.t.sd = sd(f.mse.t),
            f.mse.f.mean = mean(f.mse.f),
            f.mse.f.sd = sd(f.mse.f),
            ftheta.mse.mean = mean(ftheta.mse),
            ftheta.mse.sd = sd(ftheta.mse),
            ftheta.mse.t.mean = mean(ftheta.mse.t),
            ftheta.mse.t.sd = sd(ftheta.mse.t),
            ftheta.mse.f.mean = mean(ftheta.mse.f),
            ftheta.mse.f.sd = sd(ftheta.mse.f),
            tpr.mean = mean(tpr),
            tpr.sd = sd(tpr),
            fpr.mean = mean(fpr),
            fpr.sd = sd(fpr)) %>% 
  arrange(desc(N), R, desc(SNR)) %>% 
  mutate(N = N * R) %>%
  ungroup %>% 
  dplyr::select(-R)

# calculate aucs
# filenames2 <- paste0("result/roc/", c(
#   "gaussian_a_N20R1SNR4",
#   "gaussian_a_N20R1SNR10",
#   "gaussian_a_N40R1SNR4",
#   "gaussian_a_N40R1SNR10"
# ),"roc.txt")
# 
# dat2 <- lapply(seq_along(filenames2),
#               function(i) {
#                 res <- read.delim(filenames2[i],header = TRUE,sep=",")
#                 ids <- cumsum(c(1,!(res$method[-nrow(res)]==res$method[-1])*1))
#                 nexp <- max(unique(ids))/3
#                 if (nexp%%1!=0) stop("something wrong")
#                 tmp <- data.frame(
#                   method = rep(c("jointest","GRADE","SAODE"), nexp),
#                   N = rep(N[i],nexp*3), R = rep(R[i],nexp*3), SNR = rep(SNR[i],nexp*3)
#                 )
#                 aucs <- sapply(
#                   unique(ids),
#                   function(id) {
#                     estAUC(res$tpr[ids==id], res$fpr[ids==id])
#                   }
#                 )
#                 tmp$auc <- aucs
#                 return(tmp)
#               }) %>% do.call(rbind, .)

# tab2 <- dat2 %>%
#   mutate(method = factor(method, levels = c("jointest","GRADE","SAODE"))) %>% 
#   dplyr::filter(method %in% selected_methods) %>% 
#   group_by(method, N, R, SNR) %>% 
#   summarise(auc.mean = mean(auc),
#             auc.sd = sd(auc)) %>% 
#   arrange(desc(N), R, desc(SNR)) %>% 
#   mutate(N = N * R) %>% 
#   ungroup %>% 
#   dplyr::select(-R) %>% 
#   mutate(
#     auc.mean = round(auc.mean, digits = 3),
#     auc.sd = round(auc.sd, digits = 3)
#   )

N <- N * R

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
set_names <- paste0("N",rev(N),"SNR",rev(SNR))
expos <- c(-2, 0, 0, 0, 0, 0)
ndigits <- numfloats <- c(3,3,3,3,3,3)


# tab <- left_join(tab, tab2, by = c("method", "N", "SNR"))
tab <- tab %>% 
  select(all_of(c("N", "SNR", "method", var_names))) %>% 
  relocate(method, .after = SNR) %>% 
  mutate(
    method = case_when(
      method=="jointest" ~ "JADE",
      method=="GRADE" ~ "GRADE",
      method=="SAODE" ~ "SA-ODE"
    )
  )
for (i in seq_along(crit_names)) {
  for (ii in 1:2) {
    tab[,nsetcol+2*i-2+ii] <- sapply(tab[[nsetcol+2*i-2+ii]],
                                     format_num,
                                     expo = expos[i],
                                     digits = ndigits[i],
                                     num_float = numfloats[i])
  }
}


latex_simple <- c()
for (i in seq_len(nrow(tab))) {
  rowval <- as.character(tab[i,])
  rowlatex <- c(rowval[1:3],
                combine.mean_and_sd(rowval[4:length(rowval)],
                                    pt_reduce = 7))
  if (rowval[3]=="JADE") {
    rowlatex[1:2] <- paste0(
      r"{\multirow{3}{*}{ }",rowlatex[1:2],r"{ }}"
    )
  }
  if (rowval[3]=="GRADE") {
    rowlatex[1:2] <- " "
    rowlatex[4:5] <- paste0(
      r"{\multirow{2}{*}{ }",rowlatex[4:5],r"{ }}"
    )
  }
  if (rowval[3]=="SA-ODE") {
    rowlatex[1:2] <- " "
    rowlatex[4:5] <- " "
  }
  # rowlatex <- add.phantom(rowlatex, loc = c(4,1,1,1))
  rowlatex <- add.ampersand(rowlatex)  
  rowlatex <- paste(paste(rowlatex, collapse = " "), r"{\\}", "\n")
  latex_simple <- c(latex_simple, rowlatex)
  if (i %% 3 == 0) {
    latex_simple <- c(latex_simple, paste0(r"{\cmidrule{1-}",ncrit+nsetcol,r"{}}"," \n"))
    latex_simple <- c(latex_simple, "\n")
  }
}

cat(latex_simple)


## deprecated ========================

## format numerical values -----------------
crit_names <- c("theta.mse", "dtheta.mse", "ftheta.mse.t", "ftheta.mse.f")
ncrit <- length(crit_names)
var_names <- c(
  rbind(
    paste0(crit_names, ".mean"),
    paste0(crit_names, ".sd")
  )
)
set_names <- paste0("N",rev(N),"SNR",rev(SNR))
expos <- c(-1, -1, 0, -1)
ndigits <- numfloats <- c(3,3,3,3)


text_mat <- c()
for (meth in selected_methods) {
  # extract data
  meth.mat <- as.matrix(tab[tab$method==meth,var_names])
  rownames(meth.mat) <- set_names
  # round numbers
  meth.mat.round <- matrix(nrow = nrow(meth.mat), ncol = ncol(meth.mat))
  for (i in 1:ncrit) {
    for (ii in 1:2) {
      # meth.mat[,2*i-2+ii] <- div10Expo(meth.mat[,2*i-2+ii],
      #                                  expo = expos[i], digits = 2)
      meth.mat.round[,2*i-2+ii] <- sapply(
        meth.mat[,2*i-2+ii],
        format_num, expo = expos[i],
        digits = ndigits[i], num_float = numfloats[i]
      )
    }
  }
  meth.mat <- meth.mat.round
  # make text cells
  tmp <- matrix(t(meth.mat), nrow=2)  # 2 x (ncrit * nset )
  tmp <- makecell(tmp[1,], tmp[2,])  # (ncrit * nset)
  tmp <- t(matrix(tmp, nrow = ncrit))  #
  rownames(tmp) <- set_names
  text_mat <- rbind(text_mat, tmp)
}
text_mat <- matrix(text_mat, nrow = nset)
methods.crit.name <- c(t(sapply(selected_methods,
                            function(meth) paste(meth, crit_names, sep = "."))))
colnames(text_mat) <- methods.crit.name
rownames(text_mat) <- set_names


# merge the theta.mse and dtheta.mse of two-stage methods
if (nmethod==3) {
  colnames(text_mat)[nmethod-1] <- "spline.theta.mse"
  colnames(text_mat)[2*nmethod-1] <- "spline.dtheta.mse"
  text_mat <- text_mat[,-c(nmethod,2*nmethod)]
} else {
  stop("NOT IMPLEMENTED")
}


text_mat2 <- c()
for (meth in selected_methods) {
  # extract data
  meth.mat <- as.matrix(tab2[tab2$method==meth,c("auc.mean","auc.sd")])
  rownames(meth.mat) <- set_names
  # make text cells
  tmp <- matrix(t(meth.mat), nrow=2)  # 2 x (ncrit * nset )
  tmp <- makecell(tmp[1,], tmp[2,])  # (ncrit * nset)
  tmp <- t(matrix(tmp, nrow = 1))  #
  rownames(tmp) <- set_names
  text_mat2 <- rbind(text_mat2, tmp)
}
text_mat2 <- matrix(text_mat2, nrow = nset)
colnames(text_mat2) <-  c(t(sapply(selected_methods,
                                   function(meth) paste0(meth, ".auc"))))

text_mat <- cbind(text_mat, text_mat2)


## complete the latex table --------------

crit_width <- c(2,2,3,3,3)
setN <- c(
  r"{\multirow{2}{*}{40}}", " ",
  r"{\multirow{2}{*}{20}}", " "
)
setSNR <- as.character(c(10,4,10,4))
nsetcol <- 2



# header
latex_txt <- c(r"{\toprule}")
header <- c(
  r"{MSE $\btheta$}", r"{MSE $\rmd\btheta/\rmd t$}",
  r"{MSE $f_{jk}$ ($f_{jk}\quiv 0$)}", r"{MSE $f_{jk}$ ($f_{jk}\not\quiv 0$)}",
  r"{AUC}"
)
header[1:4] <- sapply(
  seq_along(header[1:4]),
  function(i) {
    if (expos[i]==0) return(header[i])
    return(makecell(header[i], paste0(r"{$\times 10^{}",expos[i],r"{}$}")))
  }
)
header <- paste0(
  r"{\multicolumn{}", crit_width,
  r"{}{c}{}", header, "}"
)
header <- header %>% 
  add.phantom() %>% 
  c()
  add.ampersand()
latex_txt <- c(latex_txt, paste(header, collapse = " "))
latex_txt <- c(latex_txt, )

# second header: method names
header2 <- c(
  rep(c("JADE", "smooth spline"), 2),
  rep(c("JADE", "GRADE", "SA-ODE"), 3)
)
header2 <- header2 %>% 
  add.phantom(loc = crit_width[-length(crit_width)]) %>% 
  add.ampersand()
latex_txt <- c(latex_txt, paste(header2, collapse = " "))

# header

