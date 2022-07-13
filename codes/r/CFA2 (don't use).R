install.packages("lavaan")

rm(list = ls())
library(MASS)
library(lavaan)
library(purrr) # for 'transpose'

############################ model implied cov #############################

phi<- matrix(c(1, .854, .85, .854, 1, .859, .85, .859, 1), 3, byrow=TRUE)
phi

loadingx <- matrix(0, 15, 3)
loadingx[1:5, 1] <- c(.7, .7, .75, .8, .8)
loadingx[6:10, 2] <- c(.7, .7, .75, .8, .8)
loadingx[11:15, 3] <- c(.7, .7, .75, .8, .8)
loadingx

delta <- matrix(0, 15, 15)
diag(delta) <- c(.51, .51, .4375, .36, .36, .51, .51, .4375, .36, .36, .51, .51, .4375, .36, .36)
delta

sigma<- loadingx %*% phi %*% t(loadingx) + delta
sigma


############################ model specification #############################

# the true model ----------------------------------------------------------

model1 <- '
# measurement model
xi1  =~  NA*x1 + x2 + x3 + x4 + x5 
xi2  =~  NA*x6 + x7 + x8 + x9 + x10
xi3  =~  NA*x11 + x12 + x13 + x14 + x15
# factor covariances
xi1 ~~ 1*xi1
xi2 ~~ 1*xi2
xi3 ~~ 1*xi3
xi1 ~~ xi2
xi1 ~~ xi3
xi2 ~~ xi3'


# level 1 missspecification ---------------------------------------------

model2 <- '
# measurement model
xi1  =~  NA*x1 + x2 + x3 + x4 + x5 
xi2  =~  NA*x6 + x7 + x8 + x9 + x10
xi3  =~  NA*x11 + x12 + x13 + x14 + x15
# factor covariances
xi1 ~~ 1*xi1
xi2 ~~ 1*xi2
xi3 ~~ 1*xi3
xi1 ~~ 0*xi2
xi1 ~~ xi3
xi2 ~~ xi3'


# level 2 missspecification ---------------------------------------------

model3 <- '
# measurement model
xi1  =~  NA*x1 + x2 + x3 + x4 + x5 
xi2  =~  NA*x6 + x7 + x8 + x9 + x10
xi3  =~  NA*x11 + x12 + x13 + x14 + x15
# factor covariances
xi1 ~~ 1*xi1
xi2 ~~ 1*xi2
xi3 ~~ 1*xi3
xi1 ~~ 0*xi2
xi1 ~~ 0*xi3
xi2 ~~ xi3'


############################ data generation #############################

gen_data2 <- function(model, n, ...){
  y <- mvrnorm(n, rep(0, 15), sigma)  # rep(mu = 0, n = 15)
  colnames(y) <- c(paste0("x", 1:15))
  fit <- sem(model, data = y)
  # h <- list()
  # h <- append(h, fitmeasures(fit, fit.measures = "all"))
  # data.frame(h)
  h <- as.data.frame(as.list(fitmeasures(fit)))
  
  return(h)
}


robustify <- function(fun) {
  fun <- match.fun(fun)
  
  function(...) {
    
    error <- warn <- NA_character_
    result <- withCallingHandlers(
      tryCatch(
        fun(...),
        error = function(cnd) {
          error <<- conditionMessage(cnd)
          NULL
        }
      ),
      warning = function(cnd) {
        warn <<- conditionMessage(cnd)
        invokeRestart("muffleWarning")
      }
    )
    
    list(
      result = result,
      warning = warn,
      error = error
    )
    
  }
}

robust_gen <- robustify(gen_data2)

gen_data4 <- function(master.fit, w, mtype, msize, mlevel, nsize, rtime, ...) {
  # w2 <- dplyr::bind_rows(w)
  # w2 <- do.call(rbind, w)
  # w2 <- data.table::rbindlist(w)
  w2 <- as.data.frame(lapply(w[1:500,], mean)) #only using 100 because of nonconvergent solutions
  cell <- cbind(w2, mtype, msize, mlevel, nsize, rtime)
  master.fit <- rbind(master.fit, cell)
  
  return(master.fit)
}

gen_data5 <- function(master.fit, w, mtype, msize, mlevel, nsize, rtime, ...) {
  cell <- cbind(w, mtype, msize, mlevel, nsize, rtime)
  master.fit <- rbind(master.fit, cell)
  
  return(master.fit)
}


# deprecated  --------------------------------------------------------

# install.packages("furrr")
# install.packages("purrr")
# library(furrr)
# library(purrr)
# library(tidyr)
# ?safely
# 
# safe_gen_data2 <- purrr::safely(quietly(gen_data2))
# 
# set.seed(12479)
# 
# plan(multiprocess, workers = 4L)
# 
# w4 <- future_map(1:nrep, ~ safe_gen_data2(model2, 1e2)) %>%
#   transpose()
# 
# # wtf <- w4$result %>%
# #   simplify_all()

############################ fit indices gathering #############################

cfa2.fit <- data.frame()
nrep <- 800 #take 300 reps

set.seed(12479)
w4 <- vector("list", nrep) # this is awesome!!!!!!!!!!!!!!!!!!!
for (i in 1:nrep) {
  w4[[i]] <- robust_gen(model3, 1e2)
}

w5 <- transpose(w4)
w6 <- w5$result[1:nrep]
head(w6)

w7 <- dplyr::bind_rows(w6)
# w7 <- do.call(rbind, wtf)
# w7 <- data.table::rbindlist(wtf)
nrow(w7)
w7$gamma <- ncol(sigma)/(ncol(sigma) + 2*w7[1,]$df*((w7$rmsea)^2))


cfa2.fit <- gen_data4(cfa2.fit, w7, "CFA2", "big", "lvl3", 1e2, "7m29s")
cfa2.fit


# for test ----------------------------------------------------------------

set.seed(12479)
y <- mvrnorm(1e2, rep(0, 15), sigma)  # rep(mu = 0, n = 15)
colnames(y) <- c(paste0("x", 1:15))

fit1 <- cfa(model1, data = y)
fit2 <- cfa(model2, data = y)
fit3 <- cfa(model3, data = y)
# rm(list = "fit")
fitmeasures(fit1, fit.measures = "rmsea")
fitmeasures(fit2, fit.measures = "rmsea")
fitmeasures(fit3, fit.measures = "rmsea")
# ?lavoptions
# lavInspect(fit, "std.lv")
summary(fit)

h1 <- as.data.frame(as.list(fitmeasures(fit1)))
h1$gamma <- ncol(sigma)/(ncol(sigma) + 2*h1[1,]$df*((h1$rmsea)^2))
h2 <- as.data.frame(as.list(fitmeasures(fit2)))
h2$gamma <- ncol(sigma)/(ncol(sigma) + 2*h2[1,]$df*((h2$rmsea)^2))
h3 <- as.data.frame(as.list(fitmeasures(fit3)))
h3$gamma <- ncol(sigma)/(ncol(sigma) + 2*h3[1,]$df*((h3$rmsea)^2))
h1
h2
h3

cfa2.fit <- data.frame()
cfa2.fit <- gen_data5(cfa2.fit, h3, "CFA", "big", "lvl3", 1e6, "NA")
cfa2.fit

# exporting data -------------------------------------------------------

getwd()
setwd("C:/Users/dchoi2/Box/EDPS_971/Simulation studies/actual trials")
write.table(y, file = "cfa2_data.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(cfa2.fit, file = "cfa2_trial1.csv", row.names = FALSE, col.names = TRUE, sep = ",")
