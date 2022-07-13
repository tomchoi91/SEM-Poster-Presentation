install.packages("lavaan")
?lavoptions

rm(list = ls())
library(MASS)
library(lavaan)
library(purrr) # for 'transpose'

############################ model implied cov #############################

phi<- matrix(c(1, rep(.837, 2), 1), ncol=2, byrow=TRUE)
phi

loadingx <- matrix(0, 4, 2)
loadingx[1:2, 1] <- c(.8, .8)
loadingx[2:4, 2] <- c(.335, .8, .8)
loadingx

delta <- matrix(0, 4, 4)
diag(delta) <- c(rep(.36, 4))
delta

sigma<- loadingx %*% phi %*% t(loadingx) + delta
sigma

rtsigma <- cov2cor(sigma)


############################ model specification #############################

# the true model ----------------------------------------------------------

model1 <- '
# measurement model
xi1  =~  d1*x1 + d1*x2 
xi2  =~  x2 + x3 + x4
# factor covariances
xi1 ~~ xi1
xi2 ~~ xi2
xi1 ~~ xi2'


# level 1 missspecification ---------------------------------------------

model2 <- '
# measurement model
xi1  =~  d1*x1 + d1*x2 
xi2  =~  x2 + x3 + x4
# factor covariances
xi1 ~~ xi1
xi2 ~~ xi2
xi1 ~~ 1*xi2'


# level 2 missspecification --------------------------------------------- # this one looks more severe

model3 <- '
# measurement model
xi1  =~  d1*x1 + d1*x2 
xi2  =~  0*x2 + x3 + x4
# factor covariances
xi1 ~~ xi1
xi2 ~~ xi2
xi1 ~~ 1*xi2'


############################ data generation #############################

gen_data2 <- function(model, n, ...){
  y <- mvrnorm(n, rep(0, 4), sigma)  # rep(mu = 0, n = 4)
  colnames(y) <- c("x1",  "x2",  "x3",  "x4")
  fit <- cfa(model, data = y, std.lv = TRUE)
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
        # invokeRestart("muffleWarning")
      }
    )
    if(is.na(warn)){
      list(
        result = result,
        warning = warn,
        error = error
      )
    }
    else{list(
      result = NULL,
      warning = warn,
      error = error
    )}
  }
}

robust_gen <- robustify(gen_data2)

# robustify <- function(fun) {
#   fun <- match.fun(fun)
#   
#   function(...) {
#     
#     error <- warn <- NA_character_
#     result <- withCallingHandlers(
#       tryCatch(
#         fun(...),
#         error = function(cnd) {
#           error <<- conditionMessage(cnd)
#           NULL
#         }
#       ),
#       warning = function(cnd) {
#         warn <<- conditionMessage(cnd)
#         invokeRestart("muffleWarning")
#       }
#     )
#     
#     list(
#       result = result,
#       warning = warn,
#       error = error
#     )
#     
#   }
# }


gen_data4 <- function(master.fit, w, mtype, msize, mlevel, nsize, rtime, ...) {
  # w2 <- as.data.frame(lapply(w[1:500,], mean)) #only using 100 because of nonconvergent solutions
  w2 <- w[1:500,]
  cell <- cbind(w2, mtype, msize, mlevel, nsize, rtime)
  master.fit <- rbind(master.fit, cell)
  
  return(master.fit)
}

# for gamma ---------------------------------------------------------------

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

cfa3.fit <- data.frame()
nrep <- 800 #take 800 reps

set.seed(12479)
w4 <- vector("list", nrep) # this is awesome!!!!!!!!!!!!!!!!!!!
for (i in 1:nrep) {
  w4[[i]] <- robust_gen(model3, 1e3)
}

w5 <- transpose(w4)
w6 <- w5$result[1:nrep]
# head(w6)

w7 <- dplyr::bind_rows(w6)
# w7 <- do.call(rbind, wtf)
# w7 <- data.table::rbindlist(wtf)
nrow(w7)
w7$gamma <- ncol(sigma)/(ncol(sigma) + 2*w7[1,]$df*((w7$rmsea)^2))


cfa3.fit <- gen_data4(cfa3.fit, w7, "CFA3", "small", "LV2", 1e3, "NA")
head(cfa3.fit)


# for test ----------------------------------------------------------------

set.seed(12479)
y <- mvrnorm(1e4, rep(0, 4), sigma)  # rep(mu = 0, n = 4)
colnames(y) <- c("x1",  "x2",  "x3",  "x4")

fit1 <- cfa(model1, data = y, std.lv= TRUE)
# fit2 <- cfa(model2, data = y, std.lv= TRUE)
fit3 <- cfa(model3, data = y, std.lv= TRUE)
# rm(list = "fit")
# fitmeasures(fit1, fit.measures = "rmsea")
# fitmeasures(fit2, fit.measures = "rmsea")
# fitmeasures(fit3, fit.measures = "rmsea")
# ?lavoptions
# lavInspect(fit, "std.lv")
summary(fit1)
summary(fit3)

h1 <- as.data.frame(as.list(fitmeasures(fit1)))
h1$gamma <- ncol(sigma)/(ncol(sigma) + 2*h1[1,]$df*((h1$rmsea)^2))
h2 <- as.data.frame(as.list(fitmeasures(fit2)))
h2$gamma <- ncol(sigma)/(ncol(sigma) + 2*h2[1,]$df*((h2$rmsea)^2))
h3 <- as.data.frame(as.list(fitmeasures(fit3)))
h3$gamma <- ncol(sigma)/(ncol(sigma) + 2*h3[1,]$df*((h3$rmsea)^2))
h1
h2
h3

cfa3.fit <- data.frame()
cfa3.fit <- gen_data5(cfa3.fit, h3, "CFA", "small", "lvl3", 1e6, "NA")
cfa3.fit

# exporting data -------------------------------------------------------

getwd()
setwd("C:/Users/dchoi2/Box/EDPS_971/Simulation studies/actual trials")
write.table(y, file = "cfa3_data.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(cfa3.fit, file = "cfa3_trial9.csv", row.names = FALSE, col.names = TRUE, sep = ",")
