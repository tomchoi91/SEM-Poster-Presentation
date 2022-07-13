install.packages("lavaan")

rm(list = ls())
library(MASS)
library(lavaan)
library(purrr) # for 'transpose'

############################ model implied cov #############################

loadingx <- matrix(c(1, .5))
loadingy <- matrix(c(1, .75, 0, 0, 0, 0, 1, 0.9), ncol=2)
gamma <- matrix(c(-.6, -.435))
beta <- matrix(c(rep(0,2), .6, 0), ncol=2, byrow=TRUE)
phi <- matrix(7, 1)
psi <- matrix(c(5, rep(0, 2), 4), ncol =2, byrow=TRUE)
delta <- matrix(0, 2, 2)
diag(delta) <- c(3, 2.5)
epsil <- matrix(c(0, 0, 1.16, rep(0, 4), 1.1, 1.16, rep(0, 4), 1.1, 0, 0), ncol= 4)
diag(epsil) <- c(4.75, 2.5, 4.5, 3)
ID <- matrix(0,nrow(beta),nrow(beta))
diag(ID) <- 1
a<- loadingy %*% solve(ID-beta) %*% (gamma %*% phi %*% t(gamma) + psi) %*% solve(ID-t(beta)) %*% t(loadingy) + epsil
b<- loadingy %*% solve(ID-beta) %*% gamma %*% phi %*% t(loadingx)
c<- loadingx %*%  phi %*% t(gamma) %*% solve(ID-t(beta)) %*% t(loadingy)
d<- loadingx %*% phi %*% t(loadingx) + delta
e<- cbind(a,b)
f<- cbind(c,d)
sigma<- rbind(e,f)
sigma

# e<- cbind(b,a)
# f<- cbind(d,c)
# sigma<- rbind(f,e)
# sigma
# cov2cor(sigma)

# loadingx <- matrix(c(1, .5))
# loadingy <- matrix(c(1, .95, 0, 0, 0, 0, 1, 0.9), ncol=2)
# gamma <- matrix(c(-.6, -.25))
# beta <- matrix(c(rep(0,2), .6, 0), ncol=2, byrow=TRUE)
# phi <- matrix(7, 1)
# psi <- matrix(c(5, rep(0, 2), 4), ncol =2, byrow=TRUE)
# delta <- matrix(0, 2, 2)
# diag(delta) <- c(3, 2.5)
# epsil <- matrix(c(0, 0, 1.60, rep(0, 4), 0.3, 1.60, rep(0, 4), 0.3, 0, 0), ncol= 4)
# diag(epsil) <- c(4.75, 2.5, 4.5, 3)
# ID <- matrix(0,nrow(beta),nrow(beta))
# diag(ID) <- 1
# a<- loadingy %*% solve(ID-beta) %*% (gamma %*% phi %*% t(gamma) + psi) %*% solve(ID-t(beta)) %*% t(loadingy) + epsil
# b<- loadingy %*% solve(ID-beta) %*% gamma %*% phi %*% t(loadingx)
# c<- loadingx %*%  phi %*% t(gamma) %*% solve(ID-t(beta)) %*% t(loadingy)
# d<- loadingx %*% phi %*% t(loadingx) + delta
# e<- cbind(a,b)
# f<- cbind(c,d)
# sigma<- rbind(e,f)
# sigma

# for sas code ------------------------------------------------------------

# sigma<- matrix(c(10, 3.5, -4.2, -3.99, -4.27, -3.843,
#                 3.5, 4.25, -2.1, -1.995, -2.135, -1.9215,
#                 -4.2, -2.1, 12.27, 7.144, 7.162, 5.0058,
#                 -3.99, -1.995, 7.144, 9.2868, 5.2839, 5.05551,
#                 -4.27, -2.135, 7.162, 5.2839, 12.9047, 7.56423,
#                 -3.843, -1.9215, 5.0058, 5.05551, 7.56423, 9.807807), ncol = 6)
# sigma

# sigma<- matrix(c(1, .536875, -.379164, -.414038, -.375884, -.388047,
#                  .536875, 1, -.290805, -317552, -.288290, -.297619,
#                  -.379164, -.290805, 1, .669246, .569164, .456315,
#                  -.414038, -317552, .669246, 1, .482667, .529719,
#                  -.375884, -.288290, .569164, .482667, 1, .672364,
#                  -.388047, -.297619, .456315, .529719, .672364, 1), ncol = 6)
# sigma

############################ model specification #############################

# the true model ----------------------------------------------------------

model1 <- '
# measurement model
xi1   =~  1*x1 + x2 
eta1  =~  1*y1 + y2
eta2  =~  1*y3 + y4 
# regresssions
eta1 ~ xi1
eta2 ~ xi1 + eta1 
# estimator covariances
y1 ~~ y3
y2 ~~ y4'


# level 1 missspecification ---------------------------------------------

model2 <- '
# measurement model
xi1   =~  1*x1 + x2 
eta1  =~  1*y1 + y2
eta2  =~  1*y3 + y4 
# regresssions
eta1 ~ xi1
eta2 ~ 0*xi1 + eta1 
# estimator covariances
y1 ~~ y3
y2 ~~ y4'


# level 2 missspecification ---------------------------------------------

model3 <- '
# measurement model
xi1   =~  1*x1 + d1*x2 
eta1  =~  1*y1 + d1*y2
eta2  =~  1*y3 + y4 
# regresssions
eta1 ~ xi1
eta2 ~ 0*xi1 + eta1 
# estimator covariances
y1 ~~ 0*y3
y2 ~~ 0*y4'


############################ data generation #############################

gen_data2 <- function(model, n, ...){
  y <- mvrnorm(n, rep(0, 6), sigma)  # rep(mu = 0, n = 6)
  colnames(y) <- c("y1", "y2", "y3", "y4", "x1", "x2")
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
# invokeRestart("mufflewarning")
# ?withCallingHandlers
# ?invokeRestart

robust_gen <- robustify(gen_data2)

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

sem2.fit <- data.frame()
nrep <- 800 #take 800 reps

set.seed(12479)
w4 <- vector("list", nrep) # this is awesome!!!!!!!!!!!!!!!!!!!
for (i in 1:nrep) {
  w4[[i]] <- robust_gen(model3, 1e3)
}

# lavInspect(w4[[277]], "std.lv")
w5 <- transpose(w4)
w6 <- w5$result[1:nrep]
# head(w6)

w7 <- dplyr::bind_rows(w6)
# w7 <- do.call(rbind, wtf)
# w7 <- data.table::rbindlist(wtf)
nrow(w7)
w7$gamma <- ncol(sigma)/(ncol(sigma) + 2*w7[1,]$df*((w7$rmsea)^2))
# w8 <- w7
# w9 <- rbind(w8, w7)
# nrow(w9)
# head(w8)
# head(w9)

sem2.fit <- gen_data4(sem2.fit, w7, "SEM2", "small", "LV2", 1e3, "NA")
head(sem2.fit)


# for test ----------------------------------------------------------------

set.seed(12479)
y <- mvrnorm(1e6, rep(0, 6), sigma)  # rep(mu = 0, n = 6)
# colnames(y) <- c("y1", "y2", "y3", "y4", "x1", "x2")
colnames(y) <- c(paste0("y", 1:4),paste0("x", 1:2))


fit1 <- sem(model1, data = y)
# fit2 <- sem(model2, data = y)
# lavInspect(fit3, "std.lv")
# lavInspect(fit3, "theta")
# fit3 <- sem(model3, data = y)

summary(fit1, fit.measures = TRUE)
?`summary,lavaan-method`
# resid(fit3, type = "raw")
?`fitMeasures,lavaan-method`
?`fitted,lavaan-method`
?`fitted.values,lavaan-method`
?fitmeasures

fitmeasures(fit3)
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

sem2.fit <- data.frame()
sem2.fit <- gen_data5(sem2.fit, h3, "SEM", "small", "lvl3", 1e6, "NA")
sem2.fit

# exporting data -------------------------------------------------------

getwd()
setwd("C:/Users/dchoi2/Box/EDPS_971/Simulation studies/actual trials")
write.table(y, file = "sem2_data.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sem2.fit, file = "sem2_trial9.csv", row.names = FALSE, col.names = TRUE, sep = ",")
