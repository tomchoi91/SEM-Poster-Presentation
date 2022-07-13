install.packages("lavaan")


rm(list = ls())
library(MASS)
library(lavaan)
library(purrr) # for 'transpose'

############################ model implied cov #############################

loadingx <- matrix(c(.9, .8, .8, 0, 0, 1, 0, .45, .8, 1.2, 1, 0), ncol=2)
loadingy <- matrix(c(.9, .7, .7, 0, 0, 1, 0, .45, .7, .9, 1, 0), ncol=2)
gamma <- matrix(c(1.1, 0, 0, 1), ncol=2, byrow=TRUE)
beta <- matrix(c(rep(0,4)), ncol=2, byrow=TRUE)
phi <- matrix(c(105, rep(90, 2),  115), ncol =2, byrow=TRUE)
psi <- matrix(c(25, rep(16.2, 2), 25), ncol =2, byrow =TRUE)
delta <- matrix(0,6,6)
diag(delta) <- c(30, 30, 40, 45, 20, 50)
epsil <- matrix(0,6,6)
diag(epsil) <- c(25, 40, 50, 40, 20, 70)
ID <- matrix(0,nrow(beta),nrow(beta))
diag(ID) <- 1

a<- loadingy %*% solve(ID-beta) %*% (gamma %*% phi %*% t(gamma) + psi) %*% solve(ID-t(beta)) %*% t(loadingy) + epsil
b<- loadingy %*% solve(ID-beta) %*% gamma %*% phi %*% t(loadingx)
c<- loadingx %*%  phi %*% t(gamma) %*% solve(ID-t(beta)) %*% t(loadingy)
d<- loadingx %*% phi %*% t(loadingx) + delta
e<- cbind(a,b)
f<- cbind(c,d)
sigma <- rbind(e,f)
sigma
cor<- cov2cor(sigma)

# e<- cbind(b,a)
# f<- cbind(d,c)
# sigma<- rbind(f,e)
# sigma
# cor<- cov2cor(sigma)



############################ model specification #############################

# the true model ----------------------------------------------------------

model1 <- '
# measurement model
xi1   =~  NA*x1 + x2 + x3 + 1*x6 
xi2   =~  NA*x2 + x3 + x4 + 1*x5
eta1  =~  NA*y1 + y2 + y3 + 1*y6
eta2  =~  NA*y2 + y3 + y4 + 1*y5 
# regresssions
eta1 ~ xi1
eta2 ~ xi2
# factor covariances
xi1 ~~ xi2
eta1 ~~ eta2'

# level 1 missspecification ---------------------------------------------

model2 <- '
# measurement model
xi1   =~  NA*x1 + x2 + x3 + 1*x6 
xi2   =~  0*x2 + x3 + x4 + 1*x5
eta1  =~  NA*y1 + y2 + y3 + 1*y6
eta2  =~  0*y2 + y3 + y4 + 1*y5 
# regresssions
eta1 ~ xi1
eta2 ~ xi2
# factor covariances
xi1 ~~ xi2
eta1 ~~ eta2'


# level 2 missspecification ---------------------------------------------

model3 <- '
# measurement model
xi1   =~  NA*x1 + x2 + x3 + 1*x6 
xi2   =~  0*x2 + x3 + x4 + 1*x5
eta1  =~  NA*y1 + y2 + y3 + 1*y6
eta2  =~  0*y2 + y3 + y4 + 1*y5 
# regresssions
eta1 ~ xi1
eta2 ~ xi2
# factor covariances
xi1 ~~ xi2
eta1 ~~ 0*eta2'


############################ data generation #############################

gen_data2 <- function(model, n, ...){
  y <- mvrnorm(n, rep(0, 12), sigma)  # rep(mu = 0, n = 12)
  colnames(y) <- c(paste0("y", 1:6), paste0("x", 1:6))
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
#           error <<- conditionMessage(cnd)     # when error occurs the error message gets returned
#           NULL                                # because no result is generated  
#         }                                     # so we need to put NULL here
#       ),
#       warning = function(cnd) {               # but when warning occurs generation function proceeds and  
#         warn <<- conditionMessage(cnd)        # result is generated so no warning message gets returned
#         invokeRestart("muffleWarning")
#       }
#     )
#     
#     list(
#       result = result,                        # this is where the error message is pasted by = operator
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

sem1.fit <- data.frame()
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

sem1.fit <- gen_data4(sem1.fit, w7, "SEM1", "big", "LV2", 1e3, "NA")
head(sem1.fit)
warnings()

# deprecated --------------------------------------------------------------
# 
# set.seed(12479)
# w <- vector("list", nrep) # this is awesome!!!!!!!!!!!!!!!!!!!
# for (i in 1:nrep) {
#   w[[i]] <- gen_data2(model3, n = 1e2)
# }
# w


# obtaining the "population parameters" ----------------------------------------------------------------

set.seed(12479)
y <- mvrnorm(1e6, rep(0, 12), sigma)  # rep(mu = 0, n = 12)
colnames(y) <- c(paste0("y", 1:6), paste0("x", 1:6))

fit1 <- sem(model1, data = y)
fit2 <- sem(model2, data = y)
fit3 <- sem(model3, data = y)
summary(fit1)

fitmeasures(fit1, fit.measures = "rmsea")
fitmeasures(fit2, fit.measures = "rmsea")
fitmeasures(fit3, fit.measures = "rmsea")

h1 <- as.data.frame(as.list(fitmeasures(fit1)))
h1$gamma <- ncol(sigma)/(ncol(sigma) + 2*h1[1,]$df*((h1$rmsea)^2))
h2 <- as.data.frame(as.list(fitmeasures(fit2)))
h2$gamma <- ncol(sigma)/(ncol(sigma) + 2*h2[1,]$df*((h2$rmsea)^2))
h3 <- as.data.frame(as.list(fitmeasures(fit3)))
h3$gamma <- ncol(sigma)/(ncol(sigma) + 2*h3[1,]$df*((h3$rmsea)^2))
h1
h2
h3

sem1.fit <- data.frame()
sem1.fit <- gen_data5(sem1.fit, h3, "SEM", "big", "lvl3", 1e6, "NA")
sem1.fit

# exporting data -------------------------------------------------------

getwd()
setwd("C:/Users/dchoi2/Box/EDPS_971/Simulation studies/actual trials")
write.table(y, file = "sem1_data.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sem1.fit, file = "sem1_trial9.csv", row.names = FALSE, col.names = TRUE, sep = ",")
