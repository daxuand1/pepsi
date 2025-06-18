
# Simulation/example code for the PSI estimator.

rm(list=ls(all=TRUE))
suppressPackageStartupMessages({
  library(survival)
  library(rootSolve)
  library(geepack)
  library(MASS)
  library(wgeesel)
  library(bindata)
  library(psych)
  library(tidyverse)
  library(foreach)
  library(doParallel)
  library(quadprog)
})
source("functions.R")

### main data parameter config
n=600
betaT = c(1, 1, 1, 1, 1)
sec_outcome_ind = 1:50 
n_sec = length(sec_outcome_ind)
thetaT = c( 1,  0,  0,  0,  1,  1,  1,  1,  0,  0,
            1,  0,  0,  0,  1,  1,  0,  0,  1,  1,
            0,  1,  1,  1,  0,  0,  0,  0,  1,  1,
            0,  1,  1,  1,  0,  0,  1,  1,  0,  0
)
thetaT = matrix(thetaT, nrow = 4, byrow = T)
thetaT = cbind(thetaT, matrix(1, nrow = 4, ncol = 40))
thetaT = as.matrix(thetaT[, sec_outcome_ind])
y_corr = matrix(0.5, 51, 51)
diag(y_corr) = 1
y_corr[1, 1 + 1:10] = y_corr[1 + 1:10, 1] = 
  c(0.8, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
y_corr = y_corr[c(0, sec_outcome_ind) + 1, 
                c(0, sec_outcome_ind) + 1]
x_corr = matrix(0, 4, 4)
for(i in 1:4){
  for(j in 1:4){
    x_corr[i, j] = 0.5^(abs(i-j))
  }
}
dist = "gaussian"


fn = c("pepsi_avg_simulation")
print(fn)

# parallel computation
plist = c("survival", "rootSolve", "geepack", 
          "MASS", "wgeesel",
          "bindata", "psych", "tidyverse", "quadprog")
n.cores <- 32

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)

iteration<-1e4
timestart<-Sys.time()
result = foreach(tt = 1:iteration, .combine = comb_fun,
                 .packages = plist) %dopar% 
{
  set.seed(tt)                 
  source("functions.R")
  rslt = list()                
  
  # generate main data
  x = mvrnorm(n, rep(0, 4), x_corr)
  
  mu = cbind(cbind(1, x) %*% betaT, x %*% thetaT)
  y = mu + mvrnorm(n, rep(0, n_sec+1), y_corr)
  y_main = y[, 1]
  y_sec = apply(as.matrix(y[, -1]), 2, function(x){x - mean(x)})
  
  # secondary model estimation begins here
  para_pool = rep()
  prop_score_pool = rep()
  
  threshold = 1e-3
  
  tau_list = c(0.001, seq(0.01, 0.5, by = 0.01),
               seq(0.6, 1, by = 0.1)) # grid search for the optimal tau
  
  for(i in 1:n_sec){
    
    sec_fit = glm(y_sec[, i] ~ x-1)
    para_init = sec_fit$coefficient
    
    para_list = rep()
    bic_list = rep()
    
    for(tau in tau_list){
      
      if(tau == min(tau_list) || 
         max(abs(para - para_init)) >= 2 * max(abs(para_init))){
        lambda = NULL
        para = para_init
      }
      
      para_fit = theta_find_newton(para, tau, x, y_sec[, i],
                                   lambda_init = lambda,
                                   threshold = threshold, maxit = 30)
      para = para_fit$para
      para_list = cbind(para_list, para)

      lambda=lambda_find(para, x, y_sec[, i])
      ef = elfun_glm(para, x, y_sec[, i])$elef
      bic = 2*(sum(R0der(lambda, ef))) + log(n)*sum(para != 0)
      bic_list = c(bic_list, bic)
      
      if(max(abs(para)) == 0 || bic > 2 * bic_list[1]){
        break
      }
      
    }
    ind = which.min(bic_list)
    
    para = para_list[, ind]
    para_pool = cbind(para_pool, para)
    ef = elfun_glm(para, x, y_sec[, i])$elef
    lambda = lambda_find(para, x, y_sec[, i])
    prop_score_pool = cbind(prop_score_pool, R1der(lambda, ef))
  }
  
  # secondary model estimation ends here
  
  # naive estimator
  beta_null = glm(y_main ~ x)$coefficient
  rslt$beta_null = unname(beta_null)
  
  # PEPSI estimator
  pepsi_result = pepsi(x, y_main, y_sec, para_pool, dist = dist)
  rslt$beta_pepsi = unname(pepsi_result$beta)
  rslt$Vnull = diag(pepsi_result$Vnull)
  rslt$Vpepsi = diag(pepsi_result$Vpepsi)
  
  # Avg estimator optimizing IIB
  avg_result = multi_avg(x, y_main, y_sec, 
                         para_pool, prop_score_pool, dist = dist)
  rslt$beta_avg = unname(avg_result$beta)
  rslt$Vavg = diag(avg_result$Vavg)
  
  return(rslt)
}
timeend<-Sys.time()
timeend-timestart

parallel::stopCluster(cl = my.cluster)

fdir = paste0("result/", fn, ".RData")
save.image(fdir)

beta_null = rowMeans(result$beta_null)
beta_pepsi = rowMeans(result$beta_pepsi)
beta_avg = rowMeans(result$beta_avg)

emp_var_null = diag(var(t(result$beta_null)))
emp_var_pepsi = diag(var(t(result$beta_pepsi)))
emp_var_avg = diag(var(t(result$beta_avg)))

est_var_null = rowMeans(result$Vnull)
est_var_pepsi = rowMeans(result$Vpepsi)
est_var_avg = rowMeans(result$Vavg)

cp_pepsi = sapply(1:iteration, function(i){
  return(result$beta_pepsi[, i] >= betaT - qnorm(0.975) *
           sqrt(result$Vpepsi[, i]) &
           result$beta_pepsi[, i] <= betaT + qnorm(0.975) *
           sqrt(result$Vpepsi[, i]))
})

cp_avg = sapply(1:iteration, function(i){
  return(result$beta_avg[, i] >= betaT - qnorm(0.975) *
           sqrt(result$Vavg[, i]) &
           result$beta_avg[, i] <= betaT + qnorm(0.975) *
           sqrt(result$Vavg[, i]))
})

print("Summary table:")

summarytable = data.frame(Bias_pepsi = (beta_pepsi-betaT)*100,
                          SD_pepsi = sqrt(emp_var_pepsi)*100,
                          SE_pepsi = sqrt(est_var_pepsi)*100,
                          CP_pepsi = rowMeans(cp_pepsi)*100,
                          RE_pepsi = emp_var_null/emp_var_pepsi,
                          Bias_avg = (beta_avg-betaT)*100,
                          SD_avg = sqrt(emp_var_avg)*100,
                          SE_avg = sqrt(est_var_avg)*100,
                          CP_avg = rowMeans(cp_avg)*100,
                          RE_avg = emp_var_null/emp_var_avg
)

round(summarytable[, 1:5], 2)
round(summarytable[, -c(1:5)], 2)

save.image(fdir)
