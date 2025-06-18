###basic functions to solve empirical likelihood by Owen(2001) and Han(2019)


## log EL
R0der=function(lambda, ZZ)
{
  eps = 1/ncol(ZZ)
  lambda_h = 1 + as.vector(t(ZZ) %*% lambda)
  r0der = lambda_h
  ind = lambda_h >= eps
  r0der[ind] = log(lambda_h[ind])
  r0der[!ind] = log(eps)-1.5+2*lambda_h[!ind]/eps-lambda_h[!ind]^2/eps^2/2
  return(r0der)
}

##first derivative of log EL
R1der=function(lambda,ZZ)
{
  eps = 1/ncol(ZZ)
  lambda_h = 1 + as.vector(t(ZZ) %*% lambda)
  r1der = ifelse(lambda_h >= eps, 1/lambda_h,
                 2/eps-lambda_h/eps^2)
  return(r1der)
}

##second derivative of log EL
R2der=function(lambda, ZZ)
{
  eps = 1/ncol(ZZ)
  lambda_h = 1 + as.vector(t(ZZ) %*% lambda)
  r2der = ifelse(lambda_h >= eps, -1/lambda_h^2, -1/eps^2)
  return(r2der)
}

#function to find lambda, given theta, by Han(2019)
lambda_find=function(para, x, y, ZZ = NULL)
{
  if(is.null(ZZ)){
    ZZ=elfun_glm(para, x, y)$elef
  }
  
  lambda_old = rep(0,nrow(ZZ))
  k=0
  tol=1e-6
  
  repeat{
    # step 1
    rl = colSums(R1der(lambda_old,ZZ) * t(ZZ))
    rll = ZZ %*% diag(R2der(lambda_old,ZZ)) %*% t(ZZ) 
    gamma=0
    
    # step 2
    repeat{
      update = as.vector(2^(-gamma) * ginv(rll)%*%rl)
      lambda = as.vector(lambda_old - update)
      
      
      if(sum(is.na(update)) > 0 |
         sum(update == Inf) > 0 |
         sum(update == -Inf) > 0){
        lambda = lambda_old
        break
      }else if(max(abs(update))<tol){
        break
      }
      
      # if(gamma>10){
      #   break
      # }
      
      index_1 = as.vector(1+t(ZZ) %*% lambda) <= 1/ncol(ZZ)
      index_2 = sum(R0der(lambda,ZZ))<=sum(R0der(lambda_old,ZZ))
      
      if (sum(index_1)>0 | index_2>0){
        gamma=gamma+1
      }else{
        break
      }
    }
    
    # step 3
    if(max(abs(lambda-lambda_old))<tol){
      break
    }else{
      lambda_old = lambda
    }
    
  }
  
  return(lambda)
}

# --- functions for PSI, PEPSI and AVG ---

elfun_glm=function(para, x, y, dist = "gaussian")
{
  if(!is.matrix(para)){para = as.matrix(para, nrow = length(para))}
  if(!is.matrix(y)){y = as.matrix(y, nrow = length(y))}
  
  n = nrow(x)
  p = ncol(x)
  
  elef = rep()
  delef = list()
  
  for(i in 1:(dim(y)[2]) ){
    mu = as.vector(x %*% para[, i])
    if(dist == "binomial"){
      mu = exp(mu) / (1 + exp(mu))
    }
    
    elef = rbind(elef, t(x * as.vector(y[, i] - mu)) )
    
    if(dist == "gaussian"){
      delef_i = sapply(1:n, function(j){
        - x[j,] %*% t(x[j,])
      }, simplify = "array")
    }else{
      delef_i = sapply(1:n, function(j){
        - x[j,] %*% t(x[j,]) * mu[j] * (1 - mu[j])
      }, simplify = "array")
    }
    dim(delef_i) = c(p, p*n)
    
    delef[[i]] = delef_i
  }
  
  if((dim(y)[2]) == 1){
    delef = delef[[1]]
  }

  return(list(elef = elef, delef = delef))
}

scad = function(para, t = tau, a = 3.7){
  ifelse(abs(para)<=t, t*abs(para),
         ifelse(abs(para)<=a*t, (-para^2+2*a*t*abs(para)-t^2)/2/(a-1),
                (a+1)*t^2/2))
}

scad_der = function(para, t = tau, a = 3.7){
  t * (as.numeric(abs(para) <= t) + 
           pmax(a * t - abs(para), 0) * 
           as.numeric(abs(para) > t) / (a - 1) / t)
}

theta_find_newton = function(para_init, tau, x, y, 
                                lambda_init = NULL,
                                threshold = 1e-3, maxit = 30,
                                bound = 30)
{
  #tau = 0.05;y=y_sec[,i];maxit = 30;bound=30;
  para_old = para_init
  if(is.null(lambda_init)){
    lambda = lambda_find(para_old, x, y)
  }else{
    lambda = lambda_init
  }
  k = 0
  n = nrow(x)
  
  repeat{
    
    total=elfun_glm(para_old, x, y)
    ZZ=total$elef
    ZZ_d=total$delef
    m_old = sum(R0der(lambda, ZZ)) + sum(scad(para_old, tau)) * n
    
    n = ncol(ZZ)
    p = nrow(ZZ)    
    non_zero_ind = para_old != 0 # only update non-zero coefficients
    ZZ_d[, rep(!non_zero_ind, n)] = 0
    
    scaler = R1der(lambda, ZZ)
    penalty_factor = ifelse(para_old == 0, 0,
                            scad_der(para_old, tau)/abs(para_old))
    
    # calculate gradient and hessian following Han(2019)
    # h - a pxp matrix
    h = ZZ_d
    dim(h) = c(p^2, n)
    h = t(t(h) * scaler)
    dim(h) = c(p, p, n)
    h1 = apply(h, 1:2, sum)
    m_der = t(h1) %*% lambda + penalty_factor * para_old * n
    
    m_der2 = - t(h1) %*% 
      ginv(ZZ %*% diag(R2der(lambda,ZZ)) %*% t(ZZ)) %*% (h1) + 
      diag(penalty_factor) * n
    
    if(sum(is.na(m_der2)) > 0|
       sum(m_der2 == Inf) > 0|
       sum(m_der2 == -Inf) > 0){
      direction = rep(0, length(m_der))
    }else{
      direction = ginv(m_der2) %*% m_der
    }
    
    for(sigma in 0:5){
      step_length = 2^(-sigma)
      
      para_temp = as.vector(para_old - step_length * direction)
      ZZ_temp=elfun_glm(para_temp, x, y)$elef
      #lambda_temp = lambda_find(para_temp, elfun)
      m_temp = sum(R0der(lambda, ZZ_temp)) +
        sum(scad(para_old, tau) +
              penalty_factor*(para_temp^2 - para_old^2)/2) * n
      if(m_temp <= m_old){
        break
      }
    }
    
    para = para_temp
    para = ifelse(abs(para) < threshold, 0, para)
    
    if(max(abs(para - para_old)) < threshold | k > maxit | 
       sum(para != 0) == 0 | max(abs(para)) > bound){
      break
    }else{
      para_old = para
      k = k+1
      lambda=lambda_find(para_old, x, y)
    }
    
  }
  
  #para = ifelse(abs(para) < 1e-3, 0, para)
  return(list(para = para,
              converge = max(abs(para - para_old)) < threshold |
                sum(para != 0) == 0
                ))
  #return(para = para)
}

weighted_glm = function(beta_null, x, y, weight, dist){
  
  k = 0
  beta_old = beta_null
  p = length(beta_null)
  n = nrow(x)
  
  repeat{
    total = elfun_glm(beta_old, x, y, dist)
    sf = colSums(t(total$elef) * weight)
    dsf = total$delef
    dim(dsf) = c(p*p, n)
    dsf = t(t(dsf)*weight)
    dim(dsf) = c(p, p, n)
    dsf = apply(dsf, 1:2, sum)
    
    beta = as.vector(beta_old - ginv(dsf) %*% sf)
    
    if(max(abs(beta - beta_old)) < 1e-6 | k > 20){
      break
    }else{
      beta = beta_old
      k = k+1
    }
  }
  
  return(beta)
}

psi = function(x, y_main, y_sec, para, prop_score,
                dist = "gaussian", x_sec = x){
  
  sec_total = elfun_glm(para, x_sec, y_sec, dist = dist)
  ZZ = sec_total$elef
  ZZ_d = sec_total$delef
  
  n = ncol(ZZ)
  r = nrow(ZZ)
  
  if(sum(para == 0) == 0){
    beta = glm(y_main ~ x, family = dist)$coefficients
    main_total = elfun_glm(beta, cbind(1, x), y_main, dist = dist)
    sf = main_total$elef
    dsf = main_total$delef
    t = nrow(sf)
    
    Sigma = sf %*% t(sf) / n
    Gamma = dsf
    dim(Gamma) = c(t, t, n)
    Gamma = apply(Gamma, 1:2, mean)
    
    Vnull = ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))
    
    return(list(beta = beta,
                Vnull = Vnull / n,
                Vpsi = Vnull / n,
                Lambda = sf %*% t(ZZ) / n,
                S = matrix(0, r, r),
                ZZ = ZZ))
    
  }
  
  H = diag(r)
  H = matrix(H[para == 0, ], nrow = sum(para == 0))
  
  A = ZZ %*% t(ZZ) / n
  B = ZZ_d
  dim(B) = c(r, r, n)
  B = apply(B, 1:2, mean)
  C = ginv(t(B) %*% ginv(A) %*% B)
  
  S = ginv(A) - ginv(A) %*% B %*% C %*% t(B) %*% ginv(A)
  P = ginv(A) %*% B %*% C %*% t(H) %*% ginv(H %*% C %*% t(H)) %*%
    H %*% C %*% t(B) %*% ginv(A)
  
  if(dist == "gaussian"){
    beta = glm(y_main ~ x, weights = prop_score)$coefficient
  }else{
    beta_null = glm(y_main ~ x, family = dist)$coefficient
    beta = weighted_glm(beta_null, cbind(1, x), y_main, prop_score, dist)
  }
  
  main_total = elfun_glm(beta, cbind(1, x), y_main, dist = dist)
  sf = main_total$elef
  t = nrow(sf)
  dsf = main_total$delef
  
  Lambda = sf %*% t(ZZ) / n
  Sigma = sf %*% t(sf) / n
  Gamma = dsf
  dim(Gamma) = c(t, t, n)
  Gamma = apply(Gamma, 1:2, mean)
  
  Vnull = ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))
  Vpsi = ginv(Gamma) %*% 
    (Sigma - Lambda %*% (S + P) %*% t(Lambda)) %*% t(ginv(Gamma))
  
  return(list(beta = beta,
              Vnull = Vnull / n,
              Vpsi = Vpsi / n,
              Lambda = Lambda,
              SP = S+P,
              ZZ = ZZ))
}

pepsi = function(x, y_main, y_sec, para, 
                 dist = "gaussian", x_sec = x)
{
  if(!is.matrix(para)){para = as.matrix(para, nrow = length(para))}
  
  select_ind = apply(para, 2, function(x){sum(x == 0) > 0})
  n = nrow(x)
  
  if(sum(select_ind) == 0){
    beta = glm(y_main ~ x, family = dist)$coefficient
    main_total = elfun_glm(beta, cbind(1, x), y_main, dist = dist)
    sf = main_total$elef
    dsf = main_total$delef
    t = nrow(sf)
    
    Sigma = sf %*% t(sf) / n
    Gamma = dsf
    dim(Gamma) = c(t, t, n)
    Gamma = apply(Gamma, 1:2, mean)
    
    Vnull = ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))
    
    return(list(beta = beta,
                Vnull = Vnull / n,
                Vpepsi = Vnull / n))
    
  }
  
  
  y_sec_select = as.matrix(y_sec[, select_ind])
  para_select = as.matrix(para[, select_ind])
  n_sec_select = sum(select_ind)
  
  # pca on ef
  
  sec_total = elfun_glm(para_select, x_sec, y_sec_select)
  r = nrow(para_select)
  
  R = rep()
  
  for( i in 1:n_sec_select ){
    
    ZZ = sec_total$elef[((i-1)*r+1):(i*r), ]
    
    H = diag(r)
    H = matrix(H[para_select[, i] == 0, ], nrow = sum(para_select[, i] == 0))
    
    A = ZZ %*% t(ZZ) / n
    if(n_sec_select == 1){
      B = sec_total$delef
    }else{
      B = sec_total$delef[[i]]
    }
    dim(B) = c(r, r, n)
    B = apply(B, 1:2, mean)
    C = ginv(t(B) %*% ginv(A) %*% B)
    
    S = ginv(A) - ginv(A) %*% B %*% C %*% t(B) %*% ginv(A)
    P = ginv(A) %*% B %*% C %*% t(H) %*% ginv(H %*% C %*% t(H)) %*%
      H %*% C %*% t(B) %*% ginv(A)
    
    R_i = diag(r) - B %*% C %*% t(B) %*% ginv(A) +
      B %*% C %*% t(H) %*% ginv(H %*% C %*% t(H)) %*%
      H %*% C %*% t(B) %*% ginv(A)
    R = cbind(R, rbind(matrix(0, (i-1)*r, r),
                       R_i,
                       matrix(0, (n_sec_select-i)*r, r)))
    
  }
  
  pca_fit = summary(prcomp(t(R %*% sec_total$elef)))
  # num_pc = which.max(pca_fit$importance[3, ] >= 0.9)
  # num_pc = length(as.vector(para))
  num_pc = sum(para_select == 0)
  # rslt$num_pc = num_pc
  
  
  ZZ = t(as.matrix(pca_fit$rotation[, 1:num_pc])) %*% R %*% sec_total$elef
  # lambda = lambda_find(para_select, x_sec, y_sec_select, ZZ)
  # prop_score = R1der(lambda, ZZ)
  if(num_pc == 1){
    lambda = (pca_fit$sdev[1])^(-2) * rowMeans(ZZ)
  }else{
    lambda = diag((pca_fit$sdev[1:num_pc])^(-2)) %*% rowMeans(ZZ)
  }
  
  prop_score = 1 - as.vector(t(ZZ) %*% lambda)
  

  beta_null = glm(y_main ~ x, family = dist)$coefficient
  beta = weighted_glm(beta_null, cbind(1, x), y_main, prop_score, dist)
  
  main_total = elfun_glm(beta, cbind(1, x), y_main, dist = dist)
  sf = main_total$elef
  t = nrow(sf)
  dsf = main_total$delef
  
  Lambda = sf %*% t(ZZ) / n
  Sigma = sf %*% t(sf) / n
  Gamma = dsf
  dim(Gamma) = c(t, t, n)
  Gamma = apply(Gamma, 1:2, mean)
  
  Vnull = ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))
  if(num_pc == 1){
    Vpepsi = ginv(Gamma) %*% 
      ( Sigma - (pca_fit$sdev[1])^(-2) * Lambda %*% t(Lambda) ) %*% 
      t(ginv(Gamma))
  }else{
    Vpepsi = ginv(Gamma) %*% 
      ( Sigma - 
          Lambda %*% diag((pca_fit$sdev[1:num_pc])^(-2)) %*% t(Lambda) ) %*% 
      t(ginv(Gamma))
  }
  
  return(list(beta = beta,
              Vnull = Vnull / n,
              Vpepsi = Vpepsi / n))
}

multi_avg = function(x, y_main, y_sec, para, prop_score_pool, 
                     dist = "gaussian", prior = c(0, rep(1, ncol(x))),
                     x_sec = x)
{
  # x = x_bs; y_main = y_main_bs; y_sec = y_sec_bs;
  # para = para_pool; prior = c(0, rep(1, ncol(x)));
  # x_sec = x 
  beta_null = glm(y_main ~ x, family = dist)$coefficient
  
  select_ind = apply(para, 2, function(x){sum(x == 0) > 0})
  n = nrow(x)
  t = length(beta_null)
  n_sec = ncol(y_sec)
  
  if(sum(select_ind) == 0){
    main_total = elfun_glm(beta_null, cbind(1, x), y_main, dist = dist)
    sf = main_total$elef
    dsf = main_total$delef
    
    Sigma = sf %*% t(sf) / n
    Gamma = dsf
    dim(Gamma) = c(t, t, n)
    Gamma = apply(Gamma, 1:2, mean)
    
    Vnull = ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))
    
    return(list(beta = beta_null,
                Vnull = Vnull / n,
                Vavg = Vnull / n))
    
  }
  
  if(sum(select_ind) == 1){
    psi_result = psi(x, y_main, y_sec[, select_ind], 
                     para[, select_ind], prop_score_pool[, select_ind])
    
    return(list(beta = unname(psi_result$beta),
                Vnull = psi_result$Vnul / n,
                Vavg = psi_result$Vpsi / n))
  }
  
  y_sec_select = y_sec[, select_ind]
  para_select = para[, select_ind] 
  prop_score_pool_select = prop_score_pool[, select_ind]
  n_sec_select = sum(select_ind)
  
  main_total = elfun_glm(beta_null, cbind(1, x), y_main, dist = dist)
  sf = main_total$elef
  Sigma = sf %*% t(sf) / n
  
  sec_total = lapply(1:n_sec_select, function(i){
    psi(x, y_main, y_sec_select[, i], 
        para_select[, i], prop_score_pool_select[, i],
        dist = dist, x_sec = x_sec)})
  
  # apply quadprog
  Dmat = matrix(0, n_sec_select, n_sec_select)
  dvec = rep(0, n_sec_select)
  Amat = t(rbind(rep(1, n_sec_select), diag(n_sec_select)))
  bvec = c(1, rep(0, n_sec_select))
  
  for(i in 1:n_sec_select){
    total_i = sec_total[[i]]
    
    for(j in i:n_sec_select){
      if(i == j){
        Dmat[i, i] = dvec[i] = 
          sum(diag(total_i$Lambda %*% 
                     total_i$SP %*% t(total_i$Lambda)) /
                diag(Sigma) * prior)
      }else{
        total_j = sec_total[[j]]
        Dmat[i, j] = Dmat[j, i] = 
          sum(diag(total_i$Lambda %*% total_i$SP %*% 
                     (total_i$ZZ %*% t(total_j$ZZ) / n) %*%
                     t(total_j$SP) %*% t(total_j$Lambda)) /
                diag(Sigma) * prior)
      }
    }
  }
  
  select_ind2 = diag(Dmat) > 0
  Dmat = Dmat[select_ind2, select_ind2]
  dvec = dvec[select_ind2]
  Amat = matrix(1, sum(select_ind2), 1)
  bvec = c(1)
  
  weight = rep(0, n_sec)
  w_subset = rep(0, sum(select_ind))
  w_subset[select_ind2] = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
  weight[select_ind] = w_subset
  weight = round(weight, 4)
  
  prop_score = as.vector(prop_score_pool %*% weight)
  
  beta = weighted_glm(beta_null, cbind(1, x), y_main, prop_score, dist)
  
  main_total = elfun_glm(beta, cbind(1, x), y_main, dist = dist)
  sf = main_total$elef
  dsf = main_total$delef
  
  Sigma = sf %*% t(sf) / n
  Gamma = dsf
  dim(Gamma) = c(t, t, n)
  Gamma = apply(Gamma, 1:2, mean)
  
  Vnull = ginv(Gamma) %*% Sigma %*% t(ginv(Gamma))
  
  Vavg = Sigma
  for(i in 1:n_sec_select){
    total_i = sec_total[[i]]
    
    for(j in i:n_sec_select){
      if(i == j){
        Vavg = Vavg + (weight[i]^2 - 2*weight[i]) *
          total_i$Lambda %*% total_i$SP %*% t(total_i$Lambda)
      }else{
        total_j = sec_total[[j]]
        cov_ij = total_i$Lambda %*% total_i$SP %*% 
          (total_i$ZZ %*% t(total_j$ZZ) / n) %*%
          t(total_j$SP) %*% t(total_j$Lambda)
        Vavg = Vavg + weight[i] * weight[j] * (cov_ij + t(cov_ij))
      }
    }
  }
  
  Vavg = ginv(Gamma) %*% Vavg %*% t(ginv(Gamma))
  
  return(list(beta = beta,
              Vnull = Vnull / n,
              Vavg = Vavg / n,
              w = weight))
}

# --- other functions ---

comb_fun = function(list1, list2){
  mapply(cbind, list1, list2, SIMPLIFY=F)
}

