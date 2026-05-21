##########################################################################################################################
#functions for the conditional dependence testing via                                                                    #
#adapted from the code from "A Two-Sample Conditional Distribution Test Using Conformal Prediction and Weighted Rank Sum"#
##########################################################################################################################

getUstat <- function(g1, v1, v2){
  ## to calculate \sum_j \hat{U}_j/n21
  n1 <- length(g1); n2 <- length(v2)
  zeta <- runif(n2)
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1
  inner_sum <- sapply(v1, function(x){mean(zeta*(x==v2))})
  u1 <- mean(g1*((1-ecdf1)+inner_sum))
  return(u1/mean(g1))
}

## calculate the empirical cdf under 1st sample
getEcdf1 <- function(v1, v2){
  fun <- ecdf(v2)
  ecdf1 <- fun(v1)
  ecdf2 <- sapply(v1, function(x){return(mean(x==v2))})
  return(list(ecdf1=ecdf1, ecdf2=ecdf2))
}

## calculate the empirical cdf under 2nd sample
getEcdf2 <- function(v){
  fun <- ecdf(v)
  ecdf1 <- fun(v)
  ecdf2 <- sapply(v, function(x){return(mean(x==v))})
  return(list(ecdf1=ecdf1, ecdf2=ecdf2))
}

## calculate the asymptotic variance under sample 2 (cancellation one g)
## ecdf1: P(V_2<=v), ecdf2: P(V_2=v)
getAsyVar2 <- function(g2, v2, m1){
  m2 <- length(g2)
  ecdfVal <- getEcdf2(v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  s1 <- mean(g2*(1-ecdf1)^2)-1/4
  s2 <- mean(g2*ecdf2^2)/3
  s3 <- mean(g2*(1-ecdf1)*ecdf2)
  ss <- s1+s2+s3
  ss <- ifelse(ss>0, ss, 0)
  sigma1 <- ss/m1 + 1/12/m2
  sigma2 <- (mean(g2)-1)/m1
  sigma2 <- ifelse(sigma2>0, sigma2, 0)
  sigma12 <- (mean(g2*(1-ecdf1))+mean(g2*ecdf2)/2 - 1/2)/m1
  sigma <- sigma1 + sigma2/4 - sigma12 
  return(list(sigma=sigma, sigma1=sigma1, sigma2=sigma2, sigma12=sigma12))
}


## calculate the asymptotic variance under 1st sample
getAsyVar1 <- function(g1, v1, v2){
  m1 <- length(g1); m2 <- length(v2)
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  ss <- mean(g1^2*(1-ecdf1)^2)+mean(g1^2*ecdf2^2)/3+mean(g1^2*(1-ecdf1)*ecdf2)-
    c(1/4, mean(g1*(1-ecdf1))+1/2*mean(g1*ecdf2)-1/4, (mean(g1*(1-ecdf1))+mean(g1*ecdf2)/2)^2)   
  ss <- ifelse(ss>0, ss, 0)
  sigma1 <- ss/m1 + 1/12/m2
  sigma2 <- (mean(g1^2)-c(1, 2*mean(g1)-1, mean(g1)^2))/m1
  sigma2 <- ifelse(sigma2>0, sigma2, 0)
  sigma12 <- (mean(g1^2*(1-ecdf1))+mean(g1^2*ecdf2)/2 - 
                c(1/2, mean(g1)/2+mean(g1*(1-ecdf1)+g1*ecdf2/2)-1/2, mean(g1)*mean(g1*(1-ecdf1)+g1*ecdf2/2)))/m1
  sigma <- sigma1 + sigma2/4 - sigma12 
  return(list(sigma=sigma, sigma1=sigma1, sigma2=sigma2, sigma12=sigma12))
}


## calculate the final statistic
getFinalStat <- function(g1, g2, v1, v2){
  m1 <- length(g1); m2 <- length(g2)
  U <- getUstat(g1, v1, v2)
  sigma.1 <- getAsyVar1(g1, v1, v2)$sigma # length: 3
  sigma.2 <- getAsyVar2(g2, v2, m1)$sigma
  sigma.gm <- sqrt(sigma.1*sigma.2)
  sigma.hm <- 2/(1/sigma.1+1/sigma.2)
  z1 <- (U-0.5)/sqrt(sigma.1)
  z2 <- (U-0.5)/sqrt(sigma.2)
  z.gm <- (U-0.5)/sqrt(sigma.gm)
  z.hm <- (U-0.5)/sqrt(sigma.hm)
  return(list(U=U, z1=z1, z2=z2, z.gm=z.gm, z.hm=z.hm, sigma.1=sigma.1, sigma.2=sigma.2, 
              sigma.gm=sigma.gm, sigma.hm=sigma.hm))
}


NNfun <- function(x, z, xpre1, xpre2, nnrep=10, hidden.layers=NA, acfun = 'sigmoid', optim.type='sgd',
                  n.epochs=500, learn.rates=0.001, L1=0){
  n1 <- dim(xpre1)[1]; n2 <- dim(xpre2)[1]
  prob1 <- matrix(0, nrow = n1, ncol = nnrep)
  prob2 <- matrix(0, nrow = n2, ncol = nnrep)
  for(i in 1:nnrep){
    print(paste0('Neural Net:', i))
    fit_nn <- neuralnetwork(x, z, hidden.layers = hidden.layers, optim.type = optim.type,
                            val.prop=0, learn.rates = learn.rates, L1=L1, 
                            n.epochs = n.epochs, activ.functions = acfun)
    prob1[,i] <- predict(fit_nn, xpre1)$probabilities[,2]
    prob2[,i] <- predict(fit_nn, xpre2)$probabilities[,2]
  }
  prob1.fit <- rowMeans(prob1)
  prob2.fit <- rowMeans(prob2)
  prob1.fit[prob1.fit<0.01] <- 0.01; prob1.fit[prob1.fit>0.99] <- 0.99
  prob2.fit[prob2.fit<0.01] <- 0.01; prob2.fit[prob2.fit>0.99] <- 0.99
  return(list(prob1.fit=prob1.fit, prob2.fit=prob2.fit))
}


KLR.CV <- function(x, y, xpre, ypre, lambdaseq, sigmaseq){
  lambda_len <- length(lambdaseq); sigma_len <- length(sigmaseq)
  centropy <- matrix(0, lambda_len, sigma_len)
  ypre <- as.numeric(as.character(ypre))
  minerror <- Inf
  for(i in 1:lambda_len){
    lambda <- lambdaseq[i]
    for(j in 1:sigma_len){
      print(paste0('KLR CV: lambda', i, 'sigma', j))
      sigma <- sigmaseq[j]
      cv_data <- constructData(x, y)
      cv_data <- shuffleData(cv_data)
      klr_learner <- constructKlogRegLearner()
      # rbf: exp(-sigma\|x-y\|^2)
      params <- list(kernel='rbfdot', sigma=sigma, lambda=lambda/getN(cv_data), tol=10e-6, maxiter=500)
      fit_klr <- klr_learner$learn(cv_data, params)
      K = kernelMult(fit_klr$kernel, xpre, fit_klr$data, fit_klr$alpha) # predict
      prob_pre = 1 / (1 + exp(-as.vector(K))) # predicted probabilities
      # cross entropy error
      prob_pre[prob_pre<0.01] <- 0.01; prob_pre[prob_pre>0.99] <- 0.99
      centropy[i, j] <- mean(-ypre*log(prob_pre)-(1-ypre)*log(1-prob_pre))
      if(centropy[i,j]<minerror){
        minerror <- centropy[i,j]
        iind <- i; jind <- j
        prob <- prob_pre
      }
    }
  }
  prob[prob<0.01] <- 0.01; prob[prob>0.99] <- 0.99
  lambda <- lambdaseq[iind]; sigma <- sigmaseq[jind]
  return(list(prob=prob, lambda=lambda, sigma=sigma, centropy=centropy))
}

## calculate the cross entropy error
## y:label, 0 or 1; p: estimated probability for label 1
centropy <- function(p, y){
  return(mean(-y*log(p)-(1-y)*log(1-p)))
}

gerror <- function(g1, g2){
  sum(abs(g1/sum(g1)-g2/sum(g2)))
}

## calculate relevant quantities in assumptions
getCor <- function(g1,v1,v2,gorac1,vorac1,vorac2){
  mean_hatg_minus_g <- mean(g1-gorac1)
  mean_abs_hatg_minus_g <- mean(abs(g1-gorac1))
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  mean_hatd <- mean(1-ecdf1 + ecdf2/2)
  var_hatg_minus_g <- var(g1-gorac1)
  second_hatd <- mean(1-ecdf1+ecdf2/3)
  var_hatd <- second_hatd - mean_hatd^2
  second_hatg_minus_g_hatd <- mean((g1-gorac1)*((1-ecdf1)+ecdf2/2))
  e1 <- second_hatg_minus_g_hatd - mean_hatg_minus_g*mean_hatd # covariance between \hat G- G and \hat D
  e2 <- mean_hatg_minus_g
  e3 <- mean(gorac1*((1-ecdf1)+ecdf2/2)) - mean_hatd # covariance between G and \hat D
  # correlation between \hat G- G and \hat D
  rho1 <- e1 / sqrt(var_hatg_minus_g*var_hatd)
  rho2 <- mean_hatg_minus_g/mean_abs_hatg_minus_g
  # correlation between G and \hat D
  var_g <- var(gorac1)
  rho3 <- e3/sqrt(var_hatd*var_g)
  ecdfValorac <- getEcdf1(vorac1, vorac2)
  ecdforac1 <- ecdfValorac$ecdf1; ecdforac2 <- ecdfValorac$ecdf2
  e4 <- mean(gorac1*((1-ecdf1)+ecdf2/2)) - mean(gorac1*((1-ecdforac1)+ecdforac2/2)) # mean of G(\hat D - D)
  return(c('rho1'=rho1, 'rho2'=rho2, 'rho3'=rho3, 'e1'=e1, 'e2'=e2, 'e3'=e3, 'e4'=e4))
}


LR <- function(x1, x2, y1, y2){
  # H0:a1=a2; H1:a1!=a2 (assume b1=b2, sigma2.0=sigma2.1 under either H0 or H1)
  n1 <- dim(x1)[1]; n2 <- dim(x2)[1]; p <- dim(x1)[2]; n <- n1 + n2
  ## restricted
  x <- rbind(x1, x2); y <- c(y1,y2)
  xmean <- matrix(colMeans(x), n, p, byrow = T)
  ymean <- mean(y)
  b <- solve(crossprod(x-xmean), crossprod(x-xmean, y-ymean))
  a <- ymean - sum(b*colMeans(x))
  sigma2_res <- mean((y - a - x%*%b)^2)
  ## unrestricted
  xmean1 <- matrix(colMeans(x1), n1, p, byrow = T)
  xmean2 <- matrix(colMeans(x2), n2, p, byrow = T)
  ymean1 <- mean(y1); ymean2 <- mean(y2)
  b <- solve(crossprod(x1-xmean1)+crossprod(x2-xmean2), crossprod(x1-xmean1, y1-ymean1)+
               crossprod(x2-xmean2, y2-ymean2))
  a1 <- ymean1 - sum(b*colMeans(x1)); a2 <- ymean2 - sum(b*colMeans(x2))
  sigma2_unres <- (sum((y1-a1-x1%*%b)^2) + sum((y2-a2-x2%*%b)^2)) / (n1+n2)
  stat <- n * log(sigma2_res/sigma2_unres)
  pval <- 1 - pchisq(stat, 1)
  return(pval)
}

generate_data <- function(n, mu, Sigma){
  x <- mvrnorm(2*n, mu=mu, Sigma = Sigma)
  s <- ncol(x)
  ind <- which(g.orac(x)<100&g.orac(x)>0.01)
  x <- x[ind,]
  n_res <- n - length(ind)
  while(n_res>0){
    xx <- mvrnorm(2*n_res, mu=mu, Sigma = Sigma)
    ind <- which(g.orac(xx)<100&g.orac(xx)>0.01)
    x <- rbind(x, xx[ind,])
    n_res <- n_res - length(ind)
  }
  return(x[1:n,])
}

g.orac <- function(x){
  x = matrix(x, ncol=p)
  exp(x[,1:4] %*% c(1,1,-1,-1))/exp(2)
}

sigma_corr <- function(rho, s){
  sigma <- matrix(0, s, s)
  for(i in 1:s){
    for(j in 1:s){
      sigma[i, j] <- rho^(abs(i-j))
    }
  }
  return(sigma)
}

mars_model_df_normal <- function(n, p_nuisance, sigma, mean_seq, var_seq){
  X1 <- rnorm(n, mean_seq[1], var_seq[1])
  X2 <- rnorm(n, mean_seq[2], var_seq[2])
  X3 <- rnorm(n, mean_seq[3], var_seq[3])
  X4 <- rnorm(n, mean_seq[4], var_seq[4])
  X5 <- rnorm(n, mean_seq[5], var_seq[5])
  #Y1 <- 0.1 * exp(4 * X1) +
  #  4/(1+exp(-20 * (X2 - 0.5))) + 3 * X3 +
  #  2 * X4 + X5
  Y1 <- 50 * sin(pi * X1 * X2) + 5 * (X3 - 0.05)^2 + 5 * X4 + 5 * X5
  random_noise <- rnorm(n, 0, sigma)
  Y <- Y1 + random_noise
  if(p_nuisance > 0){
    #nuisance parameters with standard normal random noise:
    X_nuiss <- matrix(NA, n, p_nuisance)
    for(i in 1:p_nuisance){
      X_nuiss[,i] <- rnorm(n, 0, 1)
    }
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, X_nuiss, Y)
    colnames(data_mat) <-
      c("Y1", paste("X", c(1:(5+as.numeric(p_nuisance))), sep = ""), "Y")
  }
  else{
    data_mat <- cbind(Y1, X1, X2, X3, X4, X5, Y)
    colnames(data_mat) <-
      c("Y1", paste("X", c(1:5), sep = ""), "Y")
  }
  return(data.frame(data_mat))
}


## functions 


getUstat <- function(g1, v1, v2){
  ## to calculate \sum_j \hat{U}_j/n21
  
  n1 <- length(g1); n2 <- length(v2)
  zeta <- runif(n2)
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1
  inner_sum <- sapply(v1, function(x){mean(zeta*(x==v2))})
  u1 <- mean(g1*((1-ecdf1)+inner_sum))
  return(u1/mean(g1))
}


## calculate the empirical cdf under 1st sample
getEcdf1 <- function(v1, v2){
  fun <- ecdf(v2)
  ecdf1 <- fun(v1)
  ecdf2 <- sapply(v1, function(x){return(mean(x==v2))})
  return(list(ecdf1=ecdf1, ecdf2=ecdf2))
}


## calculate the empirical cdf under 2nd sample
getEcdf2 <- function(v){
  fun <- ecdf(v)
  ecdf1 <- fun(v)
  ecdf2 <- sapply(v, function(x){return(mean(x==v))})
  return(list(ecdf1=ecdf1, ecdf2=ecdf2))
}


## calculate the asymptotic variance under sample 2 (cancellation one g)
## ecdf1: P(V_2<=v), ecdf2: P(V_2=v)
getAsyVar2 <- function(g2, v2, m1){
  m2 <- length(g2)
  ecdfVal <- getEcdf2(v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  s1 <- mean(g2*(1-ecdf1)^2)-1/4
  s2 <- mean(g2*ecdf2^2)/3
  s3 <- mean(g2*(1-ecdf1)*ecdf2)
  ss <- s1+s2+s3
  ss <- ifelse(ss>0, ss, 0)
  sigma1 <- ss/m1 + 1/12/m2
  sigma2 <- (mean(g2)-1)/m1
  sigma2 <- ifelse(sigma2>0, sigma2, 0)
  sigma12 <- (mean(g2*(1-ecdf1))+mean(g2*ecdf2)/2 - 1/2)/m1
  sigma <- sigma1 + sigma2/4 - sigma12 
  return(list(sigma=sigma, sigma1=sigma1, sigma2=sigma2, sigma12=sigma12))
}




## calculate the asymptotic variance under 1st sample
getAsyVar1 <- function(g1, v1, v2){
  m1 <- length(g1); m2 <- length(v2)
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  ss <- mean(g1^2*(1-ecdf1)^2)+mean(g1^2*ecdf2^2)/3+mean(g1^2*(1-ecdf1)*ecdf2)-
    c(1/4, mean(g1*(1-ecdf1))+1/2*mean(g1*ecdf2)-1/4, (mean(g1*(1-ecdf1))+mean(g1*ecdf2)/2)^2)   
  ss <- ifelse(ss>0, ss, 0)
  sigma1 <- ss/m1 + 1/12/m2
  sigma2 <- (mean(g1^2)-c(1, 2*mean(g1)-1, mean(g1)^2))/m1
  sigma2 <- ifelse(sigma2>0, sigma2, 0)
  sigma12 <- (mean(g1^2*(1-ecdf1))+mean(g1^2*ecdf2)/2 - 
                c(1/2, mean(g1)/2+mean(g1*(1-ecdf1)+g1*ecdf2/2)-1/2, mean(g1)*mean(g1*(1-ecdf1)+g1*ecdf2/2)))/m1
  sigma <- sigma1 + sigma2/4 - sigma12 
  return(list(sigma=sigma, sigma1=sigma1, sigma2=sigma2, sigma12=sigma12))
}


## calculate the final statistic
getFinalStat <- function(g1, g2, v1, v2){
  m1 <- length(g1); m2 <- length(g2)
  U <- getUstat(g1, v1, v2)
  sigma.1 <- getAsyVar1(g1, v1, v2)$sigma # length: 3
  sigma.2 <- getAsyVar2(g2, v2, m1)$sigma
  sigma.gm <- sqrt(sigma.1*sigma.2)
  sigma.hm <- 2/(1/sigma.1+1/sigma.2)
  z1 <- (U-0.5)/sqrt(sigma.1)
  z2 <- (U-0.5)/sqrt(sigma.2)
  z.gm <- (U-0.5)/sqrt(sigma.gm)
  z.hm <- (U-0.5)/sqrt(sigma.hm)
  return(list(U=U, z1=z1, z2=z2, z.gm=z.gm, z.hm=z.hm, sigma.1=sigma.1, sigma.2=sigma.2, 
              sigma.gm=sigma.gm, sigma.hm=sigma.hm))
}


NNfun <- function(x, z, xpre1, xpre2, nnrep=10, hidden.layers=NA, acfun = 'sigmoid', optim.type='sgd',
                  n.epochs=500, learn.rates=0.001, L1=0){
  
  n1 <- dim(xpre1)[1]; n2 <- dim(xpre2)[1]
  prob1 <- matrix(0, nrow = n1, ncol = nnrep)
  prob2 <- matrix(0, nrow = n2, ncol = nnrep)
  
  for(i in 1:nnrep){
    print(paste0('Neural Net:', i))
    fit_nn <- neuralnetwork(x, z, hidden.layers = hidden.layers, optim.type = optim.type,
                            val.prop=0, learn.rates = learn.rates, L1=L1, 
                            n.epochs = n.epochs, activ.functions = acfun)
    
    prob1[,i] <- predict(fit_nn, xpre1)$probabilities[,2]
    prob2[,i] <- predict(fit_nn, xpre2)$probabilities[,2]
  }
  
  prob1.fit <- rowMeans(prob1)
  prob2.fit <- rowMeans(prob2)
  prob1.fit[prob1.fit<0.01] <- 0.01; prob1.fit[prob1.fit>0.99] <- 0.99
  prob2.fit[prob2.fit<0.01] <- 0.01; prob2.fit[prob2.fit>0.99] <- 0.99
  
  return(list(prob1.fit=prob1.fit, prob2.fit=prob2.fit))
}




KLR.CV <- function(x, y, xpre, ypre, lambdaseq, sigmaseq){
  lambda_len <- length(lambdaseq); sigma_len <- length(sigmaseq)
  centropy <- matrix(0, lambda_len, sigma_len)
  ypre <- as.numeric(as.character(ypre))
  minerror <- Inf
  for(i in 1:lambda_len){
    lambda <- lambdaseq[i]
    for(j in 1:sigma_len){
      print(paste0('KLR CV: lambda', i, 'sigma', j))
      sigma <- sigmaseq[j]
      cv_data <- constructData(x, y)
      cv_data <- shuffleData(cv_data)
      klr_learner <- constructKlogRegLearner()
      # rbf: exp(-sigma\|x-y\|^2)
      params <- list(kernel='rbfdot', sigma=sigma, lambda=lambda/getN(cv_data), tol=10e-6, maxiter=500)
      fit_klr <- klr_learner$learn(cv_data, params)
      K = kernelMult(fit_klr$kernel, xpre, fit_klr$data, fit_klr$alpha) # predict
      prob_pre = 1 / (1 + exp(-as.vector(K))) # predicted probabilities
      # cross entropy error
      prob_pre[prob_pre<0.01] <- 0.01; prob_pre[prob_pre>0.99] <- 0.99
      centropy[i, j] <- mean(-ypre*log(prob_pre)-(1-ypre)*log(1-prob_pre))
      if(centropy[i,j]<minerror){
        minerror <- centropy[i,j]
        iind <- i; jind <- j
        prob <- prob_pre
      }
    }
  }
  prob[prob<0.01] <- 0.01; prob[prob>0.99] <- 0.99
  lambda <- lambdaseq[iind]; sigma <- sigmaseq[jind]
  return(list(prob=prob, lambda=lambda, sigma=sigma, centropy=centropy))
  
}


## calculate the cross entropy error
## y:label, 0 or 1; p: estimated probability for label 1
centropy <- function(p, y){
  return(mean(-y*log(p)-(1-y)*log(1-p)))
}


gerror <- function(g1, g2){
  sum(abs(g1/sum(g1)-g2/sum(g2)))
}


## calculate relevant quantities in assumptions
getCor <- function(g1,v1,v2,gorac1,vorac1,vorac2){
  mean_hatg_minus_g <- mean(g1-gorac1)
  mean_abs_hatg_minus_g <- mean(abs(g1-gorac1))
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  mean_hatd <- mean(1-ecdf1 + ecdf2/2)
  var_hatg_minus_g <- var(g1-gorac1)
  second_hatd <- mean(1-ecdf1+ecdf2/3)
  var_hatd <- second_hatd - mean_hatd^2
  second_hatg_minus_g_hatd <- mean((g1-gorac1)*((1-ecdf1)+ecdf2/2))
  
  e1 <- second_hatg_minus_g_hatd - mean_hatg_minus_g*mean_hatd # covariance between \hat G- G and \hat D
  e2 <- mean_hatg_minus_g
  e3 <- mean(gorac1*((1-ecdf1)+ecdf2/2)) - mean_hatd # covariance between G and \hat D
  
  # correlation between \hat G- G and \hat D
  rho1 <- e1 / sqrt(var_hatg_minus_g*var_hatd)
  rho2 <- mean_hatg_minus_g/mean_abs_hatg_minus_g
  # correlation between G and \hat D
  var_g <- var(gorac1)
  rho3 <- e3/sqrt(var_hatd*var_g)
  
  ecdfValorac <- getEcdf1(vorac1, vorac2)
  ecdforac1 <- ecdfValorac$ecdf1; ecdforac2 <- ecdfValorac$ecdf2
  e4 <- mean(gorac1*((1-ecdf1)+ecdf2/2)) - mean(gorac1*((1-ecdforac1)+ecdforac2/2)) # mean of G(\hat D - D)
  
  return(c('rho1'=rho1, 'rho2'=rho2, 'rho3'=rho3, 'e1'=e1, 'e2'=e2, 'e3'=e3, 'e4'=e4))
}




LR <- function(x1, x2, y1, y2){
  # H0:a1=a2; H1:a1!=a2 (assume b1=b2, sigma2.0=sigma2.1 under either H0 or H1)
  n1 <- dim(x1)[1]; n2 <- dim(x2)[1]; p <- dim(x1)[2]; n <- n1 + n2
  ## restricted
  x <- rbind(x1, x2); y <- c(y1,y2)
  xmean <- matrix(colMeans(x), n, p, byrow = T)
  ymean <- mean(y)
  b <- solve(crossprod(x-xmean), crossprod(x-xmean, y-ymean))
  a <- ymean - sum(b*colMeans(x))
  sigma2_res <- mean((y - a - x%*%b)^2)
  ## unrestricted
  xmean1 <- matrix(colMeans(x1), n1, p, byrow = T)
  xmean2 <- matrix(colMeans(x2), n2, p, byrow = T)
  ymean1 <- mean(y1); ymean2 <- mean(y2)
  b <- solve(crossprod(x1-xmean1)+crossprod(x2-xmean2), crossprod(x1-xmean1, y1-ymean1)+
               crossprod(x2-xmean2, y2-ymean2))
  a1 <- ymean1 - sum(b*colMeans(x1)); a2 <- ymean2 - sum(b*colMeans(x2))
  sigma2_unres <- (sum((y1-a1-x1%*%b)^2) + sum((y2-a2-x2%*%b)^2)) / (n1+n2)
  stat <- n * log(sigma2_res/sigma2_unres)
  pval <- 1 - pchisq(stat, 1)
  return(pval)
}


sim_fun <- function(x1, y1, x2, y2, n11, n12, n21, n22){
  n1 <- n11+n12
  n2 <- n21+n22
  index1 <- sample(1:n1, size = n11)
  x11 <- x1[index1,]; x12 <- x1[-index1,]
  y11 <- y1[index1]; y12 <- y1[-index1]
  index2 <- sample(1:n2, size = n21)
  x21 <- x2[index2,]; x22 <- x2[-index2,]
  y21 <- y2[index2]; y22 <- y2[-index2]
  g12.orac <- rep(1, n12)
  g22.orac <- rep(1, n22)
  v12.orac <- rep(1, n12)
  v22.orac <- rep(1, n22)
  temp <- getFinalStat(g12.orac, g22.orac, v12.orac, v22.orac)
  u.orac <- temp$U
  var.orac1 <- temp$sigma.1
  oracpvalue1 <- pnorm(temp$z1)
  var.orac2 <- temp$sigma.2
  oracpvalue2 <- pnorm(temp$z2)
  oracpvalue.gm <- pnorm(temp$z.gm)
  oracpvalue.hm <- pnorm(temp$z.hm)
  label.fit <- as.factor(c(rep(0,n11), rep(1,n21)))
  xy.fit <- data.frame(x.fit=rbind(x11, x21), y.fit=c(y11, y21))
  fit.joint <- glm(label.fit~., data=xy.fit, family="binomial")
  x.fit <- data.frame(x.fit=rbind(x11,x21))
  fit.marginal <- glm(label.fit~., data=x.fit, family="binomial")
  prob.marginal <- predict(fit.marginal, newdata=data.frame(x.fit=rbind(x12,x22)), type="response")
  prob.marginal[prob.marginal<0.01] <- 0.01; prob.marginal[prob.marginal>0.99] <- 0.99
  g12.est.ll <- prob.marginal[1:n12]/(1-prob.marginal[1:n12])*n11/n21
  g22.est.ll <- prob.marginal[(n12+1):(n12+n22)]/(1-prob.marginal[(n12+1):(n12+n22)])*n11/n21
  cerror12.ll.marginal <- mean(prob.marginal[1:n12]>0.5)
  cerror22.ll.marginal <- mean(prob.marginal[(n12+1):(n12+n22)]<0.5)
  centropy.ll.marginal <- centropy(prob.marginal, c(rep(0, n12), rep(1, n22)))
  error.ll <- gerror(g12.est.ll, g12.orac)
  prob.joint <- predict(fit.joint, newdata=data.frame(x.fit=rbind(x12,x22), y.fit=c(y12,y22)), type="response")
  prob.joint[prob.joint<0.01] <- 0.01; prob.joint[prob.joint>0.99] <- 0.99
  v12.est.ll <- (1-prob.joint[1:n12])/prob.joint[1:n12]*g12.est.ll
  v22.est.ll <- (1-prob.joint[(n12+1):(n12+n22)])/prob.joint[(n12+1):(n12+n22)]*g22.est.ll
  centropy.ll.joint <- centropy(prob.joint, c(rep(0, n12), rep(1, n22)))
  cerror12.ll.joint <- mean(prob.joint[1:n12]>0.5)
  cerror22.ll.joint <- mean(prob.joint[(n12+1):(n12+n22)]<0.5)
  temp <- getFinalStat(g12.est.ll, g22.est.ll, v12.est.ll, v22.est.ll)
  u.ll <- temp$U
  var.ll1 <- temp$sigma.1
  llpvalue1 <- pnorm(temp$z1)
  var.ll2 <- temp$sigma.2
  llpvalue2 <- pnorm(temp$z2)
  llpvalue.gm <- pnorm(temp$z.gm)
  llpvalue.hm <- pnorm(temp$z.hm)
  hidden.layers <- c(10,10)
  learn.rates <- 0.001
  n.epochs <- 4000
  x.fit <- data.frame(x=rbind(x11, x21))
  newdata1 <- data.frame(x=x12)
  newdata2 <- data.frame(x=x22)
  temp <- NNfun(x.fit, label.fit, newdata1, newdata2, nnrep = 10, hidden.layers = hidden.layers,
                n.epochs = n.epochs, learn.rates = learn.rates)
  g12.est.nn <- temp$prob1.fit/(1-temp$prob1.fit)*n11/n21
  g22.est.nn <- temp$prob2.fit/(1-temp$prob2.fit)*n11/n21
  cerror12.nn.marginal <- mean(temp$prob1.fit>0.5)
  cerror22.nn.marginal <- mean(temp$prob2.fit<0.5)
  centropy.nn.marginal <- centropy(c(temp$prob1.fit, temp$prob2.fit), c(rep(0, n12), rep(1, n22)))
  error.nn <- gerror(g12.est.nn, g12.orac)
  xy.fit <- data.frame(x=rbind(x11, x21), y=c(y11, y21))
  newdata1 <- data.frame(x=x12, y=y12)
  newdata2 <- data.frame(x=x22, y=y22)
  temp <- NNfun(xy.fit, label.fit, newdata1, newdata2, nnrep = 10, hidden.layers = hidden.layers,
                n.epochs = n.epochs, learn.rates = learn.rates)
  v12.est.nn <- (1-temp$prob1.fit)/temp$prob1.fit*g12.est.nn
  v22.est.nn <- (1-temp$prob2.fit)/temp$prob2.fit*g22.est.nn
  cerror12.nn.joint <- mean(temp$prob1.fit>0.5)
  cerror22.nn.joint <- mean(temp$prob2.fit<0.5)
  centropy.nn.joint <- centropy(c(temp$prob1.fit, temp$prob2.fit), c(rep(0, n12), rep(1, n22)))
  temp <- getFinalStat(g12.est.nn, g22.est.nn, v12.est.nn, v22.est.nn)
  u.nn <- temp$U
  var.nn1 <- temp$sigma.1
  nnpvalue1 <- pnorm(temp$z1)
  var.nn2 <- temp$sigma.2
  nnpvalue2 <- pnorm(temp$z2)
  nnpvalue.gm <- pnorm(temp$z.gm)
  nnpvalue.hm <- pnorm(temp$z.hm)
  result <- list(nnpvalue1=nnpvalue1, nnpvalue2=nnpvalue2,
                 nnpvalue.gm=nnpvalue.gm, nnpvalue.hm=nnpvalue.hm,
                 llpvalue1=llpvalue1, llpvalue2=llpvalue2, 
                 llpvalue.gm=llpvalue.gm, llpvalue.hm=llpvalue.hm,
                 oracpvalue1=oracpvalue1, oracpvalue2=oracpvalue2,
                 oracpvalue.gm=oracpvalue.gm, oracpvalue.hm=oracpvalue.hm,
                 u.orac=u.orac, u.ll=u.ll, u.nn=u.nn,
                 var.nn1=var.nn1, var.nn2=var.nn2,
                 var.ll1=var.ll1, var.ll2=var.ll2,
                 centropy.ll.marginal=centropy.ll.marginal, centropy.ll.joint=centropy.ll.joint,
                 centropy.nn.marginal=centropy.nn.marginal, centropy.nn.joint=centropy.nn.joint,
                 error.ll=error.ll, error.nn=error.nn,
                 cerror12.ll.marginal=cerror12.ll.marginal, cerror22.ll.marginal=cerror22.ll.marginal,
                 cerror12.nn.marginal=cerror12.nn.marginal, cerror22.nn.marginal=cerror22.nn.marginal,
                 cerror12.ll.joint=cerror12.ll.joint, cerror22.ll.joint=cerror22.ll.joint,
                 cerror12.nn.joint=cerror12.nn.joint, cerror22.nn.joint=cerror22.nn.joint)
}



LM_generation <- function(n, beta_hat, cor, n_nuisance, eps, mean_X = 0, var_X = 1){
  p <- length(beta_hat)
  corr_matrix <- matrix(NA, p, p)
  for(i in 1:p){
    for(j in 1:p){
      corr_matrix[i, j] <- (cor ^ (abs(i - j))) * var_X
    }
  }
  X_design <- mvrnorm(n, mu = rep(mean_X, p), corr_matrix)
  Y1 <- as.matrix(X_design) %*% as.matrix(beta_hat, nrow = p)
  random_error <- rnorm(n, 0, eps)
  #adding the nuisance random noise features:
  X_nuiss <- matrix(0, n, n_nuisance)
  for(i in 1:n_nuisance){
    X_nuiss[,i] <- rnorm(n, 0, 1)
  }
  Y <- Y1 + random_error
  df_return <- data.frame(cbind(Y1, X_design, X_nuiss, Y))
  ncol_df <- ncol(df_return)
  colnames(df_return) <- c("Y1", paste("X", c(1:p), sep = ""), paste("X_nuis", c(1:n_nuisance), sep = ""), "Y")
  X_return <- df_return[,2:(ncol_df - 1)]
  return(list(df_return = df_return, X_return = X_return))
}






