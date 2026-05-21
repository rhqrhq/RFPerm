#conformal test wrapper:
source('conformal_hu2023_benchmark.R')


## compute the oracle density ratio
g.orac <- function(x){
  x = matrix(x, ncol=p)
  exp(x[,1:4] %*% c(1,1,-1,-1))/exp(2)
}

conformal_test <- function(df1, df2, mc, alpha = 0.05){
  #Generate the coefficients
  beta <- as.matrix(replicate(s, sample(c(1,-1), 1, prob = c(0.5, 0.5))))
  u.orac <- matrix(0, mc, 3)
  u.ll <- matrix(0, mc, 3); u.ql <- matrix(0, mc, 3)
  u.klr <- matrix(0, mc, 3); u.nn <- matrix(0, mc, 3)
  var.orac1 <- matrix(0, mc, 9) # record results under population 1
  var.ll1 <- matrix(0, mc, 9); var.ql1 <- matrix(0, mc, 9)
  var.klr1 <- matrix(0, mc, 9); var.nn1 <- matrix(0, mc, 9)
  var.orac2 <- matrix(0, mc, 3) # record results under population 2
  var.ll2 <- matrix(0, mc, 3); var.ql2 <- matrix(0, mc, 3)
  var.klr2 <- matrix(0, mc, 3); var.nn2 <- matrix(0, mc, 3)
  var.orac.hm <- matrix(0, mc, 9) # record results under harmonic mean
  var.ll.hm <- matrix(0, mc, 9); var.ql.hm <- matrix(0, mc, 9)
  var.klr.hm <- matrix(0, mc, 9); var.nn.hm <- matrix(0, mc, 9)
  oracpvalue.hm <- matrix(0, mc, 9)
  llpvalue.hm <- matrix(0, mc, 9); qlpvalue.hm <- matrix(0, mc, 9)
  klrpvalue.hm <- matrix(0, mc, 9); nnpvalue.hm <- matrix(0, mc, 9)
  centropy.ll.joint <- matrix(0, mc, 3); centropy.ql.joint <- matrix(0, mc, 3)
  centropy.klr.joint <- matrix(0, mc, 3); centropy.nn.joint <- matrix(0, mc, 3)
  centropy.ll.marginal <- matrix(0, mc, 3); centropy.ql.marginal <- matrix(0, mc, 3)
  centropy.klr.marginal <- matrix(0, mc, 3); centropy.nn.marginal <- matrix(0, mc, 3)
  error.ll <- matrix(0, mc, 3); error.ql <- matrix(0, mc, 3)
  error.klr <- matrix(0, mc, 3); error.nn <- matrix(0, mc, 3)
  errStat.ll <- errStat.ql <- errStat.nn <- errStat.klr <- array(0, dim=c(mc,3,7))
  lrpvalue <- rep(0, mc); lrpvalue_split <- matrix(0, mc, 3)
  p_X <- ncol(df1) - 1
  n1 <- nrow(df1)
  n2 <- nrow(df2)
  prop <- c(0.3, 0.5, 0.8)
  #Doing the bootstrap resampling for re - estimate:
  for(run in 1:mc){
    df1_b <- sample(df1, n1, replace = TRUE)
    df2_b <- sample(df2, n2, replace = TRUE)
    x1 <- as.matrix(df1_b[,1:p_X])
    y1 <- as.matrix(df1_b[,(p_X + 1)])
    x2 <- as.matrix(df2_b[,1:p_X])
    y2 <- as.matrix(df2_b[,(p_X + 1)])
    lrpvalue[run] <- LR(x1, x2, y1, y2)
    for(i.prop in 1:3){
      n12 <- ceiling(n1*prop[i.prop]); n22 <- ceiling(n2*prop[i.prop])
      n11 <- n1 - n12; n21 <- n2 - n22
      x11 <- x1[1:n11, , drop=F]; x12 <- x1[-(1:n11),,drop=F]
      y11 <- y1[1:n11]; y12 <- y1[-(1:n11)]
      x21 <- x2[1:n21, , drop=F]; x22 <- x2[-(1:n21),,drop=F]
      y21 <- y2[1:n21]; y22 <- y2[-(1:n21)]  
      lrpvalue_split[run, i.prop] <- LR(x12, x22, y12, y22)
      g12.orac <- g.orac(x12); v12.orac <- rep(1, n12)
      g22.orac <- g.orac(x22); v22.orac <- rep(1, n22)
      temp <- getFinalStat(g12.orac, g22.orac, v12.orac, v22.orac)
      u.orac[run, i.prop] <- temp$U
      var.orac1[run, ((i.prop-1)*3+1):(i.prop*3)] <- temp$sigma.1
      var.orac2[run, i.prop] <- temp$sigma.2
      oracpvalue.hm[run, ((i.prop-1)*3+1):(i.prop*3)] <- pnorm(temp$z.hm)
      print(paste0('run: ', run, 'fitting: linear logistic'))
      label.fit <- as.factor(c(rep(0,n11), rep(1,n21)))
      xy.fit <- data.frame(x.fit=rbind(x11, x21), y.fit=c(y11, y21))
      fit.joint <- glm(label.fit~., data=xy.fit, family="binomial")
      x.fit <- data.frame(x.fit=rbind(x11,x21))
      fit.marginal <- glm(label.fit~., data=x.fit, family="binomial")
      print(paste0('run: ', run, 'estimating: linear logistic'))
      prob.marginal <- predict(fit.marginal, newdata=data.frame(x.fit=rbind(x12,x22)), type="response")
      prob.marginal[prob.marginal<0.01] <- 0.01; prob.marginal[prob.marginal>0.99] <- 0.99
      g12.est.ll <- prob.marginal[1:n12]/(1-prob.marginal[1:n12])*n11/n21
      g22.est.ll <- prob.marginal[(n12+1):(n12+n22)]/(1-prob.marginal[(n12+1):(n12+n22)])*n11/n21
      centropy.ll.marginal[run,i.prop] <- centropy(prob.marginal, c(rep(0, n12), rep(1, n22)))
      error.ll[run, i.prop] <- gerror(g12.est.ll, g12.orac)
      prob.joint <- predict(fit.joint, newdata=data.frame(x.fit=rbind(x12,x22), y.fit=c(y12,y22)), type="response")
      prob.joint[prob.joint<0.01] <- 0.01; prob.joint[prob.joint>0.99] <- 0.99
      v12.est.ll <- (1-prob.joint[1:n12])/prob.joint[1:n12]*g12.est.ll
      v22.est.ll <- (1-prob.joint[(n12+1):(n12+n22)])/prob.joint[(n12+1):(n12+n22)]*g22.est.ll
      centropy.ll.joint[run,i.prop] <- centropy(prob.joint, c(rep(0, n12), rep(1, n22)))
      temp <- getFinalStat(g12.est.ll, g22.est.ll, v12.est.ll, v22.est.ll)
      u.ll[run, i.prop] <- temp$U
      var.ll1[run, ((i.prop-1)*3+1):(i.prop*3)] <- temp$sigma.1
      var.ll2[run, i.prop] <- temp$sigma.2
      llpvalue.hm[run, ((i.prop-1)*3+1):(i.prop*3)] <- pnorm(temp$z.hm)
      errStat.ll[run,i.prop,] <- getCor(g12.est.ll, v12.est.ll, v22.est.ll, g12.orac, v12.orac, v22.orac)
      print(paste0('run: ', run, 'estimating: KLR'))
      xy.fit <- cbind(rbind(x11, x21), c(y11, y21))
      data.fit <- constructData(xy.fit, label.fit)
      klrlearner <- constructKlogRegLearner()
      params <- list(kernel='rbfdot', sigma=0.005, lambda=0.05/getN(data.fit), tol=10e-6, maxiter=500)
      fit.joint <- klrlearner$learn(data.fit, params)
      x.fit <- rbind(x11, x21)
      data.fit <- constructData(x.fit, label.fit)
      klrlearner <- constructKlogRegLearner()
      params <- list(kernel='rbfdot', sigma=0.005, lambda=0.05/getN(data.fit), tol=10e-6, maxiter=500)
      fit.marginal <- klrlearner$learn(data.fit, params)
      newdata <- rbind(x12, x22)
      K = kernelMult(fit.marginal$kernel, newdata, fit.marginal$data, fit.marginal$alpha)
      pi = 1 / (1 + exp(-as.vector(K))) # predicted probabilities
      pi[pi<0.01] <- 0.01; pi[pi>0.99] <- 0.99
      g12.est.klr <- pi[1:n12]/(1-pi[1:n12])*n11/n21
      g22.est.klr <- pi[(n12+1):(n12+n22)]/(1-pi[(n12+1):(n12+n22)])*n11/n21
      centropy.klr.marginal[run,i.prop] <- centropy(pi, c(rep(0, n12), rep(1, n22)))
      error.klr[run, i.prop] <- gerror(g12.est.klr, g12.orac)
      newdata <- cbind(rbind(x12, x22), c(y12, y22))
      K = kernelMult(fit.joint$kernel, newdata, fit.joint$data, fit.joint$alpha)
      pi = 1 / (1 + exp(-as.vector(K))) # predicted probabilities
      pi[pi<0.01] <- 0.01; pi[pi>0.99] <- 0.99
      v12.est.klr <- (1-pi[1:n12])/pi[1:n12]*g12.est.klr
      v22.est.klr <- (1-pi[(n12+1):(n12+n22)])/pi[(n12+1):(n12+n22)]*g22.est.klr
      centropy.klr.joint[run,i.prop] <- centropy(pi, c(rep(0, n12), rep(1, n22)))
      temp <- getFinalStat(g12.est.klr, g22.est.klr, v12.est.klr, v22.est.klr)
      u.klr[run, i.prop] <- temp$U
      var.klr1[run, ((i.prop-1)*3+1):(i.prop*3)] <- temp$sigma.1
      var.klr2[run, i.prop] <- temp$sigma.2
      klrpvalue.hm[run,((i.prop-1)*3+1):(i.prop*3)] <- pnorm(temp$z.hm)
      errStat.klr[run,i.prop,] <- getCor(g12.est.klr, v12.est.klr, v22.est.klr, g12.orac, v12.orac, v22.orac)
      print(paste0('run: ', run, 'estimating: NN'))
      hidden.layers <- c(10,10)
      learn.rates <- 0.001
      n.epochs <- 200
      x.fit <- data.frame(x=rbind(x11, x21))
      newdata1 <- data.frame(x=x12)
      newdata2 <- data.frame(x=x22)
      temp <- NNfun(x.fit, label.fit, newdata1, newdata2, nnrep = 10, hidden.layers = hidden.layers,
                    n.epochs = n.epochs, learn.rates = learn.rates)
      g12.est.nn <- temp$prob1.fit/(1-temp$prob1.fit)*n11/n21
      g22.est.nn <- temp$prob2.fit/(1-temp$prob2.fit)*n11/n21
      centropy.nn.marginal[run,i.prop] <- centropy(c(temp$prob1.fit, temp$prob2.fit), c(rep(0, n12), rep(1, n22)))
      error.nn[run, i.prop] <- gerror(g12.est.nn, g12.orac)
      hidden.layers <- c(10,10)
      learn.rates <- 0.001
      n.epochs <- 200
      xy.fit <- data.frame(x=rbind(x11, x21), y=c(y11, y21))
      newdata1 <- data.frame(x=x12, y=y12)
      newdata2 <- data.frame(x=x22, y=y22)
      temp <- NNfun(xy.fit, label.fit, newdata1, newdata2, nnrep = 10, hidden.layers = hidden.layers,
                    n.epochs = n.epochs, learn.rates = learn.rates)
      v12.est.nn <- (1-temp$prob1.fit)/temp$prob1.fit*g12.est.nn
      v22.est.nn <- (1-temp$prob2.fit)/temp$prob2.fit*g22.est.nn
      centropy.nn.joint[run,i.prop] <- centropy(c(temp$prob1.fit, temp$prob2.fit), c(rep(0, n12), rep(1, n22)))
      temp <- getFinalStat(g12.est.nn, g22.est.nn, v12.est.nn, v22.est.nn)
      u.nn[run, i.prop] <- temp$U
      var.nn1[run, ((i.prop-1)*3+1):(i.prop*3)] <- temp$sigma.1
      var.nn2[run, i.prop] <- temp$sigma.2
      nnpvalue.hm[run,((i.prop-1)*3+1):(i.prop*3)] <- pnorm(temp$z.hm)
      errStat.nn[run,i.prop,] <- getCor(g12.est.nn, v12.est.nn, v22.est.nn, g12.orac, v12.orac, v22.orac)
  	}
  }
  oracsize.hm <- colMeans(oracpvalue.hm < alpha)
  llsize.hm <- colMeans(llpvalue.hm < alpha)
  #qlsize.hm <- colMeans(qlpvalue.hm < alpha)
  klrsize.hm <- colMeans(klrpvalue.hm < alpha)
  nnsize.hm <- colMeans(nnpvalue.hm < alpha)
  lrsize_split <- colMeans(lrpvalue_split < alpha)
  lrsize <- mean(lrpvalue < alpha)
  size <- list(oracsize.hm = oracsize.hm,
               llsize.hm=llsize.hm, 
               klrsize.hm=klrsize.hm, nnsize.hm=nnsize.hm,
               lrsize=lrsize, lrsize_split=lrsize_split)
  return(size)
}











