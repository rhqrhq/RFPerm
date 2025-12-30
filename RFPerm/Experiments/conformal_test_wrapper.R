#Benchmark for the Hu2023(JASA)
library(MASS)
library(glmnet)
library(ANN2)
library(CVST)
library(kernlab)


sim_fun <- function(x1, y1, x2, y2, n11, n12, n21, n22){
  
  n1 <- n11+n12; n2 <- n21+n22
  index1 <- sample(1:n1, size = n11)
  x11 <- x1[index1,]; x12 <- x1[-index1,]
  y11 <- y1[index1]; y12 <- y1[-index1]
  index2 <- sample(1:n2, size = n21)
  x21 <- x2[index2,]; x22 <- x2[-index2,]
  y21 <- y2[index2]; y22 <- y2[-index2]
  
  #print(paste0('run: ', run, 'fitting: linear logistic'))
  label.fit <- as.factor(c(rep(0,n11), rep(1,n21)))
  xy.fit <- data.frame(x.fit=rbind(x11, x21), y.fit=c(y11, y21))
  fit.joint <- glm(label.fit~., data=xy.fit, family="binomial")
  x.fit <- data.frame(x.fit=rbind(x11,x21))
  fit.marginal <- glm(label.fit~., data=x.fit, family="binomial")
  
  #print(paste0('run: ', run, 'estimating: linear logistic'))
  prob.marginal <- predict(fit.marginal, newdata=data.frame(x.fit=rbind(x12,x22)), type="response")
  prob.marginal[prob.marginal<0.01] <- 0.01; prob.marginal[prob.marginal>0.99] <- 0.99
  g12.est.ll <- prob.marginal[1:n12]/(1-prob.marginal[1:n12])*n11/n21
  g22.est.ll <- prob.marginal[(n12+1):(n12+n22)]/(1-prob.marginal[(n12+1):(n12+n22)])*n11/n21
  cerror12.ll.marginal <- mean(prob.marginal[1:n12]>0.5)
  cerror22.ll.marginal <- mean(prob.marginal[(n12+1):(n12+n22)]<0.5)
  centropy.ll.marginal <- centropy(prob.marginal, c(rep(0, n12), rep(1, n22)))
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
  
  #print(paste0('run: ', run, 'estimating: NN'))
  #starting dimensions:
  p <- ncol(x11)
  hidden.layers <- c(max(round(sqrt(p)),5), 5)
  learn.rates <- 0.001
  n.epochs <- 2000
  
  x.fit <- data.frame(x=rbind(x11, x21))
  newdata1 <- data.frame(x=x12)
  newdata2 <- data.frame(x=x22)
  temp <- NNfun(x.fit, label.fit, newdata1, newdata2, nnrep = 5, hidden.layers = hidden.layers,
                n.epochs = n.epochs, learn.rates = learn.rates)
  g12.est.nn <- temp$prob1.fit/(1-temp$prob1.fit)*n11/n21
  g22.est.nn <- temp$prob2.fit/(1-temp$prob2.fit)*n11/n21
  cerror12.nn.marginal <- mean(temp$prob1.fit>0.5)
  cerror22.nn.marginal <- mean(temp$prob2.fit<0.5)
  centropy.nn.marginal <- centropy(c(temp$prob1.fit, temp$prob2.fit), c(rep(0, n12), rep(1, n22)))

  xy.fit <- data.frame(x=rbind(x11, x21), y=c(y11, y21))
  newdata1 <- data.frame(x=x12, y=y12)
  newdata2 <- data.frame(x=x22, y=y22)
  temp <- NNfun(xy.fit, label.fit, newdata1, newdata2, nnrep = 5, hidden.layers = hidden.layers,
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
                 u.ll=u.ll, u.nn=u.nn,
                 var.nn1=var.nn1, var.nn2=var.nn2,
                 var.ll1=var.ll1, var.ll2=var.ll2,
                 centropy.ll.marginal=centropy.ll.marginal, centropy.ll.joint=centropy.ll.joint,
                 centropy.nn.marginal=centropy.nn.marginal, centropy.nn.joint=centropy.nn.joint,
                 cerror12.ll.marginal=cerror12.ll.marginal, cerror22.ll.marginal=cerror22.ll.marginal,
                 cerror12.nn.marginal=cerror12.nn.marginal, cerror22.nn.marginal=cerror22.nn.marginal,
                 cerror12.ll.joint=cerror12.ll.joint, cerror22.ll.joint=cerror22.ll.joint,
                 cerror12.nn.joint=cerror12.nn.joint, cerror22.nn.joint=cerror22.nn.joint)

  return(result)
}





