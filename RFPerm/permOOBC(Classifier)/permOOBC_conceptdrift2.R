
source("utils.R")
set.seed(42)
cor_seq <- seq(0, 0.7, 0.1)
eps_seq <- c(0.1, 0.2, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 100)
beta_train <- c(1, 1, 1, 1, 1, 1, 1, 1)
beta_test <- c(1, 1, 1, 1, 1, 1, 6, 6)
for(k in 1:8){
  power_df <- data.frame(
  	eps = eps_seq,
  	power = numeric(14)
  	)
  for(j in 1:14){
  	pval_list <- numeric(500)
  	for(i in 1:500){
  	  df1 <- LM_generation(n = 1000, beta_hat = beta_train, cor = cor_seq[k],
        n_nuisance = 15,
  	  	eps = eps_seq[j])
  	  df2 <- LM_generation(n = 1000, beta_hat = beta_test, cor = cor_seq[k],
        n_nuisance = 15,
  	  	eps = eps_seq[j])
  	  pval_list[i] <- permoob_classification(df1, df2)
  	}
  	power_df$power[j] <- mean(pval_list < 0.05)
  }
  write.csv(power_df, paste("permOOBC_conceptdrift2_cor", cor_seq[k], ".csv", sep = "_"))
}
