#benchmark wrapper1

library(tidyr)
library(dplyr)
library(readr)
library(tidyverse)
library(Ecume)
library(highmean)
library(highDmean)
library(RandomProjectionTest)
library(mltools)
library(kernlab)

benchmark_method_compa <- function(df1, df2, B){
  p_val_oob_list <- numeric(B)
  p_val_cf_list <- numeric(B)
  p_val_kmmd_list <- numeric(B)
  p_val_c2st_list <- numeric(B)
  p_val_miles16_list <- numeric(B)
  p_val_xu16_list <- numeric(B)
  p_val_pan14_list <- numeric(B)
  p_val_cai14_list <- numeric(B)
  p_val_chen10_list <- numeric(B)
  p_val_chen14_list <- numeric(B)
  p_val_sri08_list <- numeric(B)
  p_val_bai96_list <- numeric(B)
  p_val_skk_list <- numeric(B)
  p_val_zwl_list <- numeric(B)
  p_val_zwlm_list <- numeric(B)
  p_val_miles16_list <- numeric(B)
  p_val_xgb <- numeric(B)
  conformal_result <- numeric(B)
  for(i in 1:B){
    p_val_kmmd_list[i] <- as.numeric(H0(kmmd(as.matrix(df1),
      as.matrix(df2), kernel = "rbfdot")))
    p_val_c2st_list[i] <- as.numeric(classifier_test(as.matrix(df1),
      as.matrix(df2)
      )$p.value)
    p_val_chen10_list[i] <- epval_Chen2010(as.matrix(df1), 
      as.matrix(df2))$pval
    p_val_chen14_list[i] <- epval_Chen2014(as.matrix(df1), 
      as.matrix(df2))$pval
    p_val_sri08_list[i] <- epval_Sri2008(as.matrix(df1), 
      as.matrix(df2))$pval
    p_val_bai96_list[i] <- epval_Bai1996(df1, df2)$pval
    p_val_cai14_list[i] <- epval_Cai2014(df1, df2)$pval
    p_val_skk_list[i] <- SKK_test(df1, df2)$pvalue
    p_val_zwl_list[i] <- zwl_test(as.matrix(df1), as.matrix(df2), order = 0)$pvalue
    p_val_zwlm_list[i] <- zwl_test(as.matrix(df1), as.matrix(df2), order = 2)$pvalue 
    p_val_pan14_list[i] <- min(epval_aSPU(as.matrix(df1), 
      as.matrix(df2))$pval)
    p_val_cf_list[i] <- cfperm2(as.data.frame(df1), as.data.frame(df2), n_perm = 150)
    p_val_xu16_list[i] <- min(apval_aSPU(as.matrix(df1), 
      as.matrix(df2), bandwidth = 10)$pval)
    p_val_miles16_list[i] <- random_projection_test(as.matrix(df1), as.matrix(df2))$p_value
    p_val_xgb[i] <- permXGB_nontuning(as.data.frame(df1), as.data.frame(df2), B = 5000)$p_val
  }
  return(list(
    power_kmmd = mean(p_val_kmmd_list),
    power_c2st = mean(p_val_c2st_list < 0.05),
    power_cf = mean(p_val_cf_list),
    power_oob  = mean(p_val_oob_list < 0.05),
    power_xu16  = mean(p_val_xu16_list < 0.05),
    power_chen10  = mean(p_val_chen10_list < 0.05),
    power_chen14 = mean(p_val_chen14_list < 0.05),
    power_sri08  = mean(p_val_sri08_list < 0.05),
    power_bai96 = mean(p_val_bai96_list < 0.05),
    power_cai14 = mean(p_val_skk_list < 0.05),
    power_skk = mean(p_val_zwl_list < 0.05),
    power_zwlm = mean(p_val_zwlm_list < 0.05),
    power_zwl = mean(p_val_cai14_list < 0.05),
    power_pan14 = mean(p_val_pan14_list < 0.05),
    power_xgb = mean(p_val_xgb < 0.05)
  ))
}













