# RFPerm

`RFPerm` provides permutation-based hypothesis tests for detecting dataset shift
by contrasting predictive errors (e.g., training out-of-bag error vs. test-batch
prediction error) using random forests.

## Installation (local)

From the package root (the folder containing `DESCRIPTION`):

```r
devtools::install()
```

Or from a built tarball:

```bash
R CMD INSTALL RFPerm_0.1.0.tar.gz
```

## Quick start

```r
# Regression RFPerm (OOB vs test MSE)
set.seed(1)
n_tr <- 200
n_te <- 200
df_train <- data.frame(X1=rnorm(n_tr), X2=rnorm(n_tr))
df_train$Y <- df_train$X1 - 0.5*df_train$X2 + rnorm(n_tr)

df_test <- data.frame(X1=rnorm(n_te), X2=rnorm(n_te))
df_test$Y <- df_test$X1 - 0.5*df_test$X2 + rnorm(n_te)

p <- perm_oobtest(df_train, df_test, B=500)
p
```

## Development

- Run tests: `devtools::test()`
- Full check: `devtools::check()`

## Notes

If you use **roxygen2** documentation, regenerate the NAMESPACE and man files via:

```r
devtools::document()
```
