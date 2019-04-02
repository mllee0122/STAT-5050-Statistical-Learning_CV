#Statistical Learning_Homework3
#107064522

##### 1. #####
library(MASS)
set.seed(36)
n <- 200; sigma <- 2; beta0 <- c(2, -2, 0.5, 1, -3, 0, 0, 0, 0, 0)
cormat <- diag(1, nrow = 10, ncol = 10); cormat[cormat == 0] <- 0.5
cholmat <- chol(sigma*cormat)
x <- matrix(rnorm(10*n, 0, 1), ncol = 10)%*%cholmat
err <- rnorm(n, 0, sigma)
y <- x%*%beta0+err

##### 2.(a) #####
library(glmnet)
beta_lasso = glmnet(x, y, alpha = 1,
                    lambda = seq(3, 0, length = 200),
                    standardize = F)

plot(seq(3, 0, length = 200), beta_lasso$beta[1,],
     xlim = c(0, 3),
     ylim = c(-3, 2),
     xlab = expression(lambda),
     ylab = expression(beta(lambda)), 
     main = 'Solution Path',
     type = 'l', lwd = 2, lty = 1, col = 1)
for(i in 2:5){
  lines(seq(3, 0, length = 200), beta_lasso$beta[i,], type = 'l', lwd = 2, lty = 1, col = i)
}
for(i in 6:10){
  lines(seq(3, 0, length = 200), beta_lasso$beta[i,], type = 'l', lwd = 2, lty = 2, col = i)
}
legend(2.4, 2.2, c(expression(paste(beta^"*", "= 2"), paste(beta^"*", "= -2"), paste(beta^"*", "= 0.5"), paste(beta^"*", "= 1"), paste(beta^"*", "= -3"))), 
       lwd = 2, lty = 1, col = c(1, 2, 3, 4, 5), cex = 1, bty = "n")
# AIC
sigma_hat_sq = colSums(((y-x%*%beta_lasso$beta)^2))/n
AIC = n*log(sigma_hat_sq)+2*beta_lasso$df
# BIC
BIC = n*log(sigma_hat_sq)+beta_lasso$df*log(n)
# Cross Validation
beta_cv = cv.glmnet(x, y, alpha = 1,
                    lambda = seq(3, 0, length = 200),
                    standardize = F,
                    nfolds = 5)
CV = beta_cv$cvm
# GCV
GCV = sigma_hat_sq*(1+2*beta_lasso$df/n)
# Choose the minimum value
min_AIC = which(AIC == min(AIC))
min_BIC = which(BIC == min(BIC))
min_CV = which(CV == min(CV))
min_GCV = which(GCV == min(GCV))
# lambda
lambda_AIC = seq(3, 0, length = 200)[min_AIC]
lambda_BIC = seq(3, 0, length = 200)[min_BIC]
lambda_CV = seq(3, 0, length = 200)[min_CV]
lambda_GCV = seq(3, 0, length = 200)[min_GCV]
# beta of lambda
beta_AIC = beta_lasso$beta[, min_AIC]
beta_BIC = beta_lasso$beta[, min_BIC]
beta_CV = beta_lasso$beta[, min_CV]
beta_GCV = beta_lasso$beta[, min_GCV]
# beta of lambda - beta0
err_AIC = sqrt(sum((beta_AIC-beta0)^2))
err_BIC = sqrt(sum((beta_BIC-beta0)^2))
err_CV = sqrt(sum((beta_CV-beta0)^2))
err_GCV = sqrt(sum((beta_GCV-beta0)^2))

sum(beta_AIC[1:5] != 0)
sum(beta_BIC[1:5] != 0)
sum(beta_CV[1:5] != 0)
sum(beta_GCV[1:5] != 0)

sum(beta_AIC != 0)
sum(beta_BIC != 0)
sum(beta_CV != 0)
sum(beta_GCV != 0)

# Comparison for 4.

abline(v = lambda_AIC, col = 1, lty = 6, lwd = 1)
abline(v = lambda_BIC, col = 2, lty = 6, lwd = 1)
abline(v = lambda_CV, col = 3, lty = 6, lwd = 1)
abline(v = lambda_GCV, col = 4, lty = 6, lwd = 1)

legend("bottomright", c(expression(paste("optimal ", lambda, "(AIC)"), paste("optimal ", lambda, "(BIC)"), paste("optimal ", lambda, "(CV)"), paste("optimal ", lambda, "(GCV)"))), 
       lwd = 1, lty = 6, col = c(1, 2, 3, 4), cex = 1, bty = "n")

fit = which(beta_lasso$df == 5)
lambda_fit = seq(3, 0, length = 200)[fit]
range(lambda_fit)

##### 2.(b) #####
err_AIC_b = 0; err_BIC_b = 0; err_CV_b = 0; err_GCV_b = 0; 
sure_AIC = 0; sure_BIC = 0; sure_CV = 0; sure_GCV = 0; 
exact_AIC = 0; exact_BIC = 0; exact_CV = 0; exact_GCV = 0; 
size_AIC = 0; size_BIC = 0; size_CV = 0; size_GCV = 0; 

for (i in 1:100){
  set.seed(i)
  x_b <- matrix(rnorm(10*n, 0, 1), ncol = 10)%*%cholmat
  err <- rnorm(n, 0, sigma)
  y_b <- x_b%*%beta0+err
  
  beta_lasso_b = glmnet(x_b, y_b, alpha = 1,
                      lambda = seq(3, 0, length = 200),
                      standardize = F)
  # AIC
  sigma_hat_sq_b = colSums(((y_b-x_b%*%beta_lasso_b$beta)^2))/n
  AIC_b = n*log(sigma_hat_sq_b)+2*beta_lasso_b$df
  # BIC
  BIC_b = n*log(sigma_hat_sq_b)+beta_lasso_b$df*log(n)
  # Cross Validation
  beta_cv_b = cv.glmnet(x_b, y_b, alpha = 1,
                      lambda = seq(3, 0, length = 200),
                      standardize = F,
                      nfolds = 5)
  CV_b = beta_cv_b$cvm
  # GCV
  GCV_b = sigma_hat_sq_b*(1+2*beta_lasso_b$df/n)
  # Choose the minimum value
  min_AIC_b = which(AIC_b == min(AIC_b))
  min_BIC_b = which(BIC_b == min(BIC_b))
  min_CV_b = which(CV_b == min(CV_b))
  min_GCV_b = which(GCV_b == min(GCV_b))
  # lambda
  lambda_AIC_b = seq(3, 0, length = 200)[min_AIC_b]
  lambda_BIC_b = seq(3, 0, length = 200)[min_BIC_b]
  lambda_CV_b = seq(3, 0, length = 200)[min_CV_b]
  lambda_GCV_b = seq(3, 0, length = 200)[min_GCV_b]
  # beta of lambda
  beta_AIC_b = beta_lasso_b$beta[, min_AIC_b]
  beta_BIC_b = beta_lasso_b$beta[, min_BIC_b]
  beta_CV_b = beta_lasso_b$beta[, min_CV_b]
  beta_GCV_b = beta_lasso_b$beta[, min_GCV_b]
  # beta of lambda - beta0
  err_AIC_b = err_AIC_b+sqrt(sum((beta_AIC_b-beta0)^2))
  err_BIC_b = err_BIC_b+sqrt(sum((beta_BIC_b-beta0)^2))
  err_CV_b = err_CV_b+sqrt(sum((beta_CV_b-beta0)^2))
  err_GCV_b = err_GCV_b+sqrt(sum((beta_GCV_b-beta0)^2))
  
  # Sure and Exact
  if (sum(beta_AIC_b[1:5]!=0) == 5){
    sure_AIC = sure_AIC+1
  }
  if (sum(beta_AIC_b[1:5]!=0) == 5 & sum(beta_AIC_b[6:10]!=0) == 0){
    exact_AIC = exact_AIC+1
  }
  if (sum(beta_BIC_b[1:5]!=0) == 5){
    sure_BIC = sure_BIC+1
  }
  if (sum(beta_BIC_b[1:5]!=0) == 5 & sum(beta_BIC_b[6:10]!=0) == 0){
    exact_BIC = exact_BIC+1
  }
  if (sum(beta_CV_b[1:5]!=0) == 5){
    sure_CV = sure_CV+1
  }
  if (sum(beta_CV_b[1:5]!=0) == 5 & sum(beta_CV_b[6:10]!=0) == 0){
    exact_CV = exact_CV+1
  }
  if (sum(beta_GCV_b[1:5]!=0) == 5){
    sure_GCV = sure_GCV+1
  }
  if (sum(beta_GCV_b[1:5]!=0) == 5 & sum(beta_GCV_b[6:10]!=0) == 0){
    exact_GCV = exact_GCV+1
  }
  # Size
  size_AIC = size_AIC+sum(beta_AIC_b!=0)
  size_BIC = size_BIC+sum(beta_BIC_b!=0)
  size_CV = size_CV+sum(beta_CV_b!=0)
  size_GCV = size_GCV+sum(beta_GCV_b!=0)
}
# (1) average (beta of lambda - beta0)
err_AIC_b = err_AIC_b/100
err_BIC_b = err_BIC_b/100
err_CV_b = err_CV_b/100
err_GCV_b = err_GCV_b/100

# (2) Sure
sure_AIC = sure_AIC/100
sure_BIC = sure_BIC/100
sure_CV = sure_CV/100
sure_GCV = sure_GCV/100

# (3) Exact
exact_AIC = exact_AIC/100
exact_BIC = exact_BIC/100
exact_CV = exact_CV/100
exact_GCV = exact_GCV/100

# (4) Size
size_AIC = size_AIC/100
size_BIC = size_BIC/100
size_CV = size_CV/100
size_GCV = size_GCV/100

##### 3.(a) #####
set.seed(36)
n <- 200; sigma <- 2; beta0 <- c(2, -2, 0.5, 1, -3, rep(0,995))
cormat <- diag(1, nrow = 1000, ncol = 1000); cormat[cormat == 0] <- 0.5
cholmat <- chol(sigma*cormat)
x <- matrix(rnorm(1000*n, 0, 1), ncol = 1000)%*%cholmat
err <- rnorm(n, 0, sigma)
y <- x%*%beta0+err  

beta_lasso = glmnet(x, y, alpha = 1,
                    lambda = seq(3, 0, length = 200),
                    standardize = F)

plot(seq(3, 0, length = 200), beta_lasso$beta[1,],
     xlim = c(0, 3),
     ylim = c(-3, 2),
     xlab = expression(lambda),
     ylab = expression(beta(lambda)), 
     main = 'Solution Path',
     type = 'l', lwd = 2, lty = 1, col = 1)
for(i in 2:5){
  lines(seq(3, 0, length = 200), beta_lasso$beta[i,], type = 'l', lwd = 2, lty = 1, col = i)
}
for(i in 6:1000){
  lines(seq(3, 0, length = 200), beta_lasso$beta[i,], type = 'l', lwd = 2, lty = 2, col = i)
}
legend(2.4, 2.2, c(expression(paste(beta^"*", "= 2"), paste(beta^"*", "= -2"), paste(beta^"*", "= 0.5"), paste(beta^"*", "= 1"), paste(beta^"*", "= -3"))), 
       lwd = 2, lty = 1, col = c(1, 2, 3, 4, 5), cex = 1, bty = "n")

# AIC
sigma_hat_sq = colSums(((y-x%*%beta_lasso$beta)^2))/n
AIC = n*log(sigma_hat_sq)+2*beta_lasso$df
# BIC
BIC = n*log(sigma_hat_sq)+beta_lasso$df*log(n)
# Cross Validation
beta_cv = cv.glmnet(x, y, alpha = 1,
                    lambda = seq(3, 0, length = 200),
                    standardize = F,
                    nfolds = 5)
CV = beta_cv$cvm
# GCV
GCV = sigma_hat_sq*(1+2*beta_lasso$df/n)
# Choose the minimum value
min_AIC = which(AIC == min(AIC))
min_BIC = which(BIC == min(BIC))
min_CV = which(CV == min(CV))
min_GCV = which(GCV == min(GCV))
# lambda
lambda_AIC = seq(3, 0, length = 200)[min_AIC]
lambda_BIC = seq(3, 0, length = 200)[min_BIC]
lambda_CV = seq(3, 0, length = 200)[min_CV]
lambda_GCV = seq(3, 0, length = 200)[min_GCV]
# beta of lambda
beta_AIC = beta_lasso$beta[, min_AIC]
beta_BIC = beta_lasso$beta[, min_BIC]
beta_CV = beta_lasso$beta[, min_CV]
beta_GCV = beta_lasso$beta[, min_GCV]
# beta of lambda - beta0
err_AIC = sqrt(sum((beta_AIC-beta0)^2))
err_BIC = sqrt(sum((beta_BIC-beta0)^2))
err_CV = sqrt(sum((beta_CV-beta0)^2))
err_GCV = sqrt(sum((beta_GCV-beta0)^2))

sum(beta_AIC[1:5] != 0)
sum(beta_BIC[1:5] != 0)
sum(beta_CV[1:5] != 0)
sum(beta_GCV[1:5] != 0)

sum(beta_AIC != 0)
sum(beta_BIC != 0)
sum(beta_CV != 0)
sum(beta_GCV != 0)

# Comparison for 4.

abline(v = lambda_AIC, col = 1, lty = 6, lwd = 1)
abline(v = lambda_BIC, col = 2, lty = 6, lwd = 1)
abline(v = lambda_CV, col = 3, lty = 6, lwd = 1)
abline(v = lambda_GCV, col = 4, lty = 6, lwd = 1)

legend("bottomright", c(expression(paste("optimal ", lambda, "(AIC)"), paste("optimal ", lambda, "(BIC)"), paste("optimal ", lambda, "(CV)"), paste("optimal ", lambda, "(GCV)"))), 
       lwd = 1, lty = 6, col = c(1, 2, 3, 4), cex = 1, bty = "n")

fit = which(beta_lasso$df >= 4 & beta_lasso$df <= 6)
lambda_fit = seq(3, 0, length = 200)[fit]
range(lambda_fit)

##### 3.(b) #####
err_AIC_3b = 0; err_BIC_3b = 0; err_CV_3b = 0; err_GCV_3b = 0; 
sure_AIC_3 = 0; sure_BIC_3 = 0; sure_CV_3 = 0; sure_GCV_3 = 0; 
exact_AIC_3 = 0; exact_BIC_3 = 0; exact_CV_3 = 0; exact_GCV_3 = 0; 
size_AIC_3 = 0; size_BIC_3 = 0; size_CV_3 = 0; size_GCV_3 = 0; 

for (i in 1:100){
  set.seed(i)
  x_3b <- matrix(rnorm(1000*n, 0, 1), ncol = 1000)%*%cholmat
  err <- rnorm(n, 0, sigma)
  y_3b <- x_3b%*%beta0+err
  
  beta_lasso_3b = glmnet(x_3b, y_3b, alpha = 1,
                         lambda = seq(3, 0, length = 200),
                         standardize = F)
  # AIC
  sigma_hat_sq_3b = colSums(((y_3b-x_3b%*%beta_lasso_3b$beta)^2))/n
  AIC_3b = n*log(sigma_hat_sq_3b)+2*beta_lasso_3b$df
  # BIC
  BIC_3b = n*log(sigma_hat_sq_3b)+beta_lasso_3b$df*log(n)
  # Cross Validation
  beta_cv_3b = cv.glmnet(x_3b, y_3b, alpha = 1,
                         lambda = seq(3, 0, length = 200),
                         standardize = F,
                         nfolds = 5)
  CV_3b = beta_cv_3b$cvm
  # GCV
  GCV_3b = sigma_hat_sq_3b*(1+2*beta_lasso_3b$df/n)
  # Choose the minimum value
  min_AIC_3b = which(AIC_3b == min(AIC_3b))
  min_BIC_3b = which(BIC_3b == min(BIC_3b))
  min_CV_3b = which(CV_3b == min(CV_3b))
  min_GCV_3b = which(GCV_3b == min(GCV_3b))
  # lambda
  lambda_AIC_3b = seq(3, 0, length = 200)[min_AIC_3b]
  lambda_BIC_3b = seq(3, 0, length = 200)[min_BIC_3b]
  lambda_CV_3b = seq(3, 0, length = 200)[min_CV_3b]
  lambda_GCV_3b = seq(3, 0, length = 200)[min_GCV_3b]
  # beta of lambda
  beta_AIC_3b = beta_lasso_3b$beta[, min_AIC_3b]
  beta_BIC_3b = beta_lasso_3b$beta[, min_BIC_3b]
  beta_CV_3b = beta_lasso_3b$beta[, min_CV_3b]
  beta_GCV_3b = beta_lasso_3b$beta[, min_GCV_3b]
  # beta of lambda - beta0
  err_AIC_3b = err_AIC_3b+sqrt(sum((beta_AIC_3b-beta0)^2))
  err_BIC_3b = err_BIC_3b+sqrt(sum((beta_BIC_3b-beta0)^2))
  err_CV_3b = err_CV_3b+sqrt(sum((beta_CV_3b-beta0)^2))
  err_GCV_3b = err_GCV_3b+sqrt(sum((beta_GCV_3b-beta0)^2))
  
  # Sure and Exact
  if (sum(beta_AIC_3b[1:5]!=0) == 5){
    sure_AIC_3 = sure_AIC_3+1
  }
  if (sum(beta_AIC_3b[1:5]!=0) == 5 & sum(beta_AIC_3b[6:1000]!=0) == 0){
    exact_AIC_3 = exact_AIC_3+1
  }
  if (sum(beta_BIC_3b[1:5]!=0) == 5){
    sure_BIC_3 = sure_BIC_3+1
  }
  if (sum(beta_BIC_3b[1:5]!=0) == 5 & sum(beta_BIC_3b[6:1000]!=0) == 0){
    exact_BIC_3 = exact_BIC_3+1
  }
  if (sum(beta_CV_3b[1:5]!=0) == 5){
    sure_CV_3 = sure_CV_3+1
  }
  if (sum(beta_CV_3b[1:5]!=0) == 5 & sum(beta_CV_3b[6:1000]!=0) == 0){
    exact_CV_3 = exact_CV_3+1
  }
  if (sum(beta_GCV_3b[1:5]!=0) == 5){
    sure_GCV_3 = sure_GCV_3+1
  }
  if (sum(beta_GCV_3b[1:5]!=0) == 5 & sum(beta_GCV_3b[6:1000]!=0) == 0){
    exact_GCV_3 = exact_GCV_3+1
  }
  # Size
  size_AIC_3 = size_AIC_3+sum(beta_AIC_3b!=0)
  size_BIC_3 = size_BIC_3+sum(beta_BIC_3b!=0)
  size_CV_3 = size_CV_3+sum(beta_CV_3b!=0)
  size_GCV_3 = size_GCV_3+sum(beta_GCV_3b!=0)
}
# (1) average (beta of lambda - beta0)
err_AIC_3b = err_AIC_3b/100
err_BIC_3b = err_BIC_3b/100
err_CV_3b = err_CV_3b/100
err_GCV_3b = err_GCV_3b/100

# (2) Sure
sure_AIC_3 = sure_AIC_3/100
sure_BIC_3 = sure_BIC_3/100
sure_CV_3 = sure_CV_3/100
sure_GCV_3 = sure_GCV_3/100

# (3) Exact
exact_AIC_3 = exact_AIC_3/100
exact_BIC_3 = exact_BIC_3/100
exact_CV_3 = exact_CV_3/100
exact_GCV_3 = exact_GCV_3/100

# (4) Size
size_AIC_3 = size_AIC_3/100
size_BIC_3 = size_BIC_3/100
size_CV_3 = size_CV_3/100
size_GCV_3 = size_GCV_3/100