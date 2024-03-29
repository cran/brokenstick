## ----setup, include=FALSE-----------------------------------------------------
old <- options(digits = 3)
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library("plyr")
library("mvtnorm")
make_data_generator <- function(resid_var = 1,
                                ranef_covar = diag(c(1, 1)), n = 100
                                ) {
  ni <- nrow(ranef_covar)
  generate_data <- function() {
    # sample data set under mixed effects model with random slope/intercepts 
    simulated_data <- rdply(n, {
      b <- t(rmvnorm(n = 1, sigma = ranef_covar))
      epsilon <- rnorm(n = length(b), mean = 0, sd = sqrt(resid_var))
      b + epsilon
    })
  data.frame(
    subject = rep(1:n, each = ni),
    age = rep(1:ni, n),
    simulated_data)
  }
}

## -----------------------------------------------------------------------------
resid_var <- 0
resid_var <- 1
set.seed(77711)
covar <- matrix(c(1, 0.7, 0.5, 0.3,
                  0.7, 1, 0.8, 0.5,
                  0.5, 0.8, 1, 0.6,
                  0.3, 0.5, 0.6, 1), nrow = 4)
gen_dat <- make_data_generator(n = 10000, 
                               ranef_covar = covar,
                               resid_var = resid_var)
data <- gen_dat()
head(data)

## -----------------------------------------------------------------------------
library("tidyr")
library("dplyr")
d <- as_tibble(data[,-3])
broad <- t(spread(d, subject, X1))[-1,]
cor(broad)

## -----------------------------------------------------------------------------
library("brokenstick")
knots <- 1:3
boundary <- c(1, 4)
fit <- brokenstick(X1 ~ age | subject, data, 
                   knots = knots, boundary = boundary,
                   method = "lmer")
omega <- get_omega(fit, hide = "no")
beta <- coef(fit, hide = "no")
sigma2 <- fit$sigma2
round(beta, 2)
round(sigma2, 4)

# correlation random effects
round(covar, 3)
round(omega, 2)

# covariances measured data
round(omega + diag(sigma2, 4), 3)
round(cov(broad), 3)

# convert to time-to-time correlation matrix
round(cov2cor(omega + diag(sigma2, 4)), 3)
round(cor(broad), 3)

z <- predict(fit, x = "knots", include_data = FALSE, shape = "wide")[, -1]
# off-diagonal elements of covariance of broken stick estimates approach correlation
# not enough variance in the diagonal because of smoothing
cov(z)
# correlations of broken stick estimates are inflated because of smoothing
cor(z)

## -----------------------------------------------------------------------------
cov <- get_omega(fit)
chat <- cov + diag(fit$sigma2, nrow(cov))
r <- cov2cor(chat)
r

## ----echo = FALSE-------------------------------------------------------------
options(old)

