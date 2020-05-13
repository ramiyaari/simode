context("simode call with user-defined likelihood function")

pars1 <- c('a','b','c')
vars <- c('V','R')
eq_V <- 'c*(V-V^3/3+R)'
eq_R <- '-(V-a+b*R)/c'
equations <- c(eq_V,eq_R)
names(equations) <- vars
x0 <- c(-1,1)
names(x0) <- vars
theta <- c(0.2,0.2,3)
names(theta) <- pars1

n <- 40
time <- seq(0,20,length.out=n)
model_out <- solve_ode(equations,theta,x0,time)
x_det <- model_out[,vars]

set.seed(1000)
sigma <- 0.05
obs <- list()
for(i in 1:length(vars)) {
  obs[[i]] <- x_det[,i] + rnorm(n,0,sigma)
}
names(obs) <- vars

nlin_pars <- c('c')
init_vals <-
  rnorm(length(theta[nlin_pars]),theta[nlin_pars],0.1*theta[nlin_pars])
names(init_vals) <- nlin_pars

calc_nll <- function(pars, time, obs, model_out, sigma, ...) {
  -sum(unlist(lapply(names(obs),function(var) {
      dnorm(obs[[var]],mean=model_out[,var],sd=sigma,log=TRUE)
      })))
}

fit1 <- simode(
          equations=equations, pars=pars1, time=time, obs=obs, fixed=x0,
          nlin_pars=nlin_pars, start=init_vals, calc_nll=calc_nll, sigma=sigma)


lik_pars <- 'sigma'
pars2 <- c(pars1,lik_pars)
init_vals[lik_pars] <- 0.1
lower <- NULL
lower[lik_pars] <- 0
calc_nll2 <- function(pars, time, obs, model_out, ...) {
  sigma <- pars['sigma']
  -sum(unlist(lapply(names(obs),function(var) {
    dnorm(obs[[var]],mean=model_out[,var],sd=sigma,log=TRUE)
  })))
}

fit2 <- simode(
  equations=equations, pars=pars2, time=time, obs=obs, fixed=x0,
  nlin_pars=nlin_pars, likelihood_pars=lik_pars, start=init_vals, lower=lower, calc_nll=calc_nll2)


test_that("simode call with user-defined lik", {
  fits <- list(fit1,fit2)
  pars.list <- list(pars1,pars2)
  for(i in 1:length(fits)) {
    fit <- fits[[i]]
    pars <- pars.list[[i]]
    expect_is(fit,'simode')
    expect_equal(fit$equations,equations)
    expect_equal(fit$pars,pars)
    expect_equal(fit$time,time)
    expect_equal(fit$obs,obs)
    expect_true(is.numeric(fit$im_loss))
    expect_true(is.numeric(fit$nls_loss))
    expect_equal(dim(fit$im_sol),c(n,2))
    expect_equal(colnames(fit$im_sol),vars)
    expect_true(is.numeric(fit$im_sol))
    expect_equal(dim(fit$nls_sol),c(n,2))
    expect_equal(colnames(fit$nls_sol),vars)
    expect_true(is.numeric(fit$nls_sol))
    expect_true(is.numeric(fit$im_pars_est))
    expect_true(is.numeric(fit$nls_pars_est))
    expect_equal(names(fit$im_pars_est),pars)
    expect_equal(names(fit$nls_pars_est),pars)
  }
  expect_true(all(abs((fit1$im_pars_est-theta)/theta)<0.2))
  expect_true(all(abs((fit1$nls_pars_est-theta)/theta)<0.2))
  expect_true(is.na(fit2$im_pars_est[lik_pars]))
  expect_true(all(abs((fit2$im_pars_est[pars1]-theta)/theta)<0.2))
  expect_true(all(abs((fit2$nls_pars_est-c(theta,sigma))/c(theta,sigma))<0.2))
})

