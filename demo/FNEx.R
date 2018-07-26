
## generate FitzHugh-Nagumo sample data -----------------------------

## generate equations -----------------------------------------------

pars <- c('a','b','c')
vars <- c('V','R')
eq_V <- 'c*(V-V^3/3+R)'
eq_R <- '-(V-a+b*R)/c'
equations <- c(eq_V,eq_R)
names(equations) <- vars
x0 <- c(-1,1)
names(x0) <- vars
theta <- c(0.2,0.2,3)
names(theta) <- pars

## generate observations --------------------------------------------

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


## estimate model parameters with known sigma ----------------------------

calc_nll <- function(pars, time, obs, model_out, sigma, ...) {

  -sum(unlist(lapply(names(obs),function(var) {
    dnorm(obs[[var]],mean=model_out[,var],sd=sigma,log=T)
  })))
}

lin_pars <- c('a','b')
nlin_pars <- c('c')

init_vals <-
  rnorm(length(theta[nlin_pars]),theta[nlin_pars],0.1*theta[nlin_pars])
names(init_vals) <- nlin_pars

est_fn <- simode(
  equations=equations, pars=pars, time=time, obs=obs,
  fixed=x0, nlin_pars=nlin_pars, start=init_vals,
  calc_nll=calc_nll, sigma=sigma)

est_fn
plot(est_fn, type='fit', pars_true=theta, mfrow=c(1,2), legend=T)
plot(est_fn, type='est', show='both', pars_true=theta, legend=T)

## estimate model parameters and sigma ----------------------------

calc_nll_sig <- function(pars, time, obs, model_out, ...) {
  sigma <- pars['sigma']
  -sum(unlist(lapply(names(obs),function(var) {
    dnorm(obs[[var]],mean=model_out[,var],sd=sigma,log=T)
  })))
}

names(sigma) <- 'sigma'
lik_pars <- names(sigma)
pars_fn_sig <- c(pars,lik_pars)

init_vals[names(sigma)] <- 0.3
lower <- NULL
lower[names(sigma)] <- 0

est_fn_sig <- simode(
  equations=equations, pars=pars_fn_sig, time=time, obs=obs,
  fixed=x0, nlin_pars=nlin_pars, likelihood_pars=lik_pars,
  start=init_vals, lower=lower, calc_nll=calc_nll_sig)

est_fn_sig
plot(est_fn_sig, type='fit', pars_true=c(theta,sigma), mfrow=c(1,2), legend=T)
plot(est_fn_sig, type='est', show='both', pars_true=c(theta,sigma), legend=T)
