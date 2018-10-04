
## generate Lotka-Volterra with seasonal-forcing sample data --------

## generate equations -----------------------------------------------

# X=Pred,Y=Prey

pars <- c('alpha','beta','gamma','delta','epsilon')
vars <- c('X','Y')
eq_X <- 'alpha*X-beta*(1+epsilon*seasonality)*X*Y'
eq_Y <- 'delta*(1+epsilon*seasonality)*X*Y-gamma*Y'
equations <- c(eq_X,eq_Y)
names(equations) <- vars
x0 <- c(0.9,0.9)
names(x0) <- vars
theta <- c(2/3,4/3,1,1,0.2)
names(theta) <- pars

## generate observations -----------------------------------------------

library("simode")
n <- 100
time <- seq(0,50,length.out=n)
seasonality <- rep(c(0,0,0,0,0,1,1,1,1,1),n/10)
xvars <- list(seasonality=seasonality)
model_out <- solve_ode(equations,theta,x0,time,xvars=xvars)
x_det <- model_out[,vars]
set.seed(1000)
sigma <- 0.1
obs <- list()
for(i in 1:length(vars)) {
  obs[[i]] <- x_det[,i] + rnorm(n,0,sigma)
}
names(obs) <- vars

# plot(time, x_det[,1], type='l', ylab=vars[1])
# points(time, obs[[1]], col='red')
# plot(time, x_det[,2], type='l', ylab=vars[2])
# points(time, obs[[2]], col='red')

# fit data  ----------------------------------------------------------

nlin_pars <- c('epsilon')
nlin_init <- c(0.25)
names(nlin_init) <- nlin_pars

simode_fit <- simode(
  equations=equations, pars=pars, fixed=x0,
  time=time, obs=c(obs,xvars),
  nlin_pars=nlin_pars, start=nlin_init,
  simode_ctrl = simode.control(nls_optim_method='Nelder-Mead'))

simode_fit
plot(simode_fit, type='fit', pars_true=theta, legend=T, mfrow=c(2,1))
plot(simode_fit, show='both', type='est', pars_true=theta, legend=T)

