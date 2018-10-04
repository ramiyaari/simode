
## generate Lotka-Volterra with closed-form seasonal-forcing sample data --

## generate equations -----------------------------------------------------

# X=Pred,Y=Prey

pars <- c('alpha','beta','gamma','delta','epsilon','omega')
vars <- c('X','Y')
eq_X <- 'alpha*X-beta*(1+epsilon*sin(2*pi*(t/50+omega)))*X*Y'
eq_Y <- 'delta*(1+epsilon*sin(2*pi*(t/50+omega)))*X*Y-gamma*Y'
equations <- c(eq_X,eq_Y)
names(equations) <- vars
x0 <- c(0.9,0.9)
names(x0) <- vars
theta <- c(2/3,4/3,1,1,0.2,0.5)
names(theta) <- pars

## generate observations -----------------------------------------------

library("simode")
n <- 100
time <- seq(0,50,length.out=n)
model_out <- solve_ode(equations,theta,x0,time)
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

nlin_pars <- c('epsilon','omega')
nlin_init <- c(0.3,0.4)
names(nlin_init) <- nlin_pars

pars_min <- c(0,0)
names(pars_min) <- nlin_pars
pars_max <- c(1,1)
names(pars_max) <- nlin_pars

## fit using simplex ##
simode_fit <- simode(
  equations=equations, pars=pars, fixed=c(x0), time=time, obs=obs,
  nlin_pars=nlin_pars, start=nlin_init, lower=pars_min, upper=pars_max,
  simode_ctrl=simode.control(im_optim_method='Nelder-Mead',
                             nls_optim_method='Nelder-Mead'))

simode_fit
plot(simode_fit, type='fit', pars_true=theta, mfrow=c(2,1), legend=T)
plot(simode_fit, show='both', type='est', pars_true=theta, legend=T)

