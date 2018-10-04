
## generate Lotka-Volterra with monte-carlo simulations sample data ------

## generate equations ----------------------------------------------------

# X=Pred,Y=Prey

pars <- c('alpha','beta','gamma','delta')
vars <- c('X','Y')
eq_X <- 'alpha*X-beta*X*Y'
eq_Y <- 'delta*X*Y-gamma*Y'
equations <- c(eq_X,eq_Y)
names(equations) <- vars
x0 <- c(0.9,0.9)
names(x0) <- vars
theta <- c(2/3,4/3,2,1)
names(theta) <- pars

# generate monte-carlo observations  -------------------------------------

library("simode")
n <- 100
time <- seq(0,25,length.out=n)
model_out <- solve_ode(equations,theta,x0,time)
x_det <- model_out[,vars]
N <- 5
set.seed(1000)
sigma <- 0.1
mc_obs <- list()
for(j in 1:N) {
  obs <- list()
  for(i in 1:length(vars)) {
    obs[[i]] <- rnorm(n,x_det[,i],sigma)
  }
  names(obs) <- vars
  mc_obs[[j]] <- obs
}

# fit monte-carlo observations  ------------------------------------------

simode_fits <- simode(
  equations=equations, pars=c(pars,vars), time=time, obs=mc_obs,
  obs_sets=N, simode_ctrl=simode.control(trace=1,parallel=T,save_to_log=T))

summary(simode_fits,sum_mean_sd=T,pars_true=c(theta[pars],x0),digits=2)
plot(simode_fits, type='fit', plot_mean_sd=T, pars_true=c(theta,x0), mfrow=c(2,1), legend=T)
plot(simode_fits, type='est', show='both', plot_mean_sd=T,pars_true=c(theta,x0), legend=T)

