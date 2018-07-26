
## generate S-system sample data ------------------------------------

## generate equations -----------------------------------------------

pars <- c('alpha1','g12','beta1','h11', 'alpha2','g21','beta2','h22')
vars <- paste0('x', 1:2)
eq1 <- 'alpha1*(x2^g12)-beta1*(x1^h11)'
eq2 <- 'alpha2*(x1^g21)-beta2*(x2^h22)'
equations <- c(eq1,eq2)
names(equations) <- vars
theta <- c(2,1,2.4,0.5,4,0.1,2,1)
names(theta) <- pars
x0 <- c(2,0.1)
names(x0) <- vars

equations
theta
x0

## generate observations --------------------------------------------

library("simode")
set.seed(1000)
n <- 50
time <- seq(0,10,length.out=n)
model_out <- solve_ode(equations,theta,x0,time)
x_det <- model_out[,vars]
sigma <- 0.05
obs <- list()
for(i in 1:length(vars)) {
  obs[[i]] <- x_det[,i] + rnorm(n,0,sigma)
}
names(obs) <- vars

## a) ODE linear in the parameters  ---------------------------------------

lin_pars <- c('alpha1','beta1','alpha2','beta2')
nlin_pars <- setdiff(pars,lin_pars)

est_lin <- simode(
  equations=equations, pars=lin_pars, fixed=c(x0,theta[nlin_pars]),
  time=time, obs=obs)

est_lin
plot(est_lin, type='fit', pars_true=theta[lin_pars], mfrow=c(1,2),legend=T)
plot(est_lin, type='est', show='both', pars_true=theta[lin_pars], legend=T)

step_size <- 0.01*est_lin$nls_pars_est
profile_lin <- profile(est_lin,step_size=step_size,max_steps=50,trace=1)
plot(profile_lin,mfrow=c(2,2))
confint(profile_lin,level=0.95)

## b) ODE semi-linear in the parameters ----------------------------------

nlin_init <-
  rnorm(length(theta[nlin_pars]),theta[nlin_pars],0.1*theta[nlin_pars])
names(nlin_init) <- nlin_pars

est_semilin <- simode(
  equations=equations, pars=pars, fixed=x0, time=time, obs=obs,
  nlin_pars=nlin_pars, start=nlin_init)

est_semilin
plot(est_semilin, type='fit', pars_true=theta, mfrow=c(1,2),legend=T)
plot(est_semilin, type='est', show='both', pars_true=theta, legend=T)

## c) ODE nonlinear in the parameters ----------------------------------

est_nosep <- simode(
  equations=equations, pars=pars, nlin_pars=pars, start=nlin_init,
  fixed=c(theta[lin_pars],x0),im_method = 'non-separable',
  time=time, obs=obs)

est_nosep
plot(est_nosep, type='est', show='both', pars_true=theta[nlin_pars], legend=T)


## d) unknown initial conditions ---------------------------------------

est_all <- simode(
  equations=equations, pars=c(pars,names(x0)), time=time, obs=obs,
  nlin_pars=nlin_pars, start=nlin_init)

summary(est_all)
plot(est_all, type='est', show='both', pars_true=c(theta,x0), legend=T)

