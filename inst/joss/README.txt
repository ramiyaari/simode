simode v1.1.5:
An R package for conducting statistical inference for ordinary differential equations that aims to ease the optimization process 
and provide more robust solutions to parameter estimation problems. 

Getting Started:
From an R console, type install.packages('simode') to install the package.
After download is complete type library('simode') to load the package.

Example:
Following is an example of usage of the package using the Lotka-Volterra model.

## =================================================
## Predator-Prey Lotka-Volterra model
## =================================================

## generate model equations and parameters (X=Prey,Y=Predator)
pars <- c('alpha','beta','gamma','delta')
vars <- c('X','Y')
eq_X <- 'alpha*X-beta*X*Y'
eq_Y <- 'delta*X*Y-gamma*Y'
equations <- c(eq_X,eq_Y)
names(equations) <- vars
x0 <- c(0.9,0.9)
names(x0) <- vars
theta <- c(2/3,4/3,1,1)
names(theta) <- pars

## generate observations
n <- 50
time <- seq(0,25,length.out=n)
model_out <- solve_ode(equations,theta,x0,time)
x_det <- model_out[,vars]
set.seed(1000)
sigma <- 0.05
obs <- list()
for(i in 1:length(vars)) {
  obs[[i]] <- pmax(0, rnorm(n,x_det[,i],sigma))
}
names(obs) <- vars

## estimate model parameters with known initial conditions
simode_fit1 <- simode(equations=equations, pars=pars, fixed=x0, time=time, obs=obs)
plot(simode_fit1, type='fit', time=seq(0,25,length.out=100), pars_true=theta, mfrow=c(2,1))
plot(simode_fit1, type='est', pars_true=theta)

## estimate model parameters and initial conditions
simode_fit2 <- simode(equations=equations, pars=c(pars,vars), time=time, obs=obs)
plot(simode_fit2, type='fit', time=seq(0,25,length.out=100), pars_true=c(theta,x0), mfrow=c(2,1))
plot(simode_fit2, type='est', pars_true=c(theta,x0))

profiles_fit2 <- profile(simode_fit2,step_size=0.01,max_steps=50)
plot(profiles_fit2,mfrow=c(2,3))
ci_fit2 <- confint(profiles_fit2)
ci_fit2
plot(ci_fit2,pars_true=c(theta,x0),legend=T)

## ================================================##

More examples can be found as demos inside the package (type: demo(package='simode'))

Documentation:
From an R console type ?simode for help on the package.
The package page at CRAN is at: https://cran.r-project.org/web/packages/simode/index.html
Package reference manual is at: https://cran.r-project.org/web/packages/simode/simode.pdf
Package vignette is at: https://cran.r-project.org/web/packages/simode/vignettes/R_package_simode.pdf

Support:
For support, bug reports or contributions send an email to: ramiyaari@gmail.com  

