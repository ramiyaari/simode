context("simode multi call with decoupling")

M <- 4
M_grid <- expand.grid(1:M,1:M)
a_names <-
  matrix(apply(M_grid,1,function(row) paste0('a',row[1],'_',row[2])),nrow=M)
b_names <- paste0('b',1:M)
pars <- c(a_names,b_names)
vars <- paste0('x', 1:M)

equations <- c()
for(i in 1:M) {
  a_terms <- ''
  for(j in 1:M) {
    a_terms <- paste0(a_terms,a_names[i,j],'*',vars[j])
    if(j < M)
      a_terms <- paste0(a_terms,'+')
  }
  equations[i] <- paste0(b_names[i],'*',vars[i],'*(1-(',a_terms,'))')
}
names(equations) <- vars
equations

a <- matrix(c(1,1.09,1.52,0,0,1,0.44,1.36,2.33,0,1,0.47,1.21,0.51,0.35,1),4,4,byrow=T)

b <- c(1,0.72,1.53,1.27)
theta <- c(a,b)
names(theta) <- c(a_names,b_names)
x0 <- rep(0.2,M)
names(x0) <- vars

## generate observations --------------------------------------------
n <- 100
time <- seq(1,100,length=n)
model_out <- solve_ode(equations,theta,x0,time)
x_det <- model_out[,vars]

set.seed(1000)
sigma <- 0.01
obs.list <- c()
N <- 3
x0_vals <- matrix(NA,N,M)
colnames(x0_vals) <- vars
for(j in 1:N) {
  x0 <- round(rnorm(M,0.5,0.3),2)
  x0[x0<0] <- 0
  x0[x0>1] <- 1
  x0_vals[j,] <- x0
  model_out <- solve_ode(equations,theta,x0_vals[j,],time)
  x_det <- model_out[,vars]
  obs <- list()
  for(i in 1:length(vars)) {
    obs[[i]] <- rnorm(n,x_det[,i],sigma)
  }
  names(obs) <- vars
  obs.list[[j]] <- obs
}

pars_min <- rep(0,length(c(pars,vars)))
names(pars_min) <- c(pars,vars)
pars_max <- rep(1,M)
names(pars_max) <- vars

b_init <- rep(1,M)
names(b_init) <- b_names

fits <- simode(
          equations=equations, pars=c(vars,pars), nlin_pars=b_names, start=b_init,
          time=time, obs=obs.list, obs_sets=N, lower=pars_min, upper=pars_max, decouple_equations = T,
          simode_ctrl = simode.control(optim_type='im',obs_sets_fit='separate_x0'))


test_that("simode multi call with decoupling returns expected format and values", {
  expect_is(fits,'list.simode')
  expect_equal(length(fits),N)
  for(i in 1:N) {
    fit <- fits[[i]]
    obs <- obs.list[[i]]
    x0 <- x0_vals[i,]
    expect_is(fit,'simode')
    expect_equal(fit$equations,equations)
    expect_equal(fit$pars,c(vars,pars))
    expect_equal(fit$time,time)
    expect_equal(fit$obs,obs)
    expect_true(is.numeric(fit$im_loss))
    expect_true(is.numeric(fit$im_pars_est))
    expect_equal(names(fit$im_pars_est),c(vars,pars))
  }
  x0_est <- summary(fits)$im_est[,vars]
  expect_equal(dim(x0_est),dim(x0_vals))
  expect_equal(colnames(x0_est),colnames(x0_vals))
  expect_true(mean(abs((x0_est-x0_vals)))<0.05)
  a_est <- colMeans(summary(fits)$im_est[,a_names])
  expect_true(mean(abs((a_est-theta[a_names])))<0.05)
  b_est <- colMeans(summary(fits)$im_est[,b_names])
  expect_true(mean(abs((b_est-theta[b_names])))<0.2)
})



