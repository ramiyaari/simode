context("simode basic call")

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

fit1 <- simode(equations=equations, pars=pars, fixed=x0, time=time, obs=obs)
fit2 <- simode(equations=equations, pars=c(pars,vars), time=time, obs=obs)


test_that("solve_ode returns expected format", {
  expect_is(model_out,c('deSolve','matrix'))
  expect_equal(colnames(model_out),c('time',vars))
  expect_equal(nrow(model_out),n)
  expect_equal(model_out[,1],time)
})

test_that("solve_ode2 returns expected format", {
  sol1 <- solve_ode2(fit1)
  expect_is(sol1,'list')
  expect_equal(names(sol1),c('nls','im'))
  expect_equal(sol1$nls[,1],time)
  expect_equal(sol1$im[,1],time)
  expect_equal(sol1$nls[,vars],fit1$nls_sol)
  expect_equal(sol1$im[,vars],fit1$im_sol)
})


test_that("simode returns expected format", {
  fits <- list(fit1,fit2)
  pars.list <- list(pars,c(pars,vars))
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
})

test_that("simode returns reasonable parameter estimates", {
  expect_true(all(abs((fit1$im_pars_est-theta)/theta)<0.1))
  expect_true(all(abs((fit1$nls_pars_est-theta)/theta)<0.1))
  expect_true(all(abs((fit2$im_pars_est-c(theta,x0))/c(theta,x0))<0.1))
  expect_true(all(abs((fit2$nls_pars_est-c(theta,x0))/c(theta,x0))<0.1))
})


