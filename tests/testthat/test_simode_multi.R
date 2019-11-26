context("multi simode call check")

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
n <- 50
time <- seq(0,25,length.out=n)
model_out <- solve_ode(equations,theta,x0,time)
x_det <- model_out[,vars]
N <- 3
set.seed(1000)
sigma <- 0.05
mc_obs <- list()
for(j in 1:N) {
  obs <- list()
  for(i in 1:length(vars)) {
    obs[[i]] <- rnorm(n,x_det[,i],sigma)
  }
  names(obs) <- vars
  mc_obs[[j]] <- obs
}

fits <- simode(equations=equations, pars=pars, fixed=x0, time=time, obs=mc_obs, obs_sets=N)

test_that("simode multi call returns expected format and values", {
  expect_is(fits,'list.simode')
  expect_equal(length(fits),N)
  for(i in 1:N) {
    fit <- fits[[i]]
    obs <- mc_obs[[i]]
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
    expect_true(all(abs((fit$im_pars_est-theta)/theta)<0.1))
    expect_true(all(abs((fit$nls_pars_est-theta)/theta)<0.1))
  }
})



