context("simode call arguments check")

pars <- c('alpha','beta','gamma','delta')
vars <- c('X','Y')
eq_X <- 'alpha*X-beta*X*Y'
eq_Y <- 'delta*X*Y-gamma*Y'
equations <- c(eq_X,eq_Y)
names(equations) <- vars
x0 <- c(0.9,0.9)
names(x0) <- vars
time <- 1:10
obs <- list(X=1:10,Y=1:10)


test_that("illegal equations errors are produced", {
  expect_error(simode(equations=c(1,2,3)))
  expect_error(simode(equations=c('x','y','z'),pars=c('a','b')))
})

test_that("illegal pars errors are produced", {
  expect_error(simode(equations=equations,pars=c(1,2)))
  expect_error(simode(equations=equations,pars=c(pars,'epsilon')))
})

test_that("illegal initial conditions errors are produced", {
  expect_error(simode(equations=equations,pars=pars))
  expect_error(simode(equations=equations,pars=pars,fixed=c(0.5,0.5)))
})

test_that("illegal obs/time errors are produced", {
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=c(1,2,3)))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=list(X=1,Z=3),time=1))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=list(X=1,Z=3),time=list(X=1,Y=2)))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=list(X=1,Y=2),time=c(1,2)))
})

test_that("illegal nlin_pars/likelihood_pars/start errors are produced", {
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=obs,time=time,nlin_pars='epsilon'))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=obs,time=time,nlin_pars='delta'))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=obs,time=time,start=5))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=obs,time=time,
                      start=unlist(list(epsilon=5))))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=obs,time=time,
                      likelihood_pars='delta'))
  expect_error(simode(equations=equations,pars=pars,fixed=x0,obs=obs,time=time,
                      likelihood_pars='epsilon',calc_nll=function(){return(0)}))
})

test_that("check_linearity call works", {
  expect_true(check_linearity(equations,pars))
  equations['Y'] <- 'delta*beta*X*Y-gamma*Y'
  expect_false(check_linearity(equations,pars))
  equations['Y'] <- 'delta*X*Y-Y^gamma'
  expect_false(check_linearity(equations,pars))
  equations['Y'] <- 'delta*X*Y-gamma^2*Y'
  expect_false(check_linearity(equations,pars))
})

