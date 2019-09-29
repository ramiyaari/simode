## Load multi-group SIR sample data ------------------------------

library('simode')
data(sir_example)
summary(sir_example)

equations <- sir_example$equations
S0 <- sir_example$S0
I0 <- sir_example$I0
beta <- sir_example$beta
gamma <- sir_example$gamma
kappa <- sir_example$kappa
time <- sir_example$time
I_obs <- sir_example$obs
S_obs <- lapply(1:length(S0),function(j)
  S0[j] + I0[j] - I_obs[[j]]
  - gamma*pracma::cumtrapz(time,I_obs[[j]]))
names(S_obs) <- names(S0)
obs <- c(S_obs,I_obs)

x0 <- c(S0,I0)
pars <- names(beta)
pars_min <- rep(0,length(pars))
names(pars_min) <- pars

## a) estimating beta  -------------------------------------

simode_ctrl = simode.control(trace=1, use_pars2vars_mapping=T)

est_sir_lin <- simode(
  equations=equations, pars=pars, time=time, obs=obs,
  fixed=c(gamma,kappa,x0), lower=pars_min,
  simode_ctrl = simode_ctrl)

summary(est_sir_lin)$est
x11()
plot(est_sir_lin, type='fit', which=names(I0), mfrow=c(5,2),
     time=seq(1,time[length(time)],by=0.1), pars_true=beta)
x11()
plot(est_sir_lin, type='est', show='both',pars_true=beta, legend=T)

## b) estimating beta, S0   ----------------------------------------

gen_obs <- function(equations, pars, x0, time, obs,
                    gamma, S_names, I_names, ...) {
  S0 <- x0[S_names]
  I0 <- x0[I_names]
  I_obs <- obs
  S_obs <- lapply(1:length(S0),function(i)
    S0[i]+I0[i]-I_obs[[i]]-gamma*pracma::cumtrapz(time,I_obs[[i]]))
  names(S_obs) <- S_names
  obs <- c(S_obs,I_obs)
  return (list(obs=obs, time=time))
}

pars <- c(names(beta),names(S0))
pars_min <- rep(0,length(pars))
names(pars_min) <- pars
pars_max <- rep(1,length(S0))
names(pars_max) <- names(S0)
S0_init <- rep(0.5,length(S0))
names(S0_init) <- names(S0)

est_sir_semilin <- simode(
  equations=equations, pars=pars, time=time, obs=I_obs,
  fixed=c(I0,gamma,kappa), nlin=names(S0), start=S0_init,
  lower=pars_min, upper=pars_max,
  gen_obs=gen_obs, simode_ctrl=simode_ctrl,
  gamma=gamma, S_names=names(S0), I_names=names(I0))

summary(est_sir_semilin)$est
x11()
plot(est_sir_semilin, type='est', which=names(beta), show='both',
     pars_true=beta, legend=T)
x11()
plot(est_sir_semilin, type='est', which=names(S0), show='both',
     pars_true=S0, legend=T)


## c) estimating beta, S0, gamma, kappa -------------------------

gen_obs2 <- function(equations, pars, x0, time, obs,
                     S_names, I_names, ...) {
  gamma <- pars['gamma']
  S0 <- x0[S_names]
  I0 <- x0[I_names]
  I_obs <- obs
  S_obs <- lapply(1:length(S0),function(i)
    S0[i]+I0[i]-I_obs[[i]]-gamma*pracma::cumtrapz(time,I_obs[[i]]))
  names(S_obs) <- S_names
  obs <- c(S_obs,I_obs)
  return (list(obs=obs, time=time))
}

gamma_init <- 2
names(gamma_init) <- names(gamma)
kappa_init <- rep(1,length(kappa))
names(kappa_init) <- names(kappa)

pars <- names(c(beta,gamma,kappa,S0))
nlin_pars <- names(c(gamma,kappa,S0))
start <- c(gamma_init,kappa_init,S0_init)
names(start) <- nlin_pars

pars_min <- c(rep(0,length(beta)),1.4,rep(0.25,length(kappa)),
              rep(0,length(S0)))
pars_max <- c(rep(Inf,length(beta)),3.5,rep(4,length(kappa)),
              rep(1,length(S0)))
names(pars_min) <- pars
names(pars_max) <- pars

est_sir_all <- simode(
  equations=equations, pars=pars, time=time, obs=I_obs,
  nlin_pars=nlin_pars, start=start, fixed=I0,
  lower=pars_min, upper=pars_max,
  gen_obs=gen_obs2, simode_ctrl=simode_ctrl,
  S_names=names(S0), I_names=names(I0))

summary(est_sir_all)$est
x11()
plot(est_sir_all, type='est', which=names(c(beta,gamma)), show='both',
     pars_true=c(beta,gamma), legend=T)
x11()
plot(est_sir_all, type='est', which=names(S0), show='both',
     pars_true=S0, legend=T)
x11()
plot(est_sir_all, type='est', which=names(kappa), show='both',
     pars_true=kappa, legend=T)


