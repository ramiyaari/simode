

calc_nls <- function(pars, time, obs, model_out, ...) {

  if(!is.list(time)) {
    sum(unlist(lapply(names(obs),function(var)
      (model_out[,var] - obs[[var]])^2)))
  }
  else {
    time_ode <- model_out[,1]
    sum(unlist(lapply(names(obs),function(var)
      (model_out[which(time_ode %in% time[[var]]),var] - obs[[var]])^2)))
  }
}


# Calculates the negative log-likelihood for the given model output and observations,
# assuming a normal distribution for the errors
# @param pars the evaluated parameter values
# @param time the timing of the observations
# @param obs the given observations
# @param model_out the mode output returned from a call to \code{\link{solve_ode}}
# @param sigma_name name of the parameter in 'pars' that holds the value
# of the standard deviation
# @param ... additional parameters passed from a call to \code{\link{simode}}
# @return The negative log-likelihood
# @export
calc_nll_normal <- function(pars, time, obs, model_out, sigma_name, ...) {

  sigma <- pars[sigma_name]
  if(is.vector(time)){
    -sum(unlist(lapply(names(obs),function(var) {
      stats::dnorm(obs[[var]],mean=model_out[,var],sd=sigma,log=T)
    })))
  }
  else {
    time_ode <- model_out[,1]
    -sum(unlist(lapply(names(obs),function(var) {
      stats::dnorm(obs[[var]],mean=model_out[which(time_ode %in% time[[var]]),var],
                   sd=sigma,log=T)
    })))
  }
}


# Calculates the integral-matching loss for the given parameter values,
# based on the integral matching estimates for the linear parameters.
# @param pars The parameters values. Should be named
# according to the names appearing in the equations (numeric vector)
# @param equations The equations describing the ode system.
# Each equation must have a name of the variable it describes.
# An equation can contain parameters appearing in lin_pars_names, pars names,
# or variables appearing in the equations names (character vector).
# @param x0 The initial condition values. Should be the length of equations.
# Each initial condition should be named according to the variable/equation
# name it relates to. A value inside x0 could be set to NA which means
# the value of this initial condition shoule be estimated and therefore
# should appear in pars with the same name (numeric vector)
# @param time Time points of the observations. Either a vector,
# if the same time points were used for observing all variables,
# or a list of vectors the length of obs, if different time points
# were used for observing different variables (numeric list/vector).
# @param obs The observations. A list of vectors with length of equations,
# where each list member is the length of the relevant time vector (numeric list)
# @param im_method Method to use for calculating im-loss
# ("separable"/"non-separable")
# @param lin_pars_names The names of the linear parameters to be estimated
# using integral-matching (character vector)
# @param lin_pars_min Lower bounds for the linear parameter values (numeric vector)
# @param lin_pars_max Upper bounds for the linear parameter values (numeric vector)
# @param pars_min Lower bounds for the parameter values.
# Should have the same names as pars (numeric vector)
# @param pars_max Upper bounds for the parameter values.
# Should have the same names as pars (numeric vector)
# @param gen_obs An optional user defined function
# @param simode_env Environment to store simode_fit results (environment)
# in fitted smoothed curves (numeric)
# that can be used to update missing obs according to
# a model specific logic before the call to simode_fit
# @return The integral-matching loss value
#
calc_im_loss <- function(pars, equations, x0, time, obs,
                         im_method=c("separable","non-separable"),
                         lin_pars_names, lin_pars_min=NULL, lin_pars_max=NULL,
                         pars_min=NULL, pars_max=NULL, gen_obs=NULL, scale_pars=NULL,
                         pars2vars=NULL, im_smoothing=c('splines','kernel','none'),
                         im_grid_size=0, bw_factor=1.5,
                         trace=0, save_im_trace=F, simode_env=emptyenv(), ...)
{
  im_method <- match.arg(im_method)

  if(!is.null(pars_min)) {
    stopifnot(all(names(pars)==names(pars_min)))
    pars <- pmax(pars,pars_min,na.rm=T)
  }
  if(!is.null(pars_max)) {
    stopifnot(all(names(pars)==names(pars_max)))
    pars <- pmin(pars,pars_max,na.rm=T)
  }

  if(!is.null(scale_pars))
    pars <- scale_pars(pars,...)

  x0_to_est <- names(which(is.na(x0)))
  x0_opt <- intersect(x0_to_est,names(pars))
  if(!is.null(x0_opt))
    x0[x0_opt] <- pars[x0_opt]
  x0_to_est_im <- names(which(is.na(x0)))

  if(!is.null(gen_obs)) {
    timeobs <-
      gen_obs(equations=equations, pars=pars, x0=x0, time=time, obs=obs, ...)
    time <- timeobs$time
    obs <- timeobs$obs
    check_obs_validity(obs,time,names(equations))
  }

  nlin_pars <- pars[setdiff(names(pars),c(names(x0),lin_pars_names))]
  if(!pracma::isempty(nlin_pars)) {
    equations <- fix_pars(equations, nlin_pars)
  }

  vars2update <- 1:length(equations)
  if(!is.null(pars2vars)) {
    pars_prev <- simode_env$im_est_trace[nrow(simode_env$im_est_trace),names(pars)]
    if(!is.null(pars_prev)) {
      updated_pars <- which(pars!=pars_prev)
      vars2update <- pars2vars[updated_pars]
      if(length(vars2update)>0) {
        vars2update <- unique(unlist(vars2update))
      }
    }
    simode_env$pars_prev <- pars
  }

  im_fit <- NULL
  tryCatch( {
    im_fit <- simode_fit(equations=equations, pars=lin_pars_names, time=time, obs=obs, x0=x0,
                         pars_min=lin_pars_min, pars_max=lin_pars_max,
                         im_smoothing=im_smoothing, im_grid_size=im_grid_size, bw_factor=bw_factor,
                         vars2update=vars2update, im_fit_prev=simode_env$im, trace=trace)
  },
  warning = function(w) { if(trace>2) print(w) },
  error = function(e) { if(trace>1) print(e) }
  )

  if(is.null(im_fit)) {
    return (Inf)
  }

  simode_env$im <- im_fit

  lin_pars_est <- t(1)
  if(length(lin_pars_names) > 0) {
    if(im_method=='separable')
      lin_pars_est <- c(im_fit$theta)
    else
      lin_pars_est <- pars[lin_pars_names]
  }

  resid <- apply(as.matrix(1:length(equations)), 1, function(i)
    (im_fit$im_smooth[,i] - im_fit$Z[,i] - im_fit$x0[i] - im_fit$G[[i]]%*%lin_pars_est))

  im_loss <- sum(resid^2)

  if(!is.null(pars2vars) || save_im_trace) {
    im_est <- c(lin_pars_est,im_fit$x0[x0_to_est_im],pars)
    simode_env$im_est_trace <- rbind(simode_env$im_est_trace,im_est)
    simode_env$im_loss_trace <- c(simode_env$im_loss_trace,im_loss)
  }

  if(trace > 3)
    cat(noquote(paste0("im-loss: ", im_loss, "\n")))

  if(is.na(im_loss) || is.nan(im_loss))
    return(Inf)

  return (im_loss)
}


# Calculates the nonlinear least squares loss of the given parameter values,
# @param pars The parameters values. Should be named
# according to the names appearing in the equations (numeric vector)
# @param equations The equations describing the ode system.
# Each equation must have a name of the variable it describes.
# An equation can contain parameters appearing in pars
# or variables appearing in the equations names (character vector).
# @param x0 The initial condition values. Should be the length of equations.
# Each initial condition should be named according to the variable/equation
# name it relates to. A value inside x0 could be set to NA which means
# the value of this initial condition shoule be estimated and therefore
# should appear in pars with the same name (numeric vector)
# @param time Time points of the observations. Either a vector,
# if the same time points were used for observing all variables,
# or a list of vectors the length of obs, if different time points
# were used for observing different variables (numeric list/vector).
# @param obs The observations. A list of vectors with length <= the length of equations,
# where each list member is the length of the relevant time vector.
# It is possible that not all variables are observed, that is there are no
# observations for one or more of the equations (numeric list)
# @param pars_min Lower bounds for the parameter values (numeric vector)
# @param pars_max Upper bounds for the parameter values (numeric vector)
# @param calc_nll An optional user defined function that can be used
# to implement a model-specific negative-log-likelihood calculation based
# on the output of the ode model and the given observations (character)
# @param ... Additional arguments passed to \code{solve_ode}
# @return The nonlinear least squares loss value
#
calc_nls_loss <- function(pars, equations, x0, time, obs,
                          pars_min=NULL, pars_max=NULL,
                          fixed=NULL, calc_nll=calc_nls, scale_pars=NULL,
                          trace=0, save_nls_trace=F, ode_control=NULL,
                          simode_env=emptyenv(), ...)
{

  if(!is.null(pars_min)) {
    stopifnot(all(names(pars)==names(pars_min)))
    pars <- pmax(pars,pars_min,na.rm=T)
  }
  if(!is.null(pars_max)) {
    stopifnot(all(names(pars)==names(pars_max)))
    pars <- pmin(pars,pars_max,na.rm=T)
  }

  if(!is.null(scale_pars))
    pars <- scale_pars(pars,...)

  if(any(is.na(x0))) {
    x0_to_est <- names(which(is.na(x0)))
    x0[x0_to_est] <- pars[x0_to_est]
    stopifnot(any(is.na(x0))==F)
  }
  non_x0_pars <- pars[setdiff(names(pars),names(x0))]

  time_ode <- time
  if(is.list(time)) {
    time_ode <- sort(unique(unlist(time)))
    names(time) <- names(obs)
  }

  xvars_obs <- NULL
  xvars = setdiff(names(obs),names(equations))
  if(length(xvars)>0)
    xvars_obs <- obs[xvars]

  args <- c(list(equations=equations, pars=non_x0_pars,
                 x0=x0, time=time_ode, xvars=xvars_obs, trace=trace), ode_control)
  tryCatch({
    model_out <- do.call(solve_ode, args)
  },
  error = function(e) {
    print(e)
    model_out <- NULL
    })

  if(is.null(model_out) || (length(model_out[,1]) != length(time_ode))) {
    return (Inf)
  }

  eq_obs <- obs[intersect(names(obs),names(equations))]
  nls_loss <- calc_nll(pars=c(pars,fixed), time=time, obs=eq_obs, model_out=model_out, ...)

  if(is.na(nls_loss) || is.nan(nls_loss))
    return (Inf)

  if(save_nls_trace) {
    simode_env$nls_est_trace <- rbind(simode_env$nls_est_trace,pars)
    simode_env$nls_trace <- c(simode_env$nls_trace,nls_loss)
  }

  if(trace > 3)
    cat(noquote((paste0("nls_loss: ", nls_loss, "\n"))))

  return (nls_loss)

}



