#' Class containing control parameters for a call to \code{simode}
#'
#' @param optim_type Controls what optimization will be performed: either only
#' integral-matching ('im'), only nonlinear least squares ('nls') or both
#'  (the default, i.e., first integral-matching then nonlinear least squares starting
#'   from the integral-matching estimates).
#' @param im_optim_method Method for optimization during the integral-matching stage.
#' Accepted values are any method supported by the \code{method} argument in \code{\link{optim}},
#' as well as “Rvmmin” and “Rcgmin”, if the relevant packages are installed.
#' @param nls_optim_method Method for optimization during the nonlinear least squares stage.
#' Accepted values are the same as in \code{im_optim_method}.
#' @param im_optim_control A list with control parameters for optimization during the
#' integral-matching stage. Can include anything that would appear in the \code{control} argument
#' in \code{optim}/\code{Rvmmin}/\code{Rcgmin} (depending on the choice of \code{im_optim_method}).
#' See \code{\link{optim}}, \code{\link[Rvmmin]{Rvmmin}}, \code{\link[Rcgmin]{Rcgmin}}.
#' @param nls_optim_control Control parameters for optimization during the
#' nonlinear least squares stage (as in \code{im_optim_control})
#' @param ode_control A list with control parameters for the ODE solver. Can include the argument
#' \code{method} appearing in the arguments to \code{\link[deSolve]{ode}}, as well as any other control parameters
#' accepted as additional parameters in the call to \code{\link[deSolve]{ode}}.
#' @param im_smoothing Choice of type of smoothing during the integral-matching stage (see Details).
#' @param im_grid_size Number of points used in integral-matching grid
#' (not relevant when \code{im_smoothing='kernel'}). Value <=0 means the grid size
#' will be set according to maximum number of observations for any of the
#' equations in the call to \code{simode}.
#' @param bw_factor Controls the bandwidth when \code{im_smoothing='kernel'}.
#' The bandwidth for each equation will be bw_factor*the maximum time interval
#' between two observations (should be >= 1).
#' @param use_pars2vars_mapping Whether to use pars2vars mapping (see Details).
#' @param trace Report level (0-4), with higher values producing more tracing information (see Details).
#' @param save_im_trace Whether to save trace information of integral-matching optimization,
#' which can then be plotted using \code{\link{plot_trace}}.
#' @param save_nls_trace Whether to save trace information of nonlinear least squares optimization,
#' which can then be plotted using \code{\link{plot_trace}}.
#' @param save_to_log Controls whether to redirect output to a log file.
#' If true, output will be saved to the file 'simode.log' in tempdir.
#' @param parallel Controls whether to fit multiple \code{mc_sets} in the call to \code{simode}
#' sequentially or in parallel. Fitting in parallel requires that the parallel package will be installed.
#' When running in parallel, output will not be displayed regardless of the trace level.
#' Instead, one can set \code{save_to_log} to true to save the output to a log file.
#' @details Possible values for \code{im_smoothing} are “splines” (the default),
#' in which case smoothing will be performed using \code{\link[stats]{smooth.spline}}
#' with generalized cross-validation, “kernel”, using own kernel smoother function,
#' or “none” (using the observations as is, with interpolation if necessary).
#'
#' \code{use_pars2vars_mapping} controls whether to use a mapping of which equations
#' are affected by each of the parameters. When set to true, previous matrices computed as part of
#' the integral-matching estimation are stored during the integral-matching optimization,
#' and are updated only for the equations that were affected by the change in the
#' parameter estimates from the previous iteration.
#' When the number of equations is large and some of the parameters affect only a few equations,
#' setting this option to true can significantly reduce the optimization time during
#' the integral-matching stage (while increasing the storage usage).
#' This is especially true with derivative based optimization methods (such as “BFGS” of optim)
#' which updates only one of the optimized parameters in each iteration.
#'
#' \code{trace} has 5 possible levels:\cr
#' With trace=0, there would be no output displayed if there are no errors.\cr
#' With trace=1, a message will be displayed at the beginning and end of each optimization stage.\cr
#' With trace=2, non-critical errors occurring during the optimization iterations will be displayed.\cr
#' With trace=3, non-critical warnings occurring during the optimization iterations will be displayed.\cr
#' With trace=4, the calculated loss value for each iteration of the integral-matching and
#' nonlinear least squares optimizations will be displayed.
#' @export
simode.control <- function(optim_type=c("both","im","nls"),
                          im_optim_method=
                            c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN",
                               "Brent", "Rcgmin", "Rvmmin"),
                          nls_optim_method=
                            c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN",
                              "Brent", "Rcgmin", "Rvmmin"),
                          im_optim_control=list(), nls_optim_control=list(),
                          ode_control=list(method="lsoda"),
                          im_smoothing=c('splines','kernel','none'),
                          im_grid_size=0, bw_factor=1.5,
                          use_pars2vars_mapping=F,
                          trace=0, save_im_trace=F, save_nls_trace=F,
                          save_to_log=F, parallel=F)
{

  simode_ctrl <- list()

  simode_ctrl$optim_type <- match.arg(optim_type)

  simode_ctrl$im_optim_method <- match.arg(im_optim_method)
  simode_ctrl$nls_optim_method <- match.arg(nls_optim_method)

  simode_ctrl$im_optim_control <- im_optim_control
  simode_ctrl$nls_optim_control <- nls_optim_control

  simode_ctrl$ode_control <- ode_control

  simode_ctrl$im_smoothing <- match.arg(im_smoothing)
  simode_ctrl$im_grid_size <- im_grid_size
  simode_ctrl$bw_factor <- bw_factor

  simode_ctrl$trace <- trace

  #simode_ctrl$check_linearity <- check_linearity

  simode_ctrl$use_pars2vars_mapping <- use_pars2vars_mapping

  simode_ctrl$save_im_trace <- save_im_trace
  simode_ctrl$save_nls_trace <- save_nls_trace

  simode_ctrl$save_to_log <- save_to_log

  simode_ctrl$parallel <- parallel

  class(simode_ctrl) <- "simode.control"

  return(simode_ctrl)
}

call_optim <- function(step, args, simode_ctrl) {

  if(simode_ctrl$trace > 0) {
    if(step==1)
      msg <- "Starting optimization using integral-matching"
    else
      msg <- "Starting optimization using nonlinear least squares"
    if(!is.null(simode_ctrl$job_id))
      msg <- paste0(msg, " for job [", simode_ctrl$job_id, "]")
    msg <- paste0(msg," ...\n")
    cat(noquote(msg))
  }

  if(step==1)
    optim_method_name <- simode_ctrl$im_optim_method
  else
    optim_method_name <- simode_ctrl$nls_optim_method

  if(optim_method_name=="Rcgmin" || optim_method_name=="Rvmmin" ||
     optim_method_name=="L-BFGS-B" || optim_method_name=="Brent") {

    if(!is.null(args$pars_min)){
      args$lower <- args$pars_min
      args$pars_min <- NULL
    }
    if(!is.null(args$pars_max)){
      args$upper <- args$pars_max
      args$pars_max <- NULL
    }
  }

  if(step==2) {
    args$ode_control <- simode_ctrl$ode_control
  }

  if(optim_method_name=="Rcgmin" || optim_method_name=="Rvmmin") {
    if(!is.null(args$pars_min) && !is.null(args$pars_max)) {
      args$bdmsk <- as.numeric(args$pars_min < args$pars_max)
    }
    if(optim_method_name=="Rcgmin")
      optim_method <- Rcgmin::Rcgmin
    else
      optim_method <- Rvmmin::Rvmmin
  } else {
    args$method <- optim_method_name
    optim_method <- stats::optim
  }

  optim_result <- NULL
  err_msg <- NULL
  tryCatch(
    {
      optim_runtime <- system.time( {
        optim_result <-  do.call(optim_method, args)
      })
    },
    # warning = function(w) { }
    error = function(e) { err_msg <<- conditionMessage(e) }
  )

  if(is.null(optim_result)) {
    stop(err_msg)
  }

  if(!is.null(args$pars_min))
    optim_result$par <- pmax(optim_result$par,args$pars_min,na.rm=T)
  if(!is.null(args$pars_max))
    optim_result$par <- pmin(optim_result$par,args$pars_max,na.rm=T)

  if(simode_ctrl$trace > 0) {
    if(step==1)
      msg <- "Integral-matching optimization"
    else
      msg <- "Nonlinear least squares optimization"
    if(!is.null(simode_ctrl$job_id))
      msg <- paste0(msg, " for job [", simode_ctrl$job_id, "]")
    msg <- paste0(msg," completed in [",
                  round(optim_runtime['elapsed'],2), "] sec ",
                  "with convergence code [", optim_result$convergence, "]\n")
    cat(noquote(msg))
  }

  return(optim_result)
}


simode_linear <- function(simode_obj, simode_env, ...) {

  equations <- simode_obj$equations
  pars <- simode_obj$pars
  x0 <- simode_obj$x0
  time <- simode_obj$time
  obs <- simode_obj$obs
  likelihood_pars <- simode_obj$likelihood_pars
  start <- simode_obj$start
  gen_obs <- simode_obj$gen_obs
  calc_nll <- simode_obj$calc_nll
  if(is.null(calc_nll))
    calc_nll <- calc_nls
  ctrl <- simode_obj$ctrl

  do_im_opt <- (ctrl$optim_type == 'im' || ctrl$optim_type=='both')
  do_nls_opt <- (ctrl$optim_type == 'nls' || ctrl$optim_type=='both')

  if(do_im_opt) {

    lin_pars <- setdiff(pars,c(names(x0),likelihood_pars))
    lin_pars_min <- simode_obj$lower[lin_pars]
    lin_pars_max <- simode_obj$upper[lin_pars]

    #if(ctrl$check_linearity) {
    if(!check_linearity(equations, lin_pars))
      return (NULL)
    #}

    # -----------------------------------------------------------------------------
    # estimating linear parameters using im ---------------------------------------

    im_loss <- calc_im_loss(
      pars=NULL, equations=equations, x0=x0, time=time, obs=obs,
      lin_pars_names=lin_pars,lin_pars_min=lin_pars_min,
      lin_pars_max=lin_pars_max, gen_obs=gen_obs,
      im_smoothing=ctrl$im_smoothing, im_grid_size=ctrl$im_grid_size,
      bw_factor=ctrl$bw_factor, trace=ctrl$trace+1, save_im_trace=F,
      simode_env=simode_env, ...)

    if(is.infinite(im_loss) || is.null(im_loss)) {
      cat('Error during integral-matching estimation\n')
      return (NULL)
    }

    if(ctrl$trace > 0)
      cat('Integral-matching estimation done\n')

    likelihood_pars_im_est <- NULL
    if(!is.null(likelihood_pars)) {
      likelihood_pars_im_est <- rep(NA,length(likelihood_pars))
      names(likelihood_pars_im_est) <- likelihood_pars
    }

    pars_est <- c(simode_env$im$theta ,simode_env$im$x0[which(is.na(x0))])
    if(!is.null(likelihood_pars)) {
      pars_est <- c(pars_est,start[likelihood_pars])
    }
    pars_est <- pars_est[pars]

    simode_obj$im_pars_est <- pars_est
    simode_obj$im_pars_est[likelihood_pars] <- NA
    simode_obj$im_loss <- im_loss
    simode_obj$im_smooth <- list(time=simode_env$im$t,val=simode_env$im$im_smooth)

  } else { #do_im_opt==F
    pars_est <- start[pars]
  }

  if(do_nls_opt) {
  # -----------------------------------------------------------------------------
  # running optim using full likelihood -----------------------------------------

    optim_args <-
      list(par=pars_est, fn=calc_nls_loss,
              control=ctrl$nls_optim_control,
              equations=equations, x0=x0, time=time, obs=obs,
              pars_min=simode_obj$lower, pars_max=simode_obj$upper,
              fixed=simode_obj$fixed,
              calc_nll=calc_nll, trace=ctrl$trace,
              save_nls_trace=ctrl$save_nls_trace,
              simode_env=simode_env, ...)

    optim_result <- call_optim(step=2, args=optim_args, ctrl)

    simode_obj$nls_loss <- optim_result$value
    simode_obj$nls_pars_est <- optim_result$par

    if(ctrl$save_nls_trace) {
      simode_obj$nls_trace <- simode_env$nls_trace
      simode_obj$nls_est_trace <- simode_env$nls_est_trace[,pars]
    }
  }
  return (simode_obj)
}


simode_nonlinear <- function(simode_obj, simode_env, ...) {

  equations <- simode_obj$equations
  pars <- simode_obj$pars
  x0 <- simode_obj$x0
  time <- simode_obj$time
  obs <- simode_obj$obs
  im_method <- simode_obj$im_method
  nlin_pars <- simode_obj$nlin_pars
  likelihood_pars <- simode_obj$likelihood_pars
  start <- simode_obj$start
  pars_min <- simode_obj$lower
  pars_max <- simode_obj$upper
  gen_obs <- simode_obj$gen_obs
  calc_nll <- simode_obj$calc_nll
  if(is.null(calc_nll))
    calc_nll <- calc_nls
  ctrl <- simode_obj$ctrl

  do_im_opt <- (ctrl$optim_type == 'im' || ctrl$optim_type=='both')
  do_nls_opt <- (ctrl$optim_type == 'nls' || ctrl$optim_type=='both')

  if(do_im_opt) {

    nlin_pars_init <- start[nlin_pars]
    nlin_pars_min <- pars_min[nlin_pars]
    nlin_pars_max <- pars_max[nlin_pars]

    im_pars <- setdiff(pars,likelihood_pars)
    lin_pars <- setdiff(im_pars,c(nlin_pars,names(x0)))
    lin_pars_min <- pars_min[lin_pars]
    lin_pars_max <- pars_max[lin_pars]

    #if(ctrl$check_linearity) {
    if(!pracma::isempty(lin_pars) && !check_linearity(equations, lin_pars))
      return (NULL)
    #}

    if(im_method=='separable') {
      opt_pars_est <- nlin_pars_init
      opt_pars_min <- nlin_pars_min
      opt_pars_max <- nlin_pars_max
    }
    else {
      lin_pars_init <- start[intersect(lin_pars,names(start))]
      missing_init_lin_pars <- setdiff(lin_pars,names(lin_pars_init))
      x0_init <- start[intersect(names(x0),names(start))]
      missing_init_x0 <- setdiff(names(x0),names(x0_init))
      if(!pracma::isempty(missing_init_lin_pars) ||
         !pracma::isempty(missing_init_x0)) {
        im_loss <- calc_im_loss(
          pars=nlin_pars_init, equations=equations, x0=x0, time=time, obs=obs,
          lin_pars_names=lin_pars, lin_pars_min=lin_pars_min, lin_pars_max=lin_pars_max,
          gen_obs=gen_obs, im_smoothing=ctrl$im_smoothing,
          im_grid_size=ctrl$im_grid_size, bw_factor=ctrl$bw_factor,
          trace=ctrl$trace+1, simode_env=simode_env, ...)

        if(is.infinite(im_loss) || is.null(im_loss)) {
          cat('Error during initial integral-matching estimation\n')
          return (NULL)
        }

        lin_pars_init[missing_init_lin_pars] <-
                                    simode_env$im$theta[missing_init_lin_pars]
        x0_init[missing_init_x0] <- simode_env$im$x0[missing_init_x0]
      }
      opt_pars_est <- c(lin_pars_init, nlin_pars_init, x0_init)[im_pars]
      opt_pars_min <- pars_min[im_pars]
      opt_pars_max <- pars_max[im_pars]
    }

    # -----------------------------------------------------------------------------
    # running optim using im-loss -------------------------------------------------

    pars2vars <- NULL
    if(ctrl$use_pars2vars_mapping) {
      pars2vars <- pars2vars_mapping(names(opt_pars_est), equations)
    }

    optim_args <-list(par=opt_pars_est, fn=calc_im_loss,
                   control=ctrl$im_optim_control,
                   equations=equations, x0=x0, time=time, obs=obs,
                   im_method=im_method, lin_pars_names=lin_pars,
                   lin_pars_min=lin_pars_min, lin_pars_max=lin_pars_max,
                   pars_min=opt_pars_min, pars_max=opt_pars_max,
                   gen_obs=gen_obs, pars2vars=pars2vars,
                   im_smoothing=ctrl$im_smoothing,
                   im_grid_size=ctrl$im_grid_size,
                   bw_factor=ctrl$bw_factor,
                   trace=ctrl$trace,
                   save_im_trace=ctrl$save_im_trace,
                   simode_env=simode_env, ...)

    optim1_result <- call_optim(step=1, args=optim_args, ctrl)

    im_loss <- optim1_result$value
    opt_pars_est <- optim1_result$par

    if(length(lin_pars)>0 && im_method=="separable") {
      im_loss2 <- calc_im_loss(pars=opt_pars_est, equations=equations, x0=x0,
                     time=time, obs=obs, im_method=im_method, lin_pars_names=lin_pars,
                     lin_pars_min=lin_pars_min, lin_pars_max=lin_pars_max,
                     gen_obs=gen_obs,
                     im_smoothing=ctrl$im_smoothing,
                     im_grid_size=ctrl$im_grid_size,
                     bw_factor=ctrl$bw_factor,
                     simode_env=simode_env, ...)
      pars_est <-
        c(opt_pars_est, simode_env$im$theta ,simode_env$im$x0[which(is.na(x0))])
    }
    else {
      pars_est <- c(opt_pars_est)
    }

    if(!is.null(likelihood_pars)) {
      pars_est <- c(pars_est,start[likelihood_pars])
    }
    pars_est <- pars_est[pars]

    simode_obj$im_pars_est <- pars_est
    simode_obj$im_pars_est[likelihood_pars] <- NA
    simode_obj$im_loss <- im_loss
    simode_obj$im_smooth <- list(time=simode_env$im$t,val=simode_env$im$im_smooth)

    if(ctrl$save_im_trace) {
      simode_obj$im_loss_trace <- simode_env$im_loss_trace
      simode_obj$im_est_trace <- simode_env$im_est_trace[,im_pars]
    }
  }
  else { #do_im_opt==F
    pars_est <- simode_obj$start[pars]
  }

  if(do_nls_opt) {

  # -----------------------------------------------------------------------------
  # running optim using full likelihood ------------------------------------

    optim_args <-
      list(par=pars_est, fn=calc_nls_loss,
           control=ctrl$nls_optim_control,
           equations=equations, x0=x0, time=time, obs=obs,
           pars_min=pars_min, pars_max=pars_max,
           fixed=simode_obj$fixed, calc_nll=calc_nll,
           trace=ctrl$trace, save_nls_trace=ctrl$save_nls_trace,
           simode_env=simode_env, ...)

    optim2_result <- call_optim(step=2, args=optim_args, ctrl)

    simode_obj$nls_loss <- optim2_result$value
    simode_obj$nls_pars_est <- optim2_result$par

    if(ctrl$save_nls_trace) {
      simode_obj$nls_trace <- simode_env$nls_trace
      simode_obj$nls_est_trace <- simode_env$nls_est_trace[,pars]
    }
  }

  return (simode_obj)
}


simode_create <-
  function(call, equations, pars, time, obs,
           nlin_pars, likelihood_pars,
           fixed, start, lower, upper,
           im_method, gen_obs, calc_nll, simode_ctrl) {

  stopifnot(is.character(equations), is.character(pars))

  stopifnot(!is.null(names(equations)))

  model_pars <- setdiff(pars,likelihood_pars)
  if(!check_pars_presence_in_equations(equations,model_pars))
    stop('Remove unused parameters from \'pars\'')

  x0 <- NULL
  if(!is.null(fixed)) {
    fixed_names <- names(fixed)
    stopifnot(is.numeric(fixed),!is.null(fixed_names))
    fixed_x0 <- intersect(fixed_names,names(equations))
    x0[fixed_x0] <- fixed[fixed_x0]
    fixed_pars <- fixed[setdiff(fixed_names,c(fixed_x0,likelihood_pars))]
    if(length(fixed_pars)>0) {
      equations <- fix_pars(equations, fixed_pars)
    }
    pars <- setdiff1(pars, fixed_names)
    nlin_pars <- setdiff1(nlin_pars, fixed_names)
    # likelihood_pars <- setdiff1(likelihood_pars, fixed_names)
    start <- setdiff2(start, fixed_names)
    lower <- setdiff2(lower, fixed_names)
    upper <- setdiff2(upper, fixed_names)
  }
  unknown_x0 <- setdiff(names(equations),names(x0))
  if(!pracma::isempty(setdiff(unknown_x0,pars))) {
    missing <- paste('[',setdiff(unknown_x0,pars),']',collapse ='')
    stop(paste0('unknown x0 must appear in pars - missing: ', missing))
  }
  x0[unknown_x0] <- NA
  x0 <- x0[names(equations)]

  stopifnot(is.list(obs),length(obs)<=length(equations),!is.null(names(obs)))
  if(length(obs)<length(equations) && is.null(gen_obs) && simode_ctrl$optim_type!='nls')
    stop(paste0('cannot perform integral-matching optimization - ",
                "missing observations for variables and no gen_obs method supplied'))
  if(!pracma::isempty(setdiff(names(obs),names(equations)))) {
    extra <- paste('[',setdiff(names(obs),names(equations)),']',collapse='')
    stop(paste0('observations contain variables that do not appear in equations:',
                extra))
  }
  stopifnot(is.numeric(time) || (is.list(time) && length(time)==length(obs)))
  for(i in 1:length(obs)) {
    if(is.numeric(time))
      stopifnot(length(time)==length(obs[[i]]))
    else
      stopifnot(is.numeric(time[[i]]),length(time[[i]])==length(obs[[i]]))
  }

  stopifnot(is.null(gen_obs) || is.function(gen_obs))
  stopifnot(is.null(calc_nll) || is.function(calc_nll))

  if(!is.null(nlin_pars))
    stopifnot(is.character(nlin_pars))
  if(!is.null(start))
    stopifnot(is.numeric(start),!is.null(names(start)))

  if(!pracma::isempty(setdiff(names(start),pars))) {
    extra <- paste('[',setdiff(names(start),pars),']',collapse='')
    stop(paste0("start contain values for parameters/initial conditions ",
                "that do not appear in pars: ",extra))
  }
  if(!pracma::isempty(setdiff(nlin_pars,pars))) {
    extra <- paste('[',setdiff(nlin_pars,pars),']',collapse='')
    stop(paste0("nlin_pars contain parameters that do not appear in pars: ", extra))
  }
  if(!pracma::isempty(setdiff(nlin_pars,names(start)))) {
    missing <- paste('[',setdiff(nlin_pars,names(start)),']',collapse='')
    stop(paste0("start must contain initial values for all parameters ",
                "in nlin_pars - missing: ",missing))
  }

  if(!is.null(likelihood_pars)) {
    stopifnot(is.character(likelihood_pars))
    if(is.null(calc_nll)) {
      stop('likelihood_pars should be empty if calc_nll function is not defined',immediate.=T)
    }
    if(!pracma::isempty(setdiff(likelihood_pars,pars))) {
      extra <- paste('[',setdiff(likelihood_pars,pars),']',collapse='')
      stop(paste0("likelihood_pars contain parameters that do not appear in pars: ", extra))
    }
    if(!pracma::isempty(setdiff(likelihood_pars,names(start)))) {
      missing <- paste('[',setdiff(likelihood_pars,names(start)),']',collapse='')
      stop(paste0("start must contain initial values for all parameters ",
                  "in likelihood_pars - missing: ",missing))
    }
  }

  if(simode_ctrl$optim_type!='nls') {
    if(im_method=='separable') {
      if(!pracma::isempty(setdiff(names(start),c(nlin_pars,likelihood_pars)))) {
        warning(paste0("start contain values for parameters that do not ",
                       "appear in nlin_pars or likelihood_pars - with im_method=\'separable\' ",
                       "these values will be ignored."), immediate.=T)
      }
      start <- start[c(nlin_pars,likelihood_pars)]
    }
  } else { #simode_ctrl$optim_type=='nls'
    if(!pracma::isempty(setdiff(pars,names(start)))) {
      missing <- paste('[',setdiff(pars,names(start)),']',collapse='')
      stop(paste0("when running simode with optim_type==\'nls\', start ",
                  "must contain initial values for all parameters - missing:",missing))
    }
  }

  if(!is.null(lower)) {
    stopifnot(is.numeric(lower),!is.null(names(lower)))
    # stopifnot(length(lower)<=length(pars))
    if(!pracma::isempty(setdiff(names(lower),pars))) {
      extra <- paste('[',setdiff(names(lower),pars),']',collapse='')
      warning(paste0("lower contain parameter names that do not",
                    " appear in pars and will be ignored: ",extra), immediate.=T)
    }
    lower_all <- rep(-Inf,length(pars))
    names(lower_all) <- pars
    lower_names <- intersect(pars,names(lower))
    lower_all[lower_names] <- lower[lower_names]
    lower <- lower_all
  }
  if(!is.null(upper)) {
    stopifnot(is.numeric(upper),!is.null(names(upper)))
    # stopifnot(length(upper)<=length(pars))
    if(!pracma::isempty(setdiff(names(upper),pars))) {
      extra <- paste('[',setdiff(names(upper),pars),']',collapse='')
      warning(paste0("upper contain parameter names that do not",
                   " appear in pars and will be ignored: ", extra), immediate.=T)
    }
    upper_all <- rep(Inf,length(pars))
    names(upper_all) <- pars
    upper_names <- intersect(pars,names(upper))
    upper_all[upper_names] <- upper[upper_names]
    upper <- upper_all
  }

  if(!(simode_ctrl$im_optim_method %in%
       c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
         "Brent", "Rcgmin", "Rvmmin"))) {
    warning("unknown im_optim_method selected - using default \'BFGS\'",
            immediate.=T)
    simode_ctrl$im_optim_method <- "BFGS"
  }
  if(!(simode_ctrl$nls_optim_method %in%
       c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
         "Brent", "Rcgmin", "Rvmmin"))) {
    warning("unknown nls_optim_method selected - using default \'BFGS\'",
            immediate.=T)
    simode_ctrl$nls_optim_method <- "BFGS"
  }
  if((simode_ctrl$im_optim_method=="Rvmmin" ||
      simode_ctrl$nls_optim_method=="Rvmmin") &&
     !requireNamespace("Rvmmin", quietly=T)) {
    warning("Rvmmin package not installed - using default \'BFGS\'",immediate.=T)
    if(simode_ctrl$im_optim_method=="Rvmmin")
      simode_ctrl$im_optim_method <- "BFGS"
    if(simode_ctrl$nls_optim_method=="Rvmmin")
      simode_ctrl$nls_optim_method <- "BFGS"
  }
  if((simode_ctrl$im_optim_method=="Rcgmin" ||
      simode_ctrl$nls_optim_method=="Rcgmin") &&
     !requireNamespace("Rcgmin", quietly=T)) {
    warning("Rcgmin package not installed - using default \'BFGS\'",immediate.=T)
    if(simode_ctrl$im_optim_method=="Rcgmin")
      simode_ctrl$im_optim_method <- "BFGS"
    if(simode_ctrl$nls_optim_method=="Rcgmin")
      simode_ctrl$nls_optim_method <- "BFGS"
  }

  if(!(simode_ctrl$im_smoothing %in% c("splines", "kernel", "none"))) {
    warning("unknown im_smoothing selected - using default \'splines\'",
            immediate.=T)
    simode_ctrl$im_smoothing <- "splines"
  }

  simode <- list(call=call, equations=equations, pars=pars, x0=x0,
                time=time, obs=obs, im_method=im_method, nlin_pars=nlin_pars,
                likelihood_pars=likelihood_pars, fixed=fixed,
                start=start, lower=lower, upper=upper,
                gen_obs=gen_obs, calc_nll=calc_nll, ctrl=simode_ctrl)
  class(simode) <- "simode"

  return (simode)
}


#' Statistical inference of ordinary differential equations using
#' separable integral-matching
#'
#' Estimating the parameters of an ODE system in two stages:
#' 1) Estimate the parameters using separable integral-matching,
#' 2) Estimate the parameters using nonlinear least squares
#' starting from the values obtained in stage 1.
#' @param equations Named vector. The equations describing the ODE system.
#' Each element of the vector should contain a character representation
#' of the right-hand side of an equation, and should be named according to
#' the left-hand side of the equation (i.e., the variable name).
#' An equation can contain parameters appearing in \code{pars},
#' variables appearing in the equations names and/or any function of 't',
#' which is a reserved symbol for the time domain.
#' @param pars The names of the parameters and initial conditions
#' to be estimated. An initial condition name for a certain variable
#' is the name given to the relevant equation in \code{equations}
#' (e.g., if an equation is named 'x' than its initial condition should be named 'x'
#' as well). Note: The symbol 't' is reserved for the time domain and cannot be
#' used as a parameter name.
#' @param time Time points of the observations. Either a vector,
#' if the same time points were used for observing all variables,
#' or a list of vectors the length of \code{obs}, of which each element
#' is the length of the relevant element in \code{obs}.
#' @param obs Named list. The observations. When \code{mc_sets=1}, \code{obs}
#' should contain a list of vectors with length <= equations (either partial or fully observed
#' system). Each list member should be named according to the relevant equation name and should
#' be the length of the relevant time vector. When \code{mc_sets>1}, \code{obs} should contain
#' a list the length of \code{mc_sets}, where each list member is a list that fits the
#' description in the case of \code{mc_sets=1}.
#' @param mc_sets Number of monte-carlo expriments. When \code{mc_sets>1}, the function will
#' fit each set of observations separately either sequentially or in parallel, according
#' to the value of \code{parallel} in \code{simode_ctrl}.
#' @param nlin_pars Names of parameters or initial conditions
#' that will be estimated in stage 1 using nonlinear least squares optimization.
#' The parameter names in \code{nlin_pars} must appear in \code{pars}.
#' @param likelihood_pars Names of likelihood parameters not appearing in the ODE system,
#' which are needed for the the user-defined function \code{calc_nll}.
#' The parameter names in \code{likelihood_pars} must appear in \code{pars}.
#' @param fixed Named vector. Fixed values for one or more of the ODE system parameters or
#' initial conditions. Parameters in this list will not be estimated.
#' @param start Named vector. Starting values for optimization of parameters/initial conditions.
#' Must contain starting values for all the parameters in \code{nlin_pars}
#' and \code{likelihood_pars}. If im_method="non-seperable", can optionally contain
#' starting values for any other parameter/initial condition.
#' @param lower Named vector. Lower bounds for any parameter/initial condition.
#' @param upper Named vector. Upper bounds for any parameter/initial condition.
#' @param im_method The method to use for integral-matching. Default "separable"
#' means that linear parameters are estimated directly while "non-separable"
#' means that linear parameters are estimated using nonlinear least squares optimization.
#' If none of the parameters are linear then the default can be used.
#' @param gen_obs A user-defined function for completing missing observations (see Details).
#' @param calc_nll A user-defined function for calculating negative log-likelihood for
#' the model (see Details).
#' @param simode_ctrl Various control parameters. See \code{\link{simode.control}}.
#' @param ... Additional arguments passed to \code{optim}, \code{gen_obs} and \code{calc_nll}
#' @details \code{gen_obs} can be used in cases of a partially observed system,
#' for which observations of the missing variables can be generated given values for the
#' system parameters. The function will be called during the optimization using integral-matching.\cr
#' It must be defined as
#' \code{gen_obs <- function(equations, pars, x0, time, obs, ...)}, where:
#' \itemize{
#' \item \code{equations} the ODE equations
#' \item \code{pars} the parameter values
#' \item \code{x0} the initial conditions
#' \item \code{time} the timing of the observations (vector or list)
#' \item \code{obs} the observations
#' \item \code{...} additional parameters passed from the call to \code{\link{simode}}
#' }
#' The function should return a list with two items:
#' \itemize{
#' \item \code{time} the vector or list of time points of the observations
#' \item \code{obs} the list of observations with the newly generated observations
#' }
#'
#' \code{calc_nll} allows the user to pass his own likelihood function to be used in the
#' optimization in the second stage (if not defined, the default nonlinear least squares
#' optimization will be used). The likelihood function will also be used in a following
#' call to \code{\link{profile}}, for the calculation of likelihood profiles.
#' It must be defined as
#' \code{calc_nll <- function(pars, time, obs, model_out, ...)}, where:
#' \itemize{
#' \item \code{pars} the parameter values
#' \item \code{time} the timing of the observations (vector or list)
#' \item \code{obs} the observations
#' \item \code{model_out} the model output returned from a call to \code{\link{solve_ode}}.
#' If \code{time} is a list with possibly different times for each variable
#' then \code{model_out} will contain a union of all these times.
#' \item \code{...} additional parameters passed from the call to \code{\link{simode}}
#' }
#' The function should return the negative log-likelihood.
#'
#' @return If mc_sets=1, returns a \code{simode} object containing the
#' parameter estimates after integral-matching (stage 1) and after
#' nonlinear least squares optimization (stage 2). If \code{mc_sets>1} returns a
#' \code{list.simode} object which is a list of \code{simode} objects the length of
#' \code{mc_sets}.
#' @references
#' Dattner & Klaassen (2015). Optimal Rate of Direct Estimators in Systems of Ordinary Differential Equations Linear in Functions of the Parameters,
#' Electronic Journal of Statistics, Vol. 9, No. 2, 1939-1973.
#'
#' Dattner, Miller, Petrenko, Kadouriz, Jurkevitch & Huppert (2017). Modelling and Parameter Inference of Predator-prey Dynamics in Heterogeneous Environments Using The Direct Integral Approach,
#' Journal of The Royal Society Interface 14.126: 20160525.
#' @examples
#'
#' ## =================================================
#' ## Predator-Prey Lotka-Volterra model
#' ## =================================================
#'
#' ## generate model equations and parameters (X=Prey,Y=Predator)
#' pars <- c('alpha','beta','gamma','delta')
#' vars <- c('X','Y')
#' eq_X <- 'alpha*X-beta*X*Y'
#' eq_Y <- 'delta*X*Y-gamma*Y'
#' equations <- c(eq_X,eq_Y)
#' names(equations) <- vars
#' x0 <- c(0.9,0.9)
#' names(x0) <- vars
#' theta <- c(2/3,4/3,1,1)
#' names(theta) <- pars
#'
#' ## generate observations
#' n <- 50
#' time <- seq(0,25,length.out=n)
#' model_out <- solve_ode(equations,theta,x0,time)
#' x_det <- model_out[,vars]
#' set.seed(1000)
#' sigma <- 0.05
#' obs <- list()
#' for(i in 1:length(vars)) {
#'   obs[[i]] <- pmax(0, rnorm(n,x_det[,i],sigma))
#' }
#' names(obs) <- vars
#'
#' ## estimate model parameters with known initial conditions
#' simode_fit1 <- simode(equations=equations, pars=pars, fixed=x0, time=time, obs=obs)
#' plot(simode_fit1, type='fit', time=seq(0,25,length.out=100), pars_true=theta, mfrow=c(2,1))
#' plot(simode_fit1, type='est', pars_true=theta)
#'
#' \donttest{
#' ## estimate model parameters and initial conditions
#' simode_fit2 <- simode(equations=equations, pars=c(pars,vars), time=time, obs=obs)
#' plot(simode_fit2, type='fit', time=seq(0,25,length.out=100), pars_true=c(theta,x0), mfrow=c(2,1))
#' plot(simode_fit2, type='est', pars_true=c(theta,x0))
#'
#' profiles_fit2 <- profile(simode_fit2,step_size=0.01,max_steps=50)
#' plot(profiles_fit2,mfrow=c(2,3))
#' ci_fit2 <- confint(profiles_fit2)
#' ci_fit2
#' plot(ci_fit2,pars_true=c(theta,x0),legend=T)
#' }
#'
#' @importFrom pracma interp1 trapz cumtrapz lsqlincon
#' @importFrom stats D smooth.spline predict
#' @export
#'
simode <- function(equations, pars, time, obs, mc_sets=1,
                  nlin_pars=NULL, likelihood_pars=NULL,
                  fixed=NULL, start=NULL, lower=NULL, upper=NULL,
                  im_method=c("separable","non-separable"),
                  gen_obs=NULL, calc_nll=NULL, simode_ctrl=simode.control(), ...) {

  if(is.list(equations)){
    eq_names <- names(equations)
    equations <- c(equations)
    names(equations) <- eq_names
  }

  im_method <- match.arg(im_method)

  if(is.null(simode_ctrl$job_id))
    on.exit(expr=simode_cleanup(simode_ctrl$save_to_log))

  stopifnot(mc_sets>0)
  if(mc_sets > 1) {
    stopifnot(mc_sets==length(obs))
    return (simode_multi(equations=equations, pars=pars,  time=time, obs=obs,
                        mc_sets=mc_sets, nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
                        fixed=fixed, start=start, lower=lower, upper=upper, im_method=im_method,
                        gen_obs=gen_obs, calc_nll=calc_nll, simode_ctrl=simode_ctrl, ...))
  }

  simode_obj <- simode_create(
              match.call(), equations, pars, time, obs,
              nlin_pars, likelihood_pars,
              fixed, start, lower, upper,
              im_method, gen_obs, calc_nll, simode_ctrl)

  simode_obj$extra_args <- list(...)

  if(simode_obj$ctrl$save_to_log) {
    log_file <- file.path(tempdir(),"simode.log")
    sink(file=log_file, append=F)
    cat(noquote(paste0("Call to simode on [", Sys.time(), "]:\n")))
  }

  return (simode_impl(simode_obj))
}

simode_impl <- function(x) {

  simode_env <- new.env(parent=emptyenv())
  simode_obj <- NULL
  tryCatch(
    {
      args <- c(list(simode_obj=x,simode_env=simode_env),x$extra_args)
      if(x$im_method=='separable' && pracma::isempty(x$nlin_pars))
        simode_obj <- do.call("simode_linear", args)
      else
        simode_obj <- do.call("simode_nonlinear", args)
    },
    error = function(e) { print(e) },
    finally = {
      rm(simode_env)
    }
  )
  return (simode_obj)
}

simode_cleanup <- function(save_to_log)
{
  old.options <- options(warn=-1)
  log_file <- file.path(tempdir(),"simode.log")
  if(save_to_log && file.exists(log_file)) {
    sink(NULL)
    cat(paste0('Log file available at: ',log_file))
  }
  tmpfile <- file.path(tempdir(),"simode-ode.log")
  if(file.exists(tmpfile)) {
    sink(NULL)
    file.remove(tmpfile)
  }
  options(old.options)
}


#' Plot the fit/estimates of a \code{simode} object
#'
#' Plot the fit or parameter estimates obtained from a call to \code{simode}.
#'
#' @param x \code{simode} object returned by a call to \code{\link{simode}}.
#' @param type Type of plot - 'fit' to plot the fitted variables
#' and 'est' to plot the parameter estimates.
#' @param show Whether to plot the fit/estimates obtained
#' using nonlinear least squares ('nls'), integral-matching ('im') or both ('both').
#' @param which Which variables to plot in case \code{type='fit'}, or which
#' parameters to plot in case \code{type='est'}. If empty, the plot will include
#' all of the variables/parameters in x.
#' @param pars_true The true parameter values (if are known). Should be
#' named using the parameter names. If given, the true values for
#' the variables/parameters will be added to the plot.
#' @param time The time points to use for the fitted curves
#' (relevant only for \code{type='fit'}).
#' If not given then the time points of the observations in x will be used.
#' @param plot_im_smooth Whether or not to plot the smoothed curves created and
#' used by the integral-matching procedure (relevant only for \code{type='fit'}).
#' @param legend Whether or not to add a figure legend.
#' @param mfrow A vector of the form c(nr,nc) setting the layout of
#' subplots in one plot (see also \code{\link{par}}).
#' @param cols List of colors for each element of the plot.
#' @param ... Additional argument(s) for methods.
#' @importFrom graphics plot points lines axis par
#' @importFrom pracma errorbar
#' @export
#'
plot.simode <- function(x, type=c('fit','est'), show=c('nls','im','both'),
                       which=NULL, pars_true=NULL, time=NULL,
                       plot_im_smooth=F, legend=F, mfrow=par('mfrow'),
                       cols=list(nls_fit="blue",im_fit="green", true="black",
                                 obs="red", im_smooth="magenta"), ...) {

  x_list <- list(x)
  class(x_list) <- 'list.simode'
  plot(x_list,type=type,show=show,which=which,pars_true=pars_true,
       time=time,plot_im_smooth=plot_im_smooth,legend=legend,
       mfrow=mfrow,cols=cols,...)
}

#' Plot optimization trace of a call to \code{simode}
#'
#' Plot a trace of the loss values and parameter estimates during
#' the integral-matching/nonlinear least squares optimization within a call to \code{simode}.
#' For the traces to exist, the arguments \code{save_im_trace} and/or
#' \code{save_nls_trace} in \code{\link{simode.control}} should be set to true,
#' when calling \code{simode}.
#' @param x \code{simode} object returned by a call to \code{\link{simode}}.
#' @param show Whether to plot the estimates obtained
#' using nonlinear least squares ('nls'), integral-matching ('im') or both ('both').
#' @param which Which parameters' traces to plot. If NULL, the
#' trace for all the parameters in \code{x} will be plotted.
#' @param mfrow A vector of the form c(nr,nc) setting the layout of
#' subplots in one plot (see also \code{\link{par}}).
#' @param cols List of colors for each element of the plot.
#' @param ... Additional argument(s) for methods.
#' @export
#'
plot_trace <- function(x, show=c('nls','im','both'),
                       which=NULL, mfrow=par('mfrow'),
                       cols=list(nls_fit="blue",im_fit="green"), ...) {

  show <- match.arg(show)

  if(show=="im" || show=="both") {
    if(is.null(x$im_loss_trace) || is.null(x$im_est_trace))
      stop(paste0('missing im trace in given simode object - ',
                  'run simode again with \'save_im_trace\' set to true'))
  }
  if(show=="nls" || show=="both") {
    if(is.null(x$nls_trace) || is.null(x$nls_est_trace))
      stop(paste0('no nls trace in given simode object - ',
                  'run simode again with \'save_nls_trace\' set to true'))
  }

  stopifnot(!is.null(x$pars))
  pars <- which
  if(is.null(pars))
    pars <- x$pars

  old.par <- par(mar=c(4,4,2,2),mfrow=mfrow)
  on.exit(par(old.par))

  if(show=="im" || show=="both") {
    gof_trace <- x$im_loss_trace
    ylim <- c(min(gof_trace,na.rm=T),min(max(gof_trace,na.rm=T),1.05*gof_trace[1]))
    plot(gof_trace,type='p',pch=20,ylab='im_loss',
         ylim=ylim, col=cols[['im_fit']], ...)
  }
  if(show=="nls" || show=="both") {
    gof_trace <- x$nls_trace
    ylim <- c(min(gof_trace,na.rm=T),min(max(gof_trace,na.rm=T),1.05*gof_trace[1]))
    plot(gof_trace,type='p',pch=20,ylab='nls-loss',
         ylim=ylim, col=cols[['nls_fit']], ...)
  }
  for(p in pars) {
    if(show=="im")
      plot(x$im_est_trace[,p],type='p',pch=20,ylab=p,col=cols[['im_fit']],...)
    else if(show=="nls")
      plot(x$nls_est_trace[,p],type='p',pch=20,ylab=p,col=cols[['nls_fit']],...)
    else {
      xlim <- c(0,length(x$im_est_trace[,p])+length(x$nls_est_trace[,p]))
      ylim <- c(min(c(x$im_est_trace[,p],x$nls_est_trace[,p])),
                max(c(x$im_est_trace[,p],x$nls_est_trace[,p])))
      plot(x$im_est_trace[,p],type='p',pch=20,ylab=p,
           xlim=xlim,ylim=ylim,col=cols[['im_fit']],...)
      points(seq(length(x$im_est_trace[,p])+1,
                 length(x$im_est_trace[,p])+length(x$nls_est_trace[,p])),
              x$nls_est_trace[,p],pch=20,col=cols[['nls_fit']],...)
    }
  }
}

#' Print method for \code{simode} objects
#'
#' @param x The \code{simode} object.
#' @param ... Additional argument(s) for methods.
#' @export
print.simode <- function(x, ...) {
  print(summary(x, ...), ...)
}


#' Print method for \code{summary.simode} objects
#'
#' @param x The \code{summary.simode} object.
#' @param ... Additional argument(s) for methods.
#' @export
print.summary.simode <- function(x, ...) {

  cat("call:\n")
  print(x$call)
  cat("\nequations:\n")
  print(x$equations)
  cat("\ninitial conditions:\n")
  print(x$x0)
  cat("\nparameter estimates:\n")
  print(x$est)

  cat("\nim-method: ", x$im_method,"\n")
  if(!is.null(x$im_loss)) {
    cat("\nim-loss: ", x$im_loss,"\n")
  }
  if(!is.null(x$nls_loss)) {
    cat("\nnls-loss: ", x$nls_loss,"\n")
  }
}


#' Summary method for \code{simode} objects
#'
#' @param object The \code{simode} object.
#' @param digits The number of significant digits to use.
#' @param ... Additional argument(s) for methods.
#' @export
summary.simode <- function(object, digits=max(3, getOption("digits")-3), ...) {

  summary <- list()
  summary$call <- object$call
  summary$equations <- object$equations
  summary$x0 <- object$x0

  summary$im_method <- object$im_method
  summary$im_loss <- signif(object$im_loss,digits)
  summary$nls_loss <- signif(object$nls_loss,digits)

  df <- data.frame(par=object$pars)
  df$type <- rep('linear',length(object$pars))
  df$type[which(object$pars %in% object$nlin_pars)] <- 'non-linear'
  df$type[which(object$pars %in% object$likelihood_pars)] <- 'likelihood'

  if(!is.null(object$lower)){
    lower_all <- rep(-Inf,length(object$pars))
    names(lower_all) <- object$pars
    lower_names <- intersect(object$pars,names(object$lower))
    lower_all[lower_names] <- object$lower[lower_names]
    df$lower <- lower_all
  }
  if(!is.null(object$upper)) {
    upper_all <- rep(Inf,length(object$pars))
    names(upper_all) <- object$pars
    upper_names <- intersect(object$pars,names(object$upper))
    upper_all[upper_names] <- object$upper[upper_names]
    df$upper <- upper_all
  }
  if(!is.null(object$start)) {
    start_all <- rep(NA,length(object$pars))
    names(start_all) <- object$pars
    start_names <- intersect(object$pars,names(object$start))
    start_all[start_names] <- object$start[start_names]
    df$start <- start_all
  }
  if(!is.null(object$im_pars_est))
    df$im_est <- signif(object$im_pars_est,digits)
  if(!is.null(object$nls_pars_est))
    df$nls_est <- signif(object$nls_pars_est,digits)

  summary$est <- df

  class(summary) <- "summary.simode"
  return(summary)
}

#' Example dataset for a multi-group SIR model
#'
#' A model for the spread of a seasonal influenza epidemics in two groups
#' over five seasons.
#' @format A list containing the following variables:
#' \describe{
#'   \item{equations}{The ODE equations describing the system, including
#'   ten equations for the susceptible dynamics (2 groups * 5 seasons) and
#'   ten equations for the infected dynamics.}
#'   \item{beta}{The 2x2 transmission matrix parameters (in time unit of weeks).}
#'   \item{gamma}{The recovery rate parameter (in time unit of weeks).}
#'   \item{kappa}{The relative infectiousness of seasons 2-5 compared to season 1.}
#'   \item{S0}{The initial conditions for the susceptible variables.}
#'   \item{I0}{The initial conditions for the infected variables.}
#'   \item{time}{Times in which the observations were made (in weeks).}
#'   \item{obs}{A list of observations of the infected variables,
#'   generated using Gaussian measurement error with sigma=1e-3.}
#' }
"sir_example"


