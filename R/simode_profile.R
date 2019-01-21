
#' Calculate likelihood profiles for the model parameters
#'
#' @param fitted \code{simode} object returned by a call to \code{\link{simode}}.
#' @param which Which parameters to estimate the profile for.
#' @param optim_type Whether to calculate the profiles based on
#' maximum-likelihood optimization only ('nls') or based on integral-matching
#' followed by maximum-likelihood optimization ('both').
#' @param step_size Step size for profiling (one value for all parameters or
#' a value for each parameter in \code{which}).
#' @param max_steps Maximum number of steps to take in each direction.
#' @param alpha Maximum (two-sided) likelihood ratio test confidence level to find.
#' @param skip_err Whether on not to stop the calculation if encountering a problem
#' with one point in the profile.
#' @param trace Report level (0-4), with higher values producing more tracing information.
#' @param save_to_log Whether to redirect output to log file. The log file
#' will be saved to tempdir().
#' @param ...	 Additional argument(s) for methods.
#' @details If the call to \code{\link{simode}}, which returned the fitted object given to
#' this method, included a user-defined likelihood function (with the \code{calc_nll} argument),
#' then the likelihood profiles will be calculated using this function. Otherwise,
#' the profiles will be calculated using a likelihood based on a Gaussian distribuion
#' with fixed sigma, where sigma will be estimated in the background together with
#' the rest of the model parameters.
#' @return The likelihood profiles.
#' @importFrom stats qchisq
#' @export
#'
profile.simode <- function(
  fitted, which=NULL, optim_type=c("nls","both"),
  step_size, max_steps=100, alpha=0.05, skip_err=T, trace=0, save_to_log=F, ...) {

  optim_type <- match.arg(optim_type)

  if(is.null(fitted$nls_pars_est)) {
    stop(paste0('missing nls parameter estimates in given simode object - ',
                'run simode again with \'optim_type\' set to \'nls\' or \'both\''))
  }

  stopifnot(is.null(which) || is.character(which))
  if(is.null(which))
    which <- fitted$pars

  if(!pracma::isempty(setdiff(which,fitted$pars))) {
    extra <- paste('[',setdiff(which,fitted$pars),']',collapse='')
    stop(paste0("\'which\' contain parameters that do not appear in fitted$pars: ",
                extra))
  }
  stopifnot(is.numeric(alpha),alpha>0,alpha<1)
  stopifnot(is.numeric(step_size))
  if(length(step_size)==1)
    step_size <- rep(step_size,length(which))
  stopifnot(length(step_size)==length(which))
  stopifnot(is.numeric(max_steps),max_steps>=1)

  max_delta <- qchisq(1-alpha,1)/2

  if(is.null(fitted$calc_nll)) {
    x0.na <- names(fitted$x0[is.na(fitted$x0)])
    x0_nls <- fitted$x0
    x0_nls[x0.na] <- fitted$nls_pars_est[x0.na]
    pars_est_nls <- fitted$nls_pars_est[setdiff(names(fitted$nls_pars_est),x0.na)]
    model_out <- solve_ode(fitted$equations, pars_est_nls, x0_nls, fitted$time)
    vars <- names(fitted$equations)
    model_all_vars <- unlist(lapply(vars, function(var) model_out[,var]))
    obs_all_vars <- unlist(lapply(vars, function(var) fitted$obs[[var]]))
    sigma <- stats::sd(model_all_vars-obs_all_vars)
    sigma_name <- 'lik_normal_sigma'

    # fitted$fixed[sigma_name] <- sigma

    fitted$start[sigma_name] <- sigma
    if(is.null(fitted$lower)) {
      fitted$lower[fitted$pars] <- -Inf
    }
    fitted$lower[sigma_name] <- 0
    if(is.null(fitted$upper)) {
      fitted$upper[fitted$pars] <- Inf
    }
    fitted$upper[sigma_name] <- sigma*10
    fitted$nls_pars_est[sigma_name] <- sigma
    fitted$im_pars_est[sigma_name] <- NA
    fitted$pars <- c(fitted$pars,sigma_name)

    fitted$likelihood_pars <- c(fitted$likelihood_pars,sigma_name)
    fitted$calc_nll <- calc_nll_normal
    fitted$extra_args <- c(fitted$extra_args,sigma_name=sigma_name)
  }

  if(save_to_log) {
    log_file <- file.path(tempdir(),"simode.log")
    sink(file=log_file, append=T)
    cat(noquote((paste0("Call to profile on [", Sys.time(), "]:\n"))))
  }

  on.exit(expr=simode_cleanup(save_to_log))

  res <- list(call=match.call(), alpha=alpha, pars=which,
              nls_pars_est=fitted$nls_pars_est[which],
              lower=fitted$lower[which],
              upper=fitted$upper[which])
  class(res) <- 'profile.simode'

  res$profiles <- list()
  tryCatch({
    for(i in 1:length(which)) {
      par <- which[i]
      res$profiles[[par]] <- calc_profile(fitted, par, optim_type, step_size[i],
                                          max_steps, max_delta, skip_err, trace)
    }
  },
  warning=function(w) {
    print(w)
  },
  error=function(e) {
    print(e)
    res <- NULL
  }
  )
  return (res)
}

calc_profile <- function(x, par, optim_type, step_size, max_steps,
                         max_delta, skip_err, trace) {

  if(trace > 0)
    cat(noquote((paste0("Calculating profile for parameter [", par, "]\n"))))

  is_lik_par <- F
  if(!pracma::isempty(intersect(x$likelihood_pars,par)))
    is_lik_par <- T

  x$ctrl$optim_type <- optim_type
  x$pars <- setdiff(x$pars,par)
  x$nlin_pars <- setdiff(x$nlin_pars,par)
  x$likelihood_pars <- setdiff(x$likelihood_pars,par)
  par_min <- unname(x$lower[par])
  par_max <- unname(x$upper[par])
  x$lower <- x$lower[x$pars]
  x$upper <- x$upper[x$pars]
  x$start <- x$start[setdiff(names(x$start),par)]
  x$ctrl$job_id <- NULL
  x$ctrl$trace <- trace-1

  nls_par_est <- x$nls_pars_est[par]
  par_vals <- c(unname(nls_par_est))
  init_pars_est_l <- x$nls_pars_est[x$pars]
  init_pars_est_r <- x$nls_pars_est[x$pars]
  step_res <- calc_profile_step(x, par, nls_par_est, is_lik_par,
                                init_pars_est_l, skip_err, trace)
  min_nll <- step_res$nll
  nll_vals <- min_nll

  profile_l_done <- F
  profile_r_done <- F
  i <- 0
  while(T)
  {
    i <- i+1
    if(trace > 1)
      cat(noquote((paste0("Performing step [", i, "] in calculating profile for parameter [", par ,"]\n"))))

    if(!profile_l_done) {
      par_val <- nls_par_est-i*step_size
      if(!is.null(par_min) && par_val <= par_min) {
        par_val <- par_min
        profile_l_done <- T
      }
      step_res <-
        calc_profile_step(x, par, par_val, is_lik_par,
                          init_pars_est_l, skip_err, trace)
      if(!is.null(step_res)) {
        par_vals <- c(unname(par_val),par_vals)
        nll <- step_res$nll
        nll_vals <- c(nll,nll_vals)
        if(min_nll > nll)
          min_nll <- nll
        if(nll - min_nll > max_delta)
          profile_l_done <- T
        init_pars_est_l <- step_res$nls_pars_est
        if(trace > 1)
          cat(noquote((paste0("Parameter [", par, "], Step [", i, "] L, Delta NLL [",
                              round(nll-min_nll,2) ,"]\n"))))

      }
    }
    if(!profile_r_done) {
      par_val <- nls_par_est+i*step_size
      if(!is.null(par_max) && par_val >= par_max) {
        par_val <- par_max
        profile_r_done <- T
      }
      step_res <-
        calc_profile_step(x, par, par_val, is_lik_par,
                          init_pars_est_r, skip_err, trace)
      if(!is.null(step_res)) {
        par_vals <- c(par_vals,unname(par_val))
        nll <- step_res$nll
        nll_vals <- c(nll_vals,nll)
        if(min_nll > nll)
          min_nll <- nll
        if(nll - min_nll > max_delta)
          profile_r_done <- T
        init_pars_est_r <- step_res$nls_pars_est
        if(trace > 1)
          cat(noquote((paste0("Parameter [", par, "], Step [", i, "] R, Delta NLL [",
                              round(nll-min_nll,2) ,"]\n"))))
      }
    }
    if(profile_l_done && profile_r_done)
      break
    if(i >= max_steps) {
      if(trace > 0) {
        cat(noquote((paste0("Calculating profile for parameter [", par,
                            "] - max step has been reached\n"))))
      }
      break
    }
  }

  return (list(par_vals=par_vals, nll_vals=nll_vals))
}

calc_profile_step <- function(x, par, par_val, is_lik_par,
                              init_pars_est, skip_err, trace) {

  if(is_lik_par)
    x$fixed<- c(x$fixed,par_val)
  else if(par %in% names(x$x0))
    x$x0[par] <- par_val
  else
    x$equations <- fix_pars(x$equations, par_val)

  x$start  <- init_pars_est

  if(x$ctrl$optim_type=='both' && x$decouple_equations)
    x <- simode_decouple(x)
  else
    x <- simode_impl(x)

  if(is.null(x)) {
    if(!skip_err)
      stop(paste0('Error occurred during profile calculation for parameter [', par, ']\n'))
    return (NULL)
  }

  res <- list(nll=x$nls_loss, nls_pars_est=x$nls_pars_est)
  return (res)
}

#' Plot the likelihood profiles for the model parameters
#'
#' @param x \code{profile.simode} object returned by a call to \code{\link{profile}}.
#' @param which Which parameters to plot the likelihood profiles for.
#' If empty, the plot will include all of the parameters in x.
#' @param mfrow A vector of the form c(nr,nc) setting the layout of
#' subplots in one plot (see also \code{\link{par}}).
#' @param cols List of colors for each element of the plot.
#' @param ... Additional argument(s) for methods.
#' @importFrom graphics abline
#' @export
#'
plot.profile.simode <- function(x, which=NULL, mfrow=par('mfrow'),
                               cols=list(fit="blue",threshold="red"), ...) {

  max_delta <- qchisq(1-x$alpha,1)/2

  if(is.null(which))
    which <- x$pars
  if(!pracma::isempty(setdiff(which,x$pars))) {
    extra <- paste('[',setdiff(which,x$pars),']',collapse='')
    stop(paste0("\'which\' contain parameters that do not appear in x: ",extra))
  }

  old.par <- par(mar=c(4,4,2,2),mfrow=mfrow)
  on.exit(par(old.par))

  for(i in 1:length(which)) {
    par <- which[i]
    plot(x$profiles[[par]]$par_vals,
         x$profiles[[par]]$nll_vals-min(x$profiles[[par]]$nll_vals,na.rm=T),
         xlab=par,ylab=expression(paste(Delta,' NLL')),
         col=cols[['fit']], type='o', ...)
    abline(h=max_delta, col=cols[['threshold']])
  }
}

#' Print method for \code{profile.simode} objects
#'
#' @param x The \code{profile.simode} object.
#' @param ... Additional argument(s)for methods.
#' @export
print.profile.simode <- function(x, ...) {
  cat("call:\n")
  print(x$call)
  cat("alpha:\n")
  cat(x$alpha,"\n")
  cat("profiles:\n")
  print(x$profiles)
}

#' Calculates confidence intervals for the model parameters
#'
#' Calculates confidence intervals for the model parameters,
#' based on the given likelihood profiles
#' @param object A fitted model object: \code{profile.simode} object returned by a call to \code{\link{profile}}.
#' @param parm A specification of which parameters are to be given confidence intervals (named vector).
#' If missing, all parameters are considered.
#' @param level The confidence level required.
#' @param ...	 Additional argument(s) for methods.
#' @return The confidence intervals.
#' @importFrom stats approx
#' @export
#'
confint.profile.simode<- function(object, parm=NULL, level=0.95, ...) {

  stopifnot(is.null(parm) || (is.character(parm) && is.vector(parm)))
  pars <- parm
  if(is.null(pars))
    pars <- c(object$pars)
  if(!pracma::isempty(setdiff(pars,object$pars))) {
    extra <- paste('[',setdiff(pars,object$pars),']',collapse='')
    stop(paste0("\'parm\' contain parameters that do not appear in object$pars: ",
                extra))
  }

  stopifnot(is.numeric(level),level>0,level<1)
  if(1-level<object$alpha)
    stop('1-level must be greater than value of alpha used in profile call')

  max_delta <- qchisq(level,1)/2
  res <- list(call=match.call(), level=level)
  intervals <- list()
  for(par in pars) {
    par_vals <- object$profiles[[par]]$par_vals
    nll_vals <- object$profiles[[par]]$nll_vals
    par_min <- object$lower[[par]]
    par_max <- object$upper[[par]]
    intervals[[par]] <- calc_confint(par_vals,nll_vals,max_delta,par_min,par_max)
  }
  ci_df <- data.frame(par=pars)
  if(!is.null(object$nls_pars_est))
    ci_df$nls_est <- object$nls_pars_est

  lp <- length(pars)
  ci_df$lower <- apply(as.matrix(1:lp),1,function(i) intervals[[i]]$lower)
  ci_df$upper <- apply(as.matrix(1:lp),1,function(i) intervals[[i]]$upper)

  res$intervals <- ci_df
  class(res) <- 'confint.simode'
  return(res)
}

calc_confint <- function(par_vals, nll_vals, max_delta, par_min, par_max) {

  ci_lower <- NA
  ci_upper <- NA

  delta_nll <- nll_vals-min(nll_vals)
  interp1 <- approx(x=par_vals, y=delta_nll, n=100)
  x_vals <- interp1$x
  y_vals <- interp1$y
  t_vals <- as.numeric(y_vals < max_delta)
  b_vals <- which(t_vals==1)

  if(length(b_vals) >= 1) {
    ind.min <- min(b_vals)
    ind.max <- max(b_vals)
    if(ind.min > 1)
      ci_lower <- mean(x_vals[(ind.min-1):ind.min])
    if(ind.max < length(x_vals))
      ci_upper <- mean(x_vals[ind.max:(ind.max+1)])
  }
  if(is.na(ci_lower) && !is.null(par_min) && par_vals[1]==par_min)
    ci_lower <- par_min
  if(is.na(ci_upper) && !is.null(par_max) && par_vals[length(par_vals)]==par_max)
    ci_upper <- par_max
  res <- list(lower=ci_lower,upper=ci_upper)
  return (res)
}

#' Plot confidence intervals for the model parameters
#'
#' Plot confidence intervals for the model parameters of a \code{simode} object,
#' calculated based on likelihood profiles
#' @param x \code{confint.simode} object returned by a call to \code{\link{confint}}
#' @param which Which parameters to plot the confidence intervals for.
#' If empty, the plot will include all of the parameters in x.
#' @param pars_true The true parameter values (if are known).
#' @param legend Whether or not to add a figure legend.
#' @param cols List of colors for each element of the plot.
#' @param ... Additional argument(s) for methods.
#' @export
#'
plot.confint.simode<- function(x, which=NULL, pars_true=NULL, legend=F,
                              cols=list(fit="blue",true="black"), ...) {

  pars <- which
  if(is.null(pars))
    pars <- as.vector(x$intervals$par)

  if(!is.null(pars_true))
    stopifnot(is.numeric(pars_true), length(pars_true)==length(pars))

  lp <- length(pars)
  pind <- which(as.vector(x$intervals$par) %in% pars)

  y_min <- Inf
  y_max <- -Inf
  nls_pars_est <- apply(as.matrix(1:lp),1,function(i) x$intervals[pind[i],]$nls_est)
  lower_vals <- apply(as.matrix(1:lp),1,function(i) x$intervals[pind[i],]$lower)
  upper_vals <- apply(as.matrix(1:lp),1,function(i) x$intervals[pind[i],]$upper)
  names(lower_vals) <- pars
  names(upper_vals) <- pars
  y_min <- min(c(y_min,lower_vals),na.rm=T)
  y_max <- max(c(y_max,upper_vals),na.rm=T)
  if(!is.null(pars_true)) {
    y_min <- min(c(y_min,pars_true),na.rm=T)
    y_max <- max(c(y_max,pars_true),na.rm=T)
  }

  legend_text <- c()
  legend_cols <- c()
  legend_lty <- c()
  legend_pch <- c()

  if(legend) {
    old.par <- par(xpd = T, mar = par()$mar + c(0,0,0,4))
    on.exit(par(old.par))
  }

  plot(-100, -100, xlim=c(0.5,lp+0.5), ylim=c(y_min,y_max), xaxt="n",
       xlab='parameter', ylab='value', ...)

  if(!is.null(pars_true)) {
    points(x=1:lp, pars_true, col=cols[['true']], pch=4, xaxt="n")
    legend_text <- c(legend_text,'true')
    legend_cols <- c(legend_cols,cols[['true']])
    legend_lty <- c(legend_lty,0)
    legend_pch <- c(legend_pch,4)
  }
  x_vals <- 1:lp
  points(x=x_vals, nls_pars_est, col=cols[['fit']], xaxt="n")
  for(i in 1:lp) {
    if(!is.na(lower_vals[[i]]) && !is.na(upper_vals[[i]])) {
      lines(x=c(x_vals[i],x_vals[i]),
            y=c(lower_vals[[i]], upper_vals[[i]]),
            col=cols[['fit']])
    }
  }
  legend_text <- c(legend_text,'nls_est','ci')
  legend_cols <- c(legend_cols,cols[['fit']],cols[['fit']])
  legend_lty <- c(legend_lty,0,1)
  legend_pch <- c(legend_pch,1,NA)

  axis(side=1,at=1:lp,labels=pars)

  if(legend) {
    legend(lp+1, y_max, cex=0.8, bty="n",
           legend_text ,col=legend_cols, lty=legend_lty, pch=legend_pch)
  }

}

#' Print method for \code{confint.simode} objects
#'
#' @param x The \code{confint.simode} object.
#' @param ... Additional argument(s) for methods.
#' @export
print.confint.simode <- function(x, ...) {

  cat("call:\n")
  print(x$call)
  cat("level:\n")
  cat(x$level,"\n")
  cat("intervals:\n")
  print(x$intervals)
}
