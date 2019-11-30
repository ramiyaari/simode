
simode_multi_worker <- function(...) {

  arg_list <- list(...)
  cl <- arg_list$cl
  parallel::clusterExport(cl, "arg_list", envir=environment())
  parallel::clusterExport(cl, ls("package:simode"), envir=environment())

  obs_sets <- length(arg_list$obs)
  simode_fits <- parallel::parLapply(cl, 1:obs_sets, function(i) {
    with(arg_list, {
      simode_ctrl$job_id <- i
      args <-
        c(list(equations=equations, pars=pars, time=time, obs=obs[[i]],
            nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
            fixed=fixed, start=start, lower=lower, upper=upper,
            im_method=im_method, decouple_equations=decouple_equations,
            gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
            simode_ctrl=simode_ctrl), extra_args)
      do.call("simode", args)
    })
  })
  return (simode_fits)
}



simode_multi <- function(equations, pars, time, obs,
                        nlin_pars, likelihood_pars,
                        fixed, start, lower, upper,
                        im_method, decouple_equations,
                        gen_obs, calc_nll, scale_pars, simode_ctrl, ...) {

  if(simode_ctrl$obs_sets_fit == 'together')
    return (simode_multi_together(equations=equations, pars=pars,  time=time, obs=obs,
                         nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
                         fixed=fixed, start=start, lower=lower, upper=upper,
                         im_method=im_method, decouple_equations=decouple_equations,
                         gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
                         simode_ctrl=simode_ctrl, ...))

  else if(simode_ctrl$obs_sets_fit == 'separate_x0')
    return (simode_multi_separate_x0(equations=equations, pars=pars,  time=time, obs=obs,
                                  nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
                                  fixed=fixed, start=start, lower=lower, upper=upper,
                                  im_method=im_method, decouple_equations=decouple_equations,
                                  gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
                                  simode_ctrl=simode_ctrl, ...))

  return (simode_multi_separate(equations=equations, pars=pars,  time=time, obs=obs,
                                nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
                                fixed=fixed, start=start, lower=lower, upper=upper,
                                im_method=im_method, decouple_equations=decouple_equations,
                                gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
                                simode_ctrl=simode_ctrl, ...))
}

simode_multi_together <- function(equations, pars, time, obs,
                         nlin_pars, likelihood_pars,
                         fixed, start, lower, upper,
                         im_method, decouple_equations,
                         gen_obs, calc_nll, scale_pars,
                         simode_ctrl, ...) {
  obs_mean <- list()
  obs_sets <- length(obs)
  for(i in 1:length(obs[[1]])) {
    obs1 <- lapply(1:obs_sets, function(j) obs[[j]][i])
    obs1 <- matrix(unlist(obs1), ncol = obs_sets)
    obs_mean[[i]] <- unlist(lapply(1:nrow(obs1), function(k) mean(obs1[k,])))
  }
  names(obs_mean) <- names(obs[[1]])
  obs <- obs_mean

  return (simode(equations=equations, pars=pars, time=time, obs=obs,
           nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
           fixed=fixed, start=start, lower=lower, upper=upper,
           im_method=im_method, decouple_equations=decouple_equations,
           gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
           simode_ctrl=simode_ctrl, ...))
}

simode_multi_separate_x0 <- function(equations, pars, time, obs,
                                  nlin_pars, likelihood_pars,
                                  fixed, start, lower, upper,
                                  im_method, decouple_equations,
                                  gen_obs, calc_nll, scale_pars, simode_ctrl, ...) {
  sep <- simode_ctrl$decouple_sep
  eq_names <- names(equations)
  vars <- names(obs[[1]])
  equations_m <- c()
  obs_m <- c()
  time_m <- c()
  if(!is.list(time))
    time_m <- time
  obs_sets <- length(obs)
  for(i in 1:obs_sets) {
    vars1 <- paste0(vars,sep,i)
    names(vars1) <- vars
    obs1 <- obs[[i]]
    names(obs1) <- vars1
    obs_m <- c(obs_m,obs1)
    if(is.list(time)) {
      stopifnot(length(time)==length(obs1))
      time1 <- time
      names(time1) <- vars1
      time_m <- c(time_m,time1)
    }
    equations1 <- fix_pars(equations,vars1)
    names(equations1) <- paste0(eq_names,sep,i)
    equations_m <- c(equations_m,equations1)
  }
  pars_m <- pars
  pars_x0 <- intersect(pars,eq_names)
  if(!pracma::isempty(pars_x0)) {
    pars_x0_ex <- paste0(pars_x0,sep,unlist(lapply(1:obs_sets,function(i) rep(i,length(pars_x0)))))
    pars_m <- c(setdiff(pars,pars_x0),pars_x0_ex)
  }
  nlin_pars_m <- nlin_pars
  nlin_pars_x0 <- intersect(nlin_pars,eq_names)
  if(!pracma::isempty(nlin_pars_x0)) {
    nlin_pars_x0_ex <- paste0(nlin_pars_x0,sep,unlist(lapply(1:obs_sets,function(i) rep(i,length(nlin_pars_x0)))))
    nlin_pars_m <- c(setdiff(nlin_pars,nlin_pars_x0),nlin_pars_x0_ex)
  }
  if(!is.list(fixed)) {
    fixed_m <- fixed
    fixed_x0_names <- intersect(names(fixed),eq_names)
    if(!pracma::isempty(fixed_x0_names)) {
      fixed_x0 <- rep(fixed[fixed_x0_names],obs_sets)
      names(fixed_x0) <-
        paste0(fixed_x0_names,sep,
               unlist(lapply(1:obs_sets,function(i) rep(i,length(equations)))))
      fixed_m <- c(fixed[setdiff(names(fixed),fixed_x0_names)],fixed_x0)
    }
  }
  else {
    stopifnot(length(fixed)==obs_sets)
    fixed_m <- c()
    for(i in 1:obs_sets) {
      fixed_x0_names <- intersect(names(fixed[[i]]),eq_names)
      if(!pracma::isempty(fixed_x0_names)) {
        fixed_x0 <- fixed[[i]][fixed_x0_names]
        names(fixed_x0) <- paste0(fixed_x0_names,sep,i)
        fixed_m <- c(fixed_m,fixed[[i]][setdiff(names(fixed[[i]]),fixed_x0_names)],fixed_x0)
      }
      else
        fixed_m <- c(fixed_m,fixed[[i]])
    }
  }

  start_m <- start
  start_x0_names <- intersect(names(start),eq_names)
  if(!pracma::isempty(start_x0_names)) {
    start_x0 <- rep(start[start_x0_names],obs_sets)
    names(start_x0) <-
      paste0(start_x0_names,sep,
             unlist(lapply(1:obs_sets,function(i) rep(i,length(equations)))))
    start_m <- c(start[setdiff(names(start),start_x0_names)],start_x0)
  }
  lower_m <- lower
  lower_x0_names <- intersect(names(lower),eq_names)
  if(!pracma::isempty(lower_x0_names)) {
    lower_x0 <- rep(lower[lower_x0_names],obs_sets)
    names(lower_x0) <-
      paste0(lower_x0_names,sep,
             unlist(lapply(1:obs_sets,function(i) rep(i,length(equations)))))
    lower_m <- c(lower[setdiff(names(lower),lower_x0_names)],lower_x0)
  }
  upper_m <- upper
  upper_x0_names <- intersect(names(upper),eq_names)
  if(!pracma::isempty(upper_x0_names)) {
    upper_x0 <- rep(upper[upper_x0_names],obs_sets)
    names(upper_x0) <-
      paste0(upper_x0_names,sep,
             unlist(lapply(1:obs_sets,function(i) rep(i,length(equations)))))
    upper_m <- c(upper[setdiff(names(upper),upper_x0_names)],upper_x0)
  }

  x <- simode(equations=equations_m, pars=pars_m, time=time_m, obs=obs_m,
              nlin_pars=nlin_pars_m, likelihood_pars=likelihood_pars,
              fixed=fixed_m, start=start_m, lower=lower_m, upper=upper_m,
              im_method=im_method, decouple_equations=decouple_equations,
              gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
              simode_ctrl=simode_ctrl, ...)

  if(is.null(x))
    return (NULL)

  x$equations <- equations
  x$pars <- pars
  x$nlin_pars <- nlin_pars
  x$fixed <- fixed
  x$start <- start
  x$lower <- lower
  x$upper <- upper
  x_list <- list()
  for(i in 1:obs_sets) {
    x_list[[i]] <- x
    x_list[[i]]$obs <- obs[[i]]
    x_list[[i]]$x0 <- NULL
    x_list[[i]]$x0[pars_x0] <- x$x0[paste0(pars_x0,sep,i)]
    fixed_x0_names <- intersect(names(fixed[[i]]),eq_names)
    if(!pracma::isempty(fixed_x0_names)) {
      x_list[[i]]$x0[fixed_x0_names] <- fixed[[i]][fixed_x0_names]
    }
    if(!is.null(x$im_pars_est)) {
      x_list[[i]]$im_pars_est <- x$im_pars_est[setdiff(pars,pars_x0)]
      x_list[[i]]$im_pars_est[pars_x0] <- x$im_pars_est[paste0(pars_x0,sep,i)]
      x_list[[i]]$im_pars_est <- x_list[[i]]$im_pars_est[pars]
      x_list[[i]]$im_smooth$val <- x$im_smooth$val[,paste0(eq_names,sep,i)]
      colnames(x_list[[i]]$im_smooth$val) <- eq_names
    }
    if(!is.null(x$nls_pars_est)) {
      x_list[[i]]$nls_pars_est <- x$nls_pars_est[setdiff(pars,pars_x0)]
      x_list[[i]]$nls_pars_est[pars_x0] <- x$nls_pars_est[paste0(pars_x0,sep,i)]
      x_list[[i]]$nls_pars_est <- x_list[[i]]$nls_pars_est[pars]
    }
  }
  class(x_list) <- 'list.simode'
  return (x_list)
}

simode_multi_separate <- function(equations, pars, time, obs,
                         nlin_pars, likelihood_pars,
                         fixed, start, lower, upper,
                         im_method, decouple_equations,
                         gen_obs, calc_nll, scale_pars, simode_ctrl, ...) {

  if(simode_ctrl$parallel==T && !requireNamespace("parallel", quietly=T)) {
    warning("parallel package not installed - running sequentially",
            immediate.=T)
    simode_ctrl$parallel <- F
  }

  log_file <- NULL
  if(simode_ctrl$save_to_log) {
    log_file <- file.path(tempdir(),"simode.log")
    simode_ctrl$save_to_log <- F
  }

  simode_fits <- NULL
  obs_sets <- length(obs)
  tryCatch(
    {
      if(simode_ctrl$parallel==F) {
        if(!is.null(log_file)) {
          sink(file=log_file, append=F)
          cat(noquote((paste0("Call to simode on [", Sys.time(), "]:\n"))))
        }
        simode_fits <- list()
        for(i in 1:obs_sets) {
          simode_ctrl$job_id <- i
          simode_fits[[i]] <-
              simode(equations=equations, pars=pars,  time=time, obs=obs[[i]],
                    nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
                    fixed=fixed, start=start, lower=lower, upper=upper,
                    im_method=im_method, decouple_equations=decouple_equations,
                    gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
                    simode_ctrl=simode_ctrl, ...)
        }
      }
      else {
        cl <- parallel::makeCluster(obs_sets, outfile=log_file)
        extra_args <- list(...)
        simode_fits <- simode_multi_worker(
          cl=cl, equations=equations, pars=pars, time=time, obs=obs,
          nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
          fixed=fixed, start=start, lower=lower, upper=upper,
          im_method=im_method, decouple_equations=decouple_equations,
          gen_obs=gen_obs, calc_nll=calc_nll, scale_pars=scale_pars,
          simode_ctrl=simode_ctrl, extra_args=extra_args)
      }
    },
    error = function(e) { print(e) },
    finally={
      if(simode_ctrl$parallel==T && !is.null(cl))
        parallel::stopCluster(cl)
      }
  )

  if(!is.null(simode_fits))
    class(simode_fits) <- "list.simode"

  return (simode_fits)
}

#' Summary method for \code{list.simode} objects
#'
#' @param object \code{list.simode} object returned by a call to \code{\link{simode}}
#' with \code{obs_sets>1}
#' @param sum_mean_sd Whether to calculate mean and standard deviation
#' for the parameter estimates in the fits included in the given \code{object}.
#' To be used when \code{object} is the result of fitting monte-carlo simulations.
#' @param pars_true The true parameter values (relevant only for
#' when \code{sum_mean_sd=T}). When given, the summary will also include the
#' bias and RMSE for each parameter estimate.
#' @param digits The number of significant digits to use.
#' @param ... Additional argument(s) for methods.
#' @return The mean and standard deviation for the
#' loss values and parameter estimates obtained from the integral-matching
#' and nonlinear least squares optimizations. If \code{pars_true} is given, then
#' will also calculate bias and RMSE for the parameter estimates.
#' @export
#'
summary.list.simode <- function(object, sum_mean_sd=F, pars_true=NULL,
                                digits=max(3, getOption("digits")-3), ...) {

  pars <- object[[1]]$pars
  obs_sets_fit <- object[[1]]$ctrl$obs_sets_fit


  if(!is.null(pars_true)) {
    if(!sum_mean_sd)
      warning('parameter pars_true is only relevant when sum_mean_sd==T')
    else
      stopifnot(!is.null(names(pars_true)),all(names(pars_true) %in% pars))
  }

  simode_num <- length(object)
  im_est <- matrix(NA,nrow=simode_num,ncol=length(pars))
  nls_est <- matrix(NA,nrow=simode_num,ncol=length(pars))
  im_loss_vals <- rep(NA,simode_num)
  nls_loss_vals <- rep(NA,simode_num)

  for(i in 1:simode_num) {
    if(!is.null(object[[i]]$im_pars_est))
      im_est[i,] <- object[[i]]$im_pars_est
    if(!is.null(object[[i]]$nls_pars_est))
      nls_est[i,] <- object[[i]]$nls_pars_est
    if(!is.null(object[[i]]$im_loss))
      im_loss_vals[i] <- object[[i]]$im_loss
    if(!is.null(object[[i]]$nls_loss))
      nls_loss_vals[i] <- object[[i]]$nls_loss
  }

  summary <- list()

  if(!sum_mean_sd) {
    if(!all(is.na(im_est))) {
      summary$im_est <- signif(im_est,digits)
      colnames(summary$im_est) <- pars
    }
    if(!all(is.na(nls_est))) {
      summary$nls_est <- signif(nls_est,digits)
      colnames(summary$nls_est) <- pars
    }
    if(obs_sets_fit == "separate") {
      if(any(!is.na(im_loss_vals)))
        summary$im_loss <- signif(sum(im_loss_vals),digits)
      if(any(!is.na(nls_loss_vals)))
        summary$nls_loss <- signif(sum(nls_loss_vals),digits)
    }
    else
    {
      if(any(!is.na(im_loss_vals)))
        summary$im_loss <- signif(im_loss_vals[1],digits)
      if(any(!is.na(nls_loss_vals)))
        summary$nls_loss <- signif(nls_loss_vals[1],digits)
    }
  }
  else {
    im_loss_mean <- signif(mean(im_loss_vals),digits)
    im_loss_sd <- signif(stats::sd(im_loss_vals),digits)
    nls_loss_mean <- signif(mean(nls_loss_vals),digits)
    nls_loss_sd <- signif(stats::sd(nls_loss_vals),digits)
    if(!is.na(im_loss_mean)) {
      summary$loss$im_mean <- im_loss_mean
      summary$loss$im_sd <- im_loss_sd
    }
    if(!is.na(nls_loss_mean)) {
      summary$loss$nls_mean <- nls_loss_mean
      summary$loss$nls_sd <- nls_loss_sd
    }
    summary$loss <- as.data.frame(summary$loss)

    im_est_mean <- signif(apply(im_est,2,mean),digits)
    im_est_sd <- signif(apply(im_est,2,stats::sd),digits)
    names(im_est_mean) <- pars
    nls_est_mean <- signif(apply(nls_est,2,mean),digits)
    nls_est_sd <- signif(apply(nls_est,2,stats::sd),digits)
    names(nls_est_mean) <- pars

    summary$est <- data.frame(par=pars)
    if(!is.null(pars_true)) {
      summary$est$true <- signif(pars_true[pars], digits)
      im_est_bias <- signif(im_est_mean - pars_true[pars],digits)
      nls_est_bias <- signif(nls_est_mean - pars_true[pars],digits)
      im_est_err <- im_est-matrix(rep(pars_true[pars],simode_num),nrow=simode_num,byrow=T)
      im_est_rmse <- signif(sqrt(apply(im_est_err^2,2,mean)),digits)
      nls_est_err <- nls_est-matrix(rep(pars_true[pars],simode_num),nrow=simode_num,byrow=T)
      nls_est_rmse <- signif(sqrt(apply(nls_est_err^2,2,mean)),digits)
    }

    if(any(!is.na(im_est_mean))) {
      summary$est$im_mean <- im_est_mean
      summary$est$im_sd <- im_est_sd
      if(!is.null(pars_true)) {
        summary$est$im_bias <- im_est_bias
        summary$est$im_rmse <- im_est_rmse
      }
    }
    if(any(!is.na(nls_est_mean))) {
      summary$est$nls_mean <- nls_est_mean
      summary$est$nls_sd <- nls_est_sd
      if(!is.null(pars_true)) {
        summary$est$nls_bias <- nls_est_bias
        summary$est$nls_rmse <- nls_est_rmse
      }
    }
  }

  class(summary) <- "summary.list.simode"
  return(summary)
}


#' Plot the fit/estimates of a \code{list.simode} object
#'
#' Plot the fit or parameter estimates obtained from a call to \code{simode}
#' with \code{obs_sets>1}.
#'
#' @param x \code{list.simode} object.
#' @param type Type of plot - 'fit' to plot the fitted equations
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
#' @param plot_mean_sd Plot the mean and standard deviation for the fit/estimates
#' of the \code{simode} objects in \code{x}. To be used when \code{x} is the
#' result of fitting monte-carlo simulations.
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
plot.list.simode <-
  function(x, type=c('fit','est'), show=c('nls','im','both'),
           which=NULL, pars_true=NULL, time=NULL,
           plot_mean_sd=F, plot_im_smooth=F,
           legend=F, mfrow=par('mfrow'),
           cols=list(nls_fit="blue",im_fit="green", true="black",
                     obs="red", im_smooth="magenta"), ...) {

    type <- match.arg(type)

    show <- match.arg(show)
    if((show=='nls' || show=='both') && is.null(x[[1]]$nls_pars_est)) {
      warning("no nls estimates in this simode object - will show only \'im\' fit\\estimates")
      show <- 'im'
    }
    if((show=='im' || show=='both') && is.null(x[[1]]$im_pars_est)) {
      warning("no im estimates in this simode object - will show only \'nls\' fit\\estimates")
      show <- 'nls'
    }

    if(length(x)>1 && !plot_mean_sd) {
      for(i in 1:length(x)) {
        xi <- x[i]
        class(xi) <- 'list.simode'
        if(type == 'fit')
          plot_list_simode_fit(x=xi,show=show,which=which,pars_true=pars_true,
                               time=time,plot_im_smooth=plot_im_smooth,plot_mean_sd=plot_mean_sd,
                               legend=legend, mfrow=mfrow, cols=cols, fit_num=i, ...)
        else #if(type == 'est')
          plot_list_simode_est(x=xi,show=show,which=which,pars_true=pars_true,
                               plot_mean_sd=plot_mean_sd, legend=legend, cols=cols, fit_num=i, ...)
      }
    }
    else {
      if(type == 'fit')
        plot_list_simode_fit(x=x,show=show,which=which,pars_true=pars_true,
                             time=time,plot_im_smooth=plot_im_smooth,plot_mean_sd=plot_mean_sd,
                             legend=legend, mfrow=mfrow, cols=cols, ...)
      else #if(type == 'est')
        plot_list_simode_est(x=x,show=show,which=which,pars_true=pars_true,
                             plot_mean_sd=plot_mean_sd,legend=legend, cols=cols, ...)
    }
}


plot_list_simode_est <-
        function(x, show=c('nls','im','both'),
                 which=NULL, pars_true=NULL, plot_mean_sd=F, legend=F,
                 cols=list(nls_fit="blue",im_fit="green", true="black"),
                 fit_num=NULL, ...) {

  show <- match.arg(show)
  plot_nls_est <- (show=='nls' || show=='both')
  plot_im_est <- (show=='im' || show=='both')

  pars <- which
  if(is.null(pars))
    pars <- x[[1]]$pars

  if(length(x)==1)
    plot_mean_sd <- F

  if(!is.null(pars_true)) {
      stopifnot(is.numeric(pars_true))
      stopifnot(all(pars %in% names(pars_true)))
  }

  lp <- length(pars)
  pind <- which(as.vector(x[[1]]$pars) %in% pars)

  if(plot_mean_sd) {
    sum <- summary(x, sum_mean_sd=T)
    im_est_mean <- apply(as.matrix(1:lp),1,function(i) sum$est[pind[i],]$im_mean)
    im_est_sd <- apply(as.matrix(1:lp),1,function(i) sum$est[pind[i],]$im_sd)
    nls_est_mean <- apply(as.matrix(1:lp),1,function(i) sum$est[pind[i],]$nls_mean)
    nls_est_sd <- apply(as.matrix(1:lp),1,function(i) sum$est[pind[i],]$nls_sd)
    names(im_est_mean) <- pars
    names(im_est_sd) <- pars
    names(nls_est_mean) <- pars
    names(nls_est_sd) <- pars
    im_est_sd[which(is.na(im_est_sd))] <- 0
    nls_est_sd[which(is.na(nls_est_sd))] <- 0
  }

  y_min <- Inf
  y_max <- -Inf
  if(!is.null(pars_true)) {
      y_min <- min(pars_true[pars],na.rm=T)
      y_max <- max(pars_true[pars],na.rm=T)
  }
  if(plot_mean_sd) {
    if(plot_nls_est) {
      min_nls <- nls_est_mean[pars] - nls_est_sd[pars]
      max_nls <- nls_est_mean[pars] + nls_est_sd[pars]
      y_min <- min(c(y_min,min_nls),na.rm=T)
      y_max <- max(c(y_max,max_nls),na.rm=T)
    }
    if(plot_im_est) {
      min_im <- im_est_mean[pars] - im_est_sd[pars]
      max_im <- im_est_mean[pars] + im_est_sd[pars]
      y_min <- min(c(y_min,min_im),na.rm=T)
      y_max <- max(c(y_max,max_im),na.rm=T)
    }
  }
  else {
    if(plot_nls_est) {
      min_nls <- min(unlist(lapply(1:length(x),function(i) min(x[[i]]$nls_pars_est[pars]))))
      max_nls <- max(unlist(lapply(1:length(x),function(i) max(x[[i]]$nls_pars_est[pars]))))
      y_min <- min(c(y_min,min_nls),na.rm=T)
      y_max <- max(c(y_max,max_nls),na.rm=T)
    }
    if(plot_im_est) {
      min_im <- min(unlist(lapply(1:length(x),function(i) min(x[[i]]$im_pars_est[pars]))))
      max_im <- max(unlist(lapply(1:length(x),function(i) max(x[[i]]$im_pars_est[pars]))))
      y_min <- min(c(y_min,min_im),na.rm=T)
      y_max <- max(c(y_max,max_im),na.rm=T)
    }
  }

  legend_text <- c()
  legend_cols <- c()
  legend_lty <- c()
  legend_pch <- c()

  if(legend) {
    old.par <- par(xpd=T, mar=c(4,4,2,6))
    on.exit(par(old.par))
  }

  ylab <- ifelse(is.null(fit_num),'value',paste0('value (fit-',fit_num,')'))
  plot(-100, -100, xlim=c(0.5,lp+0.5), ylim=c(y_min,y_max), xaxt="n",
       xlab='parameter', ylab=ylab, ...)

  if(!is.null(pars_true)) { # && !is.list(pars_true)) {
    points(x=1:lp, pars_true[pars], col=cols[['true']], pch=4, xaxt="n",...)
    legend_text <- c(legend_text,'true')
    legend_cols <- c(legend_cols,cols[['true']])
    legend_lty <- c(legend_lty,0)
    legend_pch <- c(legend_pch,4)
  }

  if(plot_mean_sd) {
    if(plot_im_est && !is.null(im_est_mean)) {
      x_vals <- 1:lp
      if(plot_nls_est) # || is.list(pars_true))
        x_vals <- x_vals-0.1
      points(x=x_vals, im_est_mean[pars], col=cols[['im_fit']], xaxt="n",...)
      errorbar(x=x_vals, im_est_mean[pars], yerr=im_est_sd[pars],
               bar.col=cols[['im_fit']], grid=F, add=T, xaxt="n")
      legend_text <- c(legend_text,'im_mean','im_sd')
      legend_cols <- c(legend_cols,cols[['im_fit']],cols[['im_fit']])
      legend_lty <- c(legend_lty,0,1)
      legend_pch <- c(legend_pch,1,NA)
    }
    if(plot_nls_est && !is.null(nls_est_mean)) {
      x_vals <- 1:lp
      if(plot_im_est)
        x_vals <- x_vals+0.1
      points(x=x_vals, nls_est_mean[pars], col=cols[['nls_fit']], xaxt="n",...)
      errorbar(x=x_vals, nls_est_mean[pars], yerr=nls_est_sd[pars],
               bar.col=cols[['nls_fit']], grid=F, add=T, xaxt="n")
      legend_text <- c(legend_text,'nls_mean','nls_sd')
      legend_cols <- c(legend_cols,cols[['nls_fit']],cols[['nls_fit']])
      legend_lty <- c(legend_lty,0,1)
      legend_pch <- c(legend_pch,1,NA)
    }
  }
  else {
    for(i in 1:length(x)) {
      if(plot_im_est) {
        x_vals <- 1:lp
        if(plot_nls_est)
          x_vals <- x_vals-0.1
        points(x=x_vals, x[[i]]$im_pars_est[pars], col=cols[['im_fit']], xaxt="n",...)
        if(i==1) {
          legend_text <- c(legend_text,'im_est')
          legend_cols <- c(legend_cols,cols[['im_fit']])
          legend_lty <- c(legend_lty,0)
          legend_pch <- c(legend_pch,1)
        }
      }
      if(plot_nls_est) {
        x_vals <- 1:lp
        if(plot_im_est)
          x_vals <- x_vals+0.1
        points(x=x_vals, x[[i]]$nls_pars_est[pars], col=cols[['nls_fit']], xaxt="n",...)
        if(i==1) {
          legend_text <- c(legend_text,'nls_est')
          legend_cols <- c(legend_cols,cols[['nls_fit']])
          legend_lty <- c(legend_lty,0)
          legend_pch <- c(legend_pch,1)
        }
      }
    }
  }
  axis(side=1,at=1:lp,labels=pars)

  if(legend) {
    legend(lp+1, y_max, cex=0.8, bty="n",
           legend_text ,col=legend_cols, lty=legend_lty, pch=legend_pch)
  }

}


plot_list_simode_fit <-
  function(x, show=c('nls','im','both'),
           which=NULL, pars_true=NULL, time=NULL,
           plot_mean_sd=F, plot_im_smooth=F, legend=F, mfrow=par('mfrow'),
           cols=list(nls_fit="blue",im_fit="green", true="black",
                     obs="red", im_smooth="magenta"),
           fit_num=NULL, main=NULL,...) {


  show <- match.arg(show)
  plot_nls_fit <- (show=='nls' || show=='both')
  plot_im_fit <- (show=='im' || show=='both')

  x0.na <- names(x[[1]]$x0[is.na(x[[1]]$x0)])
  xvars <- setdiff(names(x[[1]]$obs),names(x[[1]]$equations))

  vars <- which
  if(is.null(vars))
    vars <- setdiff(names(x[[1]]$obs),xvars)

  if(!is.null(main)) {
    if(length(main)==1) {
      main <- rep(main,length(vars))
    }
    else {
      stopifnot(length(main)==length(vars))
    }
    names(main) <- vars
  }

  if(length(x)==1)
    plot_mean_sd <- F

  if(!is.null(pars_true)) {
      stopifnot(is.numeric(pars_true))
      stopifnot(all(setdiff(x[[1]]$pars,x[[1]]$likelihood_pars) %in% names(pars_true)))
  }

  if(!is.list(x[[1]]$time)) {
    if(is.null(time))
      time <- x[[1]]$time
    time_obs <- rep(list(x[[1]]$time),length(x[[1]]$obs[vars]))
  }
  else {
    if(is.null(time))
      time <- sort(unique(unlist(x[[1]]$time)))
    time_obs <- x[[1]]$time[which(names(x[[1]]$obs) %in% vars)]
  }
  names(time_obs) <- vars

  y_min <- unlist(lapply(vars, function(i) min(x[[1]]$obs[[i]])))
  y_max <- unlist(lapply(vars, function(i) max(x[[1]]$obs[[i]])))
  names(y_min) <- names(x[[1]]$obs[vars])
  names(y_max) <- names(x[[1]]$obs[vars])

  legend_text <- c('obs')
  legend_cols <- cols[['obs']]

  if(!is.null(pars_true)) {
    x0_true <- x[[1]]$x0
    x0_true[x0.na] <- pars_true[x0.na]
    pars_true <- pars_true[setdiff(names(pars_true),x0.na)]
    model_true <- solve_ode(x[[1]]$equations, pars_true, x0_true, time, x[[1]]$obs[xvars])
    if(is.null(model_true)) {
      stop('Error running model with pars_true values')
    }
    legend_text <- c(legend_text, 'true')
    legend_cols <- c(legend_cols, cols[['true']])
    y_min1 <- apply(as.matrix(model_true[,vars]), 2, function(col) min(col))
    y_max1 <- apply(as.matrix(model_true[,vars]), 2, function(col) max(col))
    y_min <- pmin(y_min,y_min1)
    y_max <- pmax(y_max,y_max1)
  }

  model_nls_est <- list()
  if(plot_nls_fit) {
    for(i in 1:length(x)) {
      simode_obj <- x[[i]]
      simode_obj$time <- time
      model_nls_est[[i]] <- solve_ode2(simode_obj,'nls')$nls
      if(is.null(model_nls_est[[i]])) {
        fn <- ifelse(!is.null(fit_num),fit_num,i)
        stop(paste0('Error running model with nls estimates of fit #',fn))
      }
    }
    model_nls_est_mean <- matrix(0,nrow=length(time), ncol=length(vars))
    model_nls_est_sd <- matrix(0,nrow=length(time), ncol=length(vars))
    colnames(model_nls_est_mean) <- vars
    colnames(model_nls_est_sd) <- vars
    for(var in vars) {
      nls_est_mat <- matrix(0,nrow=length(time), ncol=length(x))
      for(i in 1:length(x)) {
        nls_est_mat[,i] <- model_nls_est[[i]][,var]
      }
      model_nls_est_mean[,var] <- apply(nls_est_mat,1,mean)
      if(length(x)>1)
        model_nls_est_sd[,var] <- apply(nls_est_mat,1,stats::sd)
    }
    text <- ifelse(plot_mean_sd,'nls_mean_sd','nls_fit')
    legend_text <- c(legend_text, text)
    legend_cols <- c(legend_cols, cols[['nls_fit']])
    y_min1 <- apply(as.matrix(model_nls_est_mean[,vars]-model_nls_est_sd[,var]), 2, function(col) min(col))
    y_max1 <- apply(as.matrix(model_nls_est_mean[,vars]+model_nls_est_sd[,var]), 2, function(col) max(col))
    y_min <- pmin(y_min,y_min1)
    y_max <- pmax(y_max,y_max1)
  }

  model_im_est <- list()
  if(plot_im_fit) {
    for(i in 1:length(x)) {
      simode_obj <- x[[i]]
      simode_obj$time <- time
      model_im_est[[i]] <- solve_ode2(simode_obj,'im')$im
      if(is.null(model_im_est[[i]])) {
        fn <- ifelse(!is.null(fit_num),fit_num,i)
        stop(paste0('Error running model with im estimates of fit #',fn))
      }
    }
    model_im_est_mean <- matrix(0,nrow=length(time), ncol=length(vars))
    model_im_est_sd <- matrix(0,nrow=length(time), ncol=length(vars))
    colnames(model_im_est_mean) <- vars
    colnames(model_im_est_sd) <- vars
    for(var in vars) {
      im_est_mat <- matrix(0,nrow=length(time), ncol=length(x))
      for(i in 1:length(x)) {
        im_est_mat[,i] <- model_im_est[[i]][,var]
      }
      model_im_est_mean[,var] <- apply(im_est_mat,1,mean)
      if(length(x)>1)
        model_im_est_sd[,var] <- apply(im_est_mat,1,stats::sd)
    }
    text <- ifelse(plot_mean_sd,'im_mean_sd','im_fit')
    legend_text <- c(legend_text, text)
    legend_cols <- c(legend_cols, cols[['im_fit']])
    y_min1 <- apply(model_im_est_mean[,vars]-model_im_est_sd[,var], 2, function(col) min(col))
    y_max1 <- apply(model_im_est_mean[,vars]+model_im_est_sd[,var], 2, function(col) max(col))
    y_min <- pmin(y_min,y_min1)
    y_max <- pmax(y_max,y_max1)
  }

  if(plot_im_smooth) {
    im_time <- x[[1]]$im_smooth$time
    model_im_smooth_mean <- matrix(0,nrow=length(im_time), ncol=length(vars))
    model_im_smooth_sd <- matrix(0,nrow=length(im_time), ncol=length(vars))
    colnames(model_im_smooth_mean) <- vars
    colnames(model_im_smooth_sd) <- vars
    im_smooth_mat <- matrix(0,nrow=length(im_time), ncol=length(x))
    for(var in vars) {
      for(i in 1:length(x)) {
        im_smooth_mat[,i] <- x[[i]]$im_smooth$val[,var]
      }
      model_im_smooth_mean[,var] <- apply(im_smooth_mat,1,mean)
      if(length(x)>1)
        model_im_smooth_sd[,var] <- apply(im_smooth_mat,1,stats::sd)
    }
    text <- ifelse(plot_mean_sd,'im_smooth_mean_sd','im_smooth')
    legend_text <- c(legend_text, text)
    legend_cols <- c(legend_cols, cols[['im_smooth']])
  }

  if(legend)
    old.par <- par(mfrow=mfrow, mar=c(4,4,2,6), xpd=T)
  else
    old.par <- par(mfrow=mfrow, mar=c(4,4,2,6))

  on.exit(par(old.par))

  for(var in vars) {

    ylab <- ifelse(is.null(fit_num),var,paste0(var,' (fit ',fit_num,')'))

    time_var <- time_obs[[var]]
    plot(-100,-100,xlab='time',ylab=ylab,
         xlim=c(time_var[1],time_var[length(time_var)]),
         ylim=c(y_min[var],y_max[var]), main=main[var], ...)

    for(i in 1:length(x)) {
      points(time_var,x[[i]]$obs[[var]],col=cols[['obs']],...)
    }

    if(!is.null(pars_true)) {
      lines(time,model_true[,var],col=cols[['true']],...)
    }

    if(plot_mean_sd) {
      if(plot_im_smooth) {
        lines(im_time,model_im_smooth_mean[,var],col=cols[['im_smooth']],...)
        errorbar(x=im_time, model_im_smooth_mean[,var], yerr=model_im_smooth_sd[,var],
                 bar.col=cols[['im_smooth']], grid=F, add=T, xaxt="n",...)
      }
      if(plot_im_fit) {
        lines(time,model_im_est_mean[,var],col=cols[['im_fit']],...)
        errorbar(x=time, model_im_est_mean[,var], yerr=model_im_est_sd[,var],
                 bar.col=cols[['im_fit']], grid=F, add=T, xaxt="n")
      }
      if(plot_nls_fit) {
        lines(time,model_nls_est_mean[,var],col=cols[['nls_fit']],...)
        errorbar(x=time, model_nls_est_mean[,var], yerr=model_nls_est_sd[,var],
                 bar.col=cols[['nls_fit']], grid=F, add=T, xaxt="n")
      }
    }
    else {
      for(i in 1:length(x)) {
        if(plot_im_smooth)
          lines(im_time,x[[i]]$im_smooth$val[,var],col=cols[['im_smooth']],...)
        if(plot_im_fit)
          lines(time,model_im_est[[i]][,var],col=cols[['im_fit']],...)
        if(plot_nls_fit)
          lines(time,model_nls_est[[i]][,var],col=cols[['nls_fit']],...)
      }
    }

    if(legend) {
      legend(1.05*time[length(time)],
             y_max[[var]],
             cex=0.8,
             bty="n",
             legend_text ,col=legend_cols,
             lty=c(0,rep(1,length(legend_text)-1)),
             pch=c(1,rep(NA,length(legend_text)-1)))

    }
  }
}
