
simode_multi_worker <- function(...) {

  arg_list <- list(...)
  cl <- arg_list$cl
  parallel::clusterExport(cl, "arg_list", envir=environment())
  parallel::clusterExport(cl, ls("package:simode"), envir=environment())

  simode_fits <- parallel::parLapply(cl, 1:arg_list$mc_sets, function(i) {
    with(arg_list, {
      simode_ctrl$job_id <- i
      args <-
        c(list(equations=equations, pars=pars, time=time, obs=obs[[i]],
            mc_sets=1, nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
            fixed=fixed, start=start, lower=lower, upper=upper,
            im_method=im_method, gen_obs=gen_obs, calc_nll=calc_nll,
            simode_ctrl=simode_ctrl), extra_args)
      do.call("simode", args)
    })
  })
  return (simode_fits)
}



simode_multi <- function(equations, pars, time, obs, mc_sets,
                        nlin_pars, likelihood_pars,
                        fixed, start, lower, upper,
                        im_method, gen_obs, calc_nll, simode_ctrl, ...) {

  if(simode_ctrl$parallel==T && !requireNamespace("parallel", quietly=T)) {
    warning("parallel package not installed - running sequentially",
            immediate.=T)
    simode_ctrl$parallel <- F
  }

  log_file <- NULL
  if(simode_ctrl$save_to_log) {
    log_file <- paste0(tempdir(),"\\simode.log")
    simode_ctrl$save_to_log <- F
  }

  simode_fits <- NULL
  tryCatch(
    {
      if(simode_ctrl$parallel==F) {
        if(!is.null(log_file)) {
          sink(file=log_file, append=F)
          cat(noquote((paste0("Call to simode on [", Sys.time(), "]:\n"))))
        }
        simode_fits <- list()
        for(i in 1:mc_sets) {
          simode_ctrl$job_id <- i
          simode_fits[[i]] <-
              simode(equations=equations, pars=pars,  time=time, obs=obs[[i]],
                    mc_sets=1, nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
                    fixed=fixed, start=start, lower=lower, upper=upper, im_method=im_method,
                    gen_obs=gen_obs, calc_nll=calc_nll, simode_ctrl=simode_ctrl, ...)
        }
      }
      else {
        cl <- parallel::makeCluster(mc_sets, outfile=log_file)
        extra_args <- list(...)
        simode_fits <- simode_multi_worker(
          cl=cl, equations=equations, pars=pars, time=time, obs=obs,
          mc_sets=mc_sets, nlin_pars=nlin_pars, likelihood_pars=likelihood_pars,
          fixed=fixed, start=start, lower=lower, upper=upper,
          im_method=im_method, gen_obs=gen_obs, calc_nll=calc_nll,
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
#' with \code{mc_sets>1}
#' @param pars_true The true parameter values (if known).
#' @param digits The number of significant digits to use.
#' @param ... Additional argument(s) for methods.
#' @return The mean and standard deviation for the
#' loss values and parameter estimates obtained from the integral-matching
#' and nonlinear least squares optimizations. If \code{pars_true} is given, then
#' will also calculate bias and RMSE for the parameter estimates.
#' @export
#'
summary.list.simode <- function(object, pars_true=NULL,
                                digits=max(3, getOption("digits")-3), ...) {

  pars <- object[[1]]$pars

  if(!is.null(pars_true))
    stopifnot(!is.null(names(pars_true)),all(names(pars_true)==pars))

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
  im_loss_mean <- signif(mean(im_loss_vals),digits)
  im_loss_sd <- signif(stats::sd(im_loss_vals),digits)
  nls_loss_mean <- signif(mean(nls_loss_vals),digits)
  nls_loss_sd <- signif(stats::sd(nls_loss_vals),digits)

  summary$loss <- data.frame(im_mean=im_loss_mean,im_sd=im_loss_sd,
                             nls_mean=nls_loss_mean,nls_sd=nls_loss_sd)

  im_est_mean <- signif(apply(im_est,2,mean),digits)
  im_est_sd <- signif(apply(im_est,2,stats::sd),digits)
  names(im_est_mean) <- pars

  nls_est_mean <- signif(apply(nls_est,2,mean),digits)
  nls_est_sd <- signif(apply(nls_est,2,stats::sd),digits)
  names(nls_est_mean) <- pars

  if(is.null(pars_true)) {
    summary$est <- data.frame(im_mean=im_est_mean, im_sd=im_est_sd,
                              nls_mean=nls_est_mean, nls_sd=nls_est_sd)
  }
  else {
    im_est_bias <- signif(im_est_mean - pars_true,digits)
    nls_est_bias <- signif(nls_est_mean - pars_true,digits)
    im_est_err <- im_est-matrix(rep(pars_true,simode_num),nrow=simode_num,byrow=T)
    im_est_rmse <- signif(sqrt(apply(im_est_err^2,2,mean)),digits)
    nls_est_err <- nls_est-matrix(rep(pars_true,simode_num),nrow=simode_num,byrow=T)
    nls_est_rmse <- signif(sqrt(apply(nls_est_err^2,2,mean)),digits)
    summary$est <- data.frame(true=signif(pars_true, digits),
                              im_mean=im_est_mean, im_sd=im_est_sd,
                              im_bias=im_est_bias, im_rmse = im_est_rmse,
                              nls_mean=nls_est_mean, nls_sd=nls_est_sd,
                              nls_bias=nls_est_bias, nls_rmse = nls_est_rmse)
  }



  class(summary) <- "summary.list.simode"
  return(summary)
}


## @rdname plot.simode
#' Plot the fit/estimates of a \code{list.simode} object
#'
#' Plot the fit or parameter estimates obtained from a call to \code{simode}
#' with \code{mc_sets>1}. Plots the mean and standard deviation obtained from
#' the multiple fits.
#'
#' @param x \code{simode} object returned by a call to \code{\link{simode}}.
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
           plot_im_smooth=F, legend=F, mfrow=par('mfrow'),
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

    if(type == 'fit')
      plot_list_simode_fit(x=x,show=show,which=which,pars_true=pars_true,
                           time=time,plot_im_smooth=plot_im_smooth,legend=legend,
                           mfrow=mfrow, cols=cols, ...)
    else #if(type == 'est')
      plot_list_simode_est(x=x,show=show,which=which,pars_true=pars_true,
                           legend=legend, cols=cols, ...)
}


plot_list_simode_est <-
        function(x, show=c('nls','im','both'),
                 which=NULL, pars_true=NULL, legend=F,
                 cols=list(nls_fit="blue",im_fit="green", true="black"), ...) {

  show <- match.arg(show)
  plot_nls_est <- (show=='nls' || show=='both')
  plot_im_est <- (show=='im' || show=='both')

  pars <- which
  if(is.null(pars))
    pars <- x[[1]]$pars

  if(!is.null(pars_true)) {
    stopifnot(is.numeric(pars_true))
    stopifnot(all(pars %in% names(pars_true)))
  }

  lp <- length(pars)
  pind <- which(as.vector(x[[1]]$pars) %in% pars)

  sum <- summary(x)
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

  y_min <- Inf
  y_max <- -Inf
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
  if(!is.null(pars_true)) {
    y_min <- min(c(y_min,pars_true[pars]),na.rm=T)
    y_max <- max(c(y_max,pars_true[pars]),na.rm=T)
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
    points(x=1:lp, pars_true[pars], col=cols[['true']], pch=4, xaxt="n")
    legend_text <- c(legend_text,'true')
    legend_cols <- c(legend_cols,cols[['true']])
    legend_lty <- c(legend_lty,0)
    legend_pch <- c(legend_pch,4)
  }

  if(plot_im_est && !is.null(im_est_mean)) {
    x_vals <- 1:lp
    if(plot_nls_est)
      x_vals <- x_vals-0.1
    points(x=x_vals, im_est_mean[pars], col=cols[['im_fit']], xaxt="n")
    if(length(x)>1) {
      errorbar(x=x_vals, im_est_mean[pars], yerr=im_est_sd[pars],
               bar.col=cols[['im_fit']], grid=F, add=T, xaxt="n")
      legend_text <- c(legend_text,'im_mean','im_sd')
      legend_cols <- c(legend_cols,cols[['im_fit']],cols[['im_fit']])
      legend_lty <- c(legend_lty,0,1)
      legend_pch <- c(legend_pch,1,NA)
    }
    else {
      legend_text <- c(legend_text,'im_est')
      legend_cols <- c(legend_cols,cols[['im_fit']])
      legend_lty <- c(legend_lty,0)
      legend_pch <- c(legend_pch,1)
    }
  }

  if(plot_nls_est && !is.null(nls_est_mean)) {
    x_vals <- 1:lp
    if(plot_im_est)
      x_vals <- x_vals+0.1
    points(x=x_vals, nls_est_mean[pars], col=cols[['nls_fit']], xaxt="n")
    if(length(x) > 1) {
      errorbar(x=x_vals, nls_est_mean[pars], yerr=nls_est_sd[pars],
               bar.col=cols[['nls_fit']], grid=F, add=T, xaxt="n")
      legend_text <- c(legend_text,'nls_mean','nls_sd')
      legend_cols <- c(legend_cols,cols[['nls_fit']],cols[['nls_fit']])
      legend_lty <- c(legend_lty,0,1)
      legend_pch <- c(legend_pch,1,NA)
    }
    else {
      legend_text <- c(legend_text,'nls_est')
      legend_cols <- c(legend_cols,cols[['nls_fit']])
      legend_lty <- c(legend_lty,0)
      legend_pch <- c(legend_pch,1)
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
           plot_im_smooth=F, legend=F, mfrow=par('mfrow'),
           cols=list(nls_fit="blue",im_fit="green", true="black",
                     obs="red", im_smooth="magenta"), ...) {


  show <- match.arg(show)
  plot_nls_fit <- (show=='nls' || show=='both')
  plot_im_fit <- (show=='im' || show=='both')

  vars <- which
  if(is.null(vars))
    vars <- names(x[[1]]$equations)

  if(!is.null(pars_true)) {
    stopifnot(is.numeric(pars_true))
    stopifnot(all(setdiff(x[[1]]$pars,x[[1]]$likelihood_pars) %in% names(pars_true)))
  }

  if(!is.list(x[[1]]$time)) {
    if(is.null(time))
      time <- x[[1]]$time
    time_obs <- rep(list(x[[1]]$time),length(x[[1]]$obs))
    names(time_obs) <- names(x[[1]]$obs)
  }
  else {
    if(is.null(time))
      time <- sort(unique(unlist(x[[1]]$time)))
    time_obs <- x[[1]]$time
  }

  y_min <- unlist(lapply(vars, function(i) min(x[[1]]$obs[[i]])))
  y_max <- unlist(lapply(vars, function(i) max(x[[1]]$obs[[i]])))
  names(y_min) <- names(x[[1]]$obs[vars])
  names(y_max) <- names(x[[1]]$obs[vars])

  legend_text <- c('obs')
  legend_cols <- cols[['obs']]

  x0.na <- names(x[[1]]$x0[is.na(x[[1]]$x0)])

  if(!is.null(pars_true)) {
    x0_true <- x[[1]]$x0
    x0_true[x0.na] <- pars_true[x0.na]
    pars_true <- pars_true[setdiff(names(pars_true),x0.na)]
    model_true <- solve_ode(x[[1]]$equations, pars_true, x0_true, time)
    legend_text <- c(legend_text, 'true')
    legend_cols <- c(legend_cols, cols[['true']])
    y_min1 <- apply(model_true[,vars], 2, function(col) min(col))
    y_max1 <- apply(model_true[,vars], 2, function(col) max(col))
    y_min <- pmin(y_min,y_min1)
    y_max <- pmax(y_max,y_max1)
  }
  if(plot_nls_fit) {
    model_nls_est <- list()
    for(i in 1:length(x)) {
      if(!is.null(x[[i]]$nls_pars_est)) {
        x0_nls <- x[[i]]$x0
        x0_nls[x0.na] <- x[[i]]$nls_pars_est[x0.na]
        pars_est_nls <- x[[i]]$nls_pars_est[setdiff(names(x[[i]]$nls_pars_est),x0.na)]
        model_nls_est[[i]] <- solve_ode(x[[i]]$equations, pars_est_nls, x0_nls, time)
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
      if(length(model_nls_est)>1)
        model_nls_est_sd[,var] <- apply(nls_est_mat,1,stats::sd)
    }
    legend_text <- c(legend_text, 'nls_fit')
    legend_cols <- c(legend_cols, cols[['nls_fit']])
    y_min1 <- apply(model_nls_est_mean[,vars]-model_nls_est_sd[,var], 2, function(col) min(col))
    y_max1 <- apply(model_nls_est_mean[,vars]+model_nls_est_sd[,var], 2, function(col) max(col))
    y_min <- pmin(y_min,y_min1)
    y_max <- pmax(y_max,y_max1)
  }
  if(plot_im_fit) {
    model_im_est <- list()
    for(i in 1:length(x)) {
      if(!is.null(x[[i]]$im_pars_est)) {
        x0_im <- x[[i]]$x0
        x0_im[x0.na] <- x[[i]]$im_pars_est[x0.na]
        pars_est_im <- x[[i]]$im_pars_est[setdiff(names(x[[i]]$im_pars_est),x0.na)]
        model_im_est[[i]] <- solve_ode(x[[i]]$equations, pars_est_im, x0_im, time)
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
      if(length(model_im_est)>1)
        model_im_est_sd[,var] <- apply(im_est_mat,1,stats::sd)
    }
    legend_text <- c(legend_text, 'im_fit')
    legend_cols <- c(legend_cols, cols[['im_fit']])
    y_min1 <- apply(model_im_est_mean[,vars]-model_im_est_sd[,var], 2, function(col) min(col))
    y_max1 <- apply(model_im_est_mean[,vars]+model_im_est_sd[,var], 2, function(col) max(col))
    y_min <- pmin(y_min,y_min1)
    y_max <- pmax(y_max,y_max1)
  }
  if(plot_im_smooth ) {
    im_time <- x[[1]]$im_smooth$time
    model_im_smooth_mean <- matrix(0,nrow=length(im_time), ncol=length(vars))
    model_im_smooth_sd <- matrix(0,nrow=length(im_time), ncol=length(vars))
    colnames(model_im_smooth_mean) <- vars
    colnames(model_im_smooth_sd) <- vars
    for(var in vars) {
      im_smooth_mat <- matrix(0,nrow=length(im_time), ncol=length(x))
      for(i in 1:length(x)) {
        im_smooth_mat[,i] <- x[[i]]$im_smooth$val[,var]
      }
      model_im_smooth_mean[,var] <- apply(im_smooth_mat,1,mean)
      if(length(x)>1)
        model_im_smooth_sd[,var] <- apply(im_smooth_mat,1,stats::sd)
    }
    legend_text <- c(legend_text, 'im_smooth')
    legend_cols <- c(legend_cols, cols[['im_smooth']])
  }

  if(legend)
    old.par <- par(mar=c(4,4,2,6), mfrow=mfrow, xpd = T)
  else
    old.par <- par(mar=c(4,4,2,2), mfrow=mfrow)

  on.exit(par(old.par))

  for(var in vars) {

    plot(-100,-100,xlab='time',ylab=var,
         xlim=c(time[1],time[length(time)]),
         ylim=c(y_min[var],y_max[var]), ...)

    for(i in 1:length(x)) {
      points(time_obs[[var]],x[[i]]$obs[[var]],col=cols[['obs']])
    }

    if(!is.null(pars_true)) {
      lines(time,model_true[,var],col=cols[['true']])
    }
    if(plot_im_smooth) {
      lines(im_time,model_im_smooth_mean[,var],col=cols[['im_smooth']])
      if(length(x)>1)
        errorbar(x=im_time, model_im_smooth_mean[,var], yerr=model_im_smooth_sd[,var],
                 bar.col=cols[['im_smooth']], grid=F, add=T, xaxt="n")
    }
    if(plot_im_fit) {
      lines(time,model_im_est_mean[,var],col=cols[['im_fit']])
      if(length(x)>1)
        errorbar(x=time, model_im_est_mean[,var], yerr=model_im_est_sd[,var],
                 bar.col=cols[['im_fit']], grid=F, add=T, xaxt="n")
    }
    if(plot_nls_fit) {
      lines(time,model_nls_est_mean[,var],col=cols[['nls_fit']])
      if(length(x)>1)
        errorbar(x=time, model_nls_est_mean[,var], yerr=model_nls_est_sd[,var],
                   bar.col=cols[['nls_fit']], grid=F, add=T, xaxt="n")
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
