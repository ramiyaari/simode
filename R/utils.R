

setdiff1 <- function(A, B) {
  res <- setdiff(A,B)
  if(length(res)==0)
    res <- NULL
  return (res)
}

setdiff2 <- function(A, B) {
  if(is.null(names(A)))
    return (A)
  res <- A[setdiff(names(A),B)]
  if(length(res)==0)
    res <- NULL
  return (res)
}


# Fix parameter values in the equations
#
# Substitues the parameter names with the parameter values in the
# given equations.
# @param equations the equations describing the ode system.
# Each equation must be named according to the variable it describes.
# @param pars the parameter values. Each parameter must be named
# according to the parameter name as it appears in the equations.
# @return The equations with the fixed parameter values.
fix_pars <- function(equations, pars) {

  vars <- names(equations)
  par_names <- names(pars)
  for(i in 1:length(pars))
    equations <- gsub(paste0('\\<',par_names[i],'\\>'),pars[i],equations)

  names(equations) <- vars
  return (equations)
}

check_pars_presence_in_equations <- function(equations, pars) {

  result <- T
  for(i in 1:length(pars)) {
    found <- grep(pars[i],equations)
    if(is.null(found) || length(found)==0) {
      cat(noquote(paste0('Parameter [', pars[i],
                         '] is not used in equations and is not in likelihood_pars.\n')))
      result <- F
    }
  }
  return (result)
}


# Check linearity of parameters
#
# Check if the given parameters are linear in the equations.
# This function is called within the call to \link{simode} to
# verify that the given linear parameters are indeed linear in the
# given system.
#
# @param equations The equations describing the ode system. See \link{simode}.
# @param pars The parameter names which should be linear.
# @return TRUE if the parameters are linear in the system and
# FALSE otherwise.
# @export
#
check_linearity <- function(equations, pars) {

  linear <- T
  eq_names <- names(equations)
  equations <- sapply(1:length(equations), function(i) parse(text=equations[i]))
  for(k in 1:length(equations)) {
    for(i in 1:length(pars)) {
      eq.deriv <- D(equations[k], pars[i])
      eq.deriv2 <- deparse(D(eq.deriv, pars[i]))
      if(eq.deriv2[1] != "0")
      {
        cat(noquote(paste0('Problem in eq.', k, ' [', eq_names[k],
                           '] - parameter [', pars[i],
                           '] should be set as non-linear\n')))

        linear <- F
      }
    }
  }
  if(!linear)
    return (F)

  if(length(pars) > 1) {
    for(k in 1:length(equations)) {
      for(i in 1:(length(pars)-1)) {
        eq.deriv <- D(equations[k], pars[i])
        for(j in (i+1):length(pars)) {
          eq.deriv2 <- deparse(D(eq.deriv, pars[j]))
          if(eq.deriv2[1] != "0")
          {
            cat(noquote(paste0('Problem in eq.', k, ' [', eq_names[k],
                               '] - parameter [', pars[i], '] or [', pars[j],
                               '] should be set as non-linear\n')))
            linear <- F
          }
        }
      }
    }
  }

  return (linear)
}


check_obs_validity <- function(obs,time,vars) {

  stopifnot(length(vars)<=length(obs),all(vars %in% names(obs)))
  stopifnot(is.numeric(time) || (is.list(time) && length(time)==length(obs)))
  for(i in 1:length(obs)) {
    if(is.numeric(time))
      stopifnot(length(time)==length(obs[[i]]))
    else
      stopifnot(is.numeric(time[[i]]),length(time[[i]])==length(obs[[i]]))
  }
}

# Creates a mapping of the parameters and the equations (variables)
# they affect.
# @param pars names of the parameters appearing in the equations.
# @param equations the equations describing the ode system.
# Each equation must have a name of the variable it describes.
# @param named whether the mapping should be based on the equation names or indexes.
# @return A mapping of the given parameters to the variables they affect
# @export
#
pars2vars_mapping <- function(pars, equations, named=F) {

  if(pracma::isempty(pars))
    return (c())

  vars <- names(equations)
  d <- length(vars)
  p <- length(pars)

  vars2vars <- list()
  for(i in 1:d) {
    vars1 <-
      vars[which(unlist(lapply(1:d,
           function(j) regexpr(pattern=paste0('\\<',vars[i],'\\>'),text=equations[j])!=-1)))]
    vars2vars[[i]] <- vars1
  }
  names(vars2vars) <- vars

  pars2vars <- list()
  for(i in 1:p) {
    if(pars[i] %in% vars)
      vars1 <- vars2vars[[pars[i]]]
    else
      vars1 <-
        vars[which(unlist(lapply(1:d,
        function(j) regexpr(pattern=paste0('\\<',pars[i],'\\>'),text=equations[j])!=-1)))]
    pars2vars[[i]] <- vars1
  }
  names(pars2vars) <- pars
  for(i in 1:p) {
    for(var in pars2vars[[i]]) {
      pars2vars <- add_linked_vars(pars[i],var,pars2vars,vars2vars)
    }
  }
  if(!named) {
    for(i in 1:p)
      pars2vars[[i]] <- which(vars %in% pars2vars[[i]])
  }
  else
  {
    for(i in 1:p)
      pars2vars[[i]] <- vars[which(vars %in% pars2vars[[i]])]
  }
  return (pars2vars)
}

add_linked_vars <- function(par, var, pars2vars, vars2vars)
{
  new_vars <- setdiff(vars2vars[[var]],pars2vars[[par]])
  if(!pracma::isempty(new_vars)) {
    pars2vars[[par]] <- c(pars2vars[[par]],new_vars)
    for(var2 in new_vars) {
      pars2vars <- add_linked_vars(par,var2,pars2vars,vars2vars)
    }
  }
  return (pars2vars)
}


#' Ordinary differential equations solver
#'
#' A wrapper for the \link[deSolve]{ode} function that solves a system of
#' ordinary differential equations described using symbolic equations.
#' @param equations The equations describing the ODE system. See \link{simode}.
#' @param pars The parameter values. Named according to their names in \code{equations}.
#' @param x0 The initial conditions. Named accroding to the names of \code{equations}.
#' @param time The time points for which the ODE variables' values will be computed.
#' @param xvars External observations of time-dependant variables refered to in \code{equations}
#' @param ... Additional argument(s) for methods.
#' @return A matrix whose first column contains the given time points
#' and subsequent columns hold the computed ODE equations' values at these time points.
#' @importFrom deSolve ode
#' @export
#'
solve_ode <- function(equations, pars, x0, time, xvars=NULL, ...)
{
  stopifnot(length(equations)==length(x0),
            !is.null(names(equations)),
            !is.null(names(x0)),
            all(names(x0)==names(equations)))

  equations <- fix_pars(equations, pars)
  if(!is.null(xvars) && length(xvars)>0) {
    stopifnot(is.list(xvars))
    for(j in 1:length(xvars)) {
      stopifnot(!is.null(names(xvars[j])))
      xvar_name <- names(xvars[j])
      equations <- gsub(paste0('\\<',xvar_name,'\\>'),
                        paste0(xvar_name,"[which.min(abs(times2-t))]"),equations)
    }
  }
  equations <- lapply(1:length(equations), function(i) parse(text=equations[i]))
  names(equations) <- names(x0)

  arg_list <- list(...)
  trace <- arg_list$trace
  arg_list$trace <- NULL
  if(is.null(trace))
    trace <- 3

  out <- NULL
  tryCatch(
  {
      if(trace<3){
        tmpfile <- file.path(tempdir(),'simode-ode.log')
        sink(file=tmpfile,append=F)
      }
      args <- c(list(y=x0, times=time, func=step_ode_model,
                     parms=pars, equations=equations, xvars=xvars, times2=time),arg_list)
      out <- do.call("ode", args)
      if(trace<3)
        sink(NULL)
    },
    warning = function(w) {
      if(trace<3)
        sink(NULL)
      if(trace>2)
        print(w)
    },
    error = function(e) {
      if(trace<3)
        sink(NULL)
      if(trace>1)
        print(e)
    }
  )
  return (out)
}

#' Ordinary differential equations solver using a \code{simode} object
#'
#' A wrapper for the \link[deSolve]{ode} function that solves a system of
#' ordinary differential equations described using symbolic equations.
#' @param x A simode object returned from a call to \code{\link{simode}}.
#' @param type Which solution to generate ('both'\\'nls'\\'im').
#' @return A matrix whose first column contains the given time points
#' and subsequent columns hold the computed ODE equations' values at these time points.
#' @export
#'
solve_ode2 <- function(x,type=c("both","im","nls")) {
  stopifnot(class(x)=='simode')
  type <- match.arg(type)
  x0.na <- names(x$x0[is.na(x$x0)])
  xvars <- setdiff(names(x$obs),names(x$equations))
  time <- x$time
  if(is.list(time)) {
    time <- sort(unique(unlist(time)))
  }
  solution <- list()
  nls_pars_est <- x$nls_pars_est
  if(type!='im' && !is.null(nls_pars_est)) {
    if(!is.null(x$scale_pars)) {
      args <- c(list(pars=nls_pars_est), x$extra_args)
      nls_pars_est <- do.call(x$scale_pars, args)
    }
    x0 <- x$x0
    x0[x0.na] <- nls_pars_est[x0.na]
    nls_pars_est <- nls_pars_est[setdiff(names(nls_pars_est),x0.na)]
    nls_solution <- solve_ode(x$equations, nls_pars_est, x0, time, x$obs[xvars])
    solution$nls <- nls_solution
  }
  im_pars_est <- x$im_pars_est
  if(type!='nls' && !is.null(im_pars_est)) {
    if(!is.null(x$scale_pars)) {
      args <- c(list(pars=im_pars_est), x$extra_args)
      im_pars_est <- do.call(x$scale_pars, args)
    }
    x0 <- x$x0
    x0[x0.na] <- im_pars_est[x0.na]
    im_pars_est <- im_pars_est[setdiff(names(im_pars_est),x0.na)]
    im_solution <- solve_ode(x$equations, im_pars_est, x0, x$time, x$obs[xvars])
    solution$im <- im_solution
  }
  return (solution)
}


step_ode_model <- function(t, states, pars, equations, xvars, times2)
{
  with(as.list(c(states,xvars)), {
    delta_states <- unlist(lapply(1:length(equations),function(i) eval(equations[[i]])))
    return (list(delta_states))
  })
}

smooth_obs <- function(im_smoothing=c('splines','kernel','none'),obs,time,t,bw_factor) {

  if(im_smoothing=='splines') {
    fit <- smooth.spline(time,obs,cv=F)
    sobs <- predict(fit,x=t)$y
  }
  else if(im_smoothing=='kernel') {
    N <- length(t)
    n <- length(time)
    b <- max(1,bw_factor)*max(diff(time))
    sobs <- c()
    for (k in 1:N) {
      ker <- fun_kernel(time, t[k], n, b)
      U <- calc_U(time, t[k], n, b)
      B <- calc_B(n, b, U, ker)
      W <- calc_W(n, b, U, B, ker)
      sobs[k] <- sum(obs*W)
    }
  }
  else {
    sobs <- interp1(time,obs,xi=t,method='linear') #'spline')
  }
  return (sobs)
}

############ Kernel related functions #########

fun_kernel <- function(time, t, n, b) {
  s <- (time - t)/b
  ker <- 0.75*(1-s^2)*(abs(s)<=1)
  return (ker)
}

calc_U <- function(time, t, n, b){
  U <- matrix(0, nrow=2, ncol = n)
  U[1,] <- 1
  U[2,] <- (time - t)/b
  return(U)
}

calc_B <- function(n, b, U, ker) {
  B <- matrix(0, nrow = 2, ncol = 2)
  for(i in 1:n){
    B <- B + U[,i]%*%t(U[,i])*ker[i]
  }
  B <- B/(n*b)
  return(B)
}

calc_W <- function(n, b, U, B, ker) {
  B_inv <- solve(B)
  W <- apply(matrix(1:n),1, function(i) t(U[,i])%*%B_inv%*%c(1,0)*ker[i])/(n*b)
  return(W)
}

