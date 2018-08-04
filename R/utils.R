

setdiff1 <- function(A, B) {
  res <- setdiff(A,B)
  if(length(res)==0)
    res <- NULL
  return (res)
}

setdiff2 <- function(A, B) {
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

  return (linear)
}


check_obs_validity <- function(obs,time,vars) {

  stopifnot(length(vars)==length(obs),all(vars==names(obs)))
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
#' @param ... Additional argument(s) for methods.
#' @return A matrix whose first column contains the given time points
#' and subsequent columns hold the computed ODE equations' values at these time points.
#' @importFrom deSolve ode
#' @export
#'
solve_ode <- function(equations, pars, x0, time, ...)
{
  stopifnot(length(equations)==length(x0),
            !is.null(names(equations)),
            !is.null(names(x0)),
            all(names(x0)==names(equations)))

  equations <- fix_pars(equations, pars)
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
                     parms=pars, equations=equations),arg_list)
      out <- do.call("ode", args)
      #out <- ode(y=x0, times=time, func=step_ode_model,
      #           parms=pars, equations=equations, ...)
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

step_ode_model <- function(t, states, pars, equations)
{
  with(as.list(states), {
    delta_states <- unlist(lapply(1:length(equations),function(i) eval(equations[[i]])))
    return (list(delta_states))
  })
}


