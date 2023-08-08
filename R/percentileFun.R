#' Survival functions of pivots of estimators of locations.
#'
#' @description
#' Survival functions of pivots of estimators of locations \eqn{\mu_1} and \eqn{\mu_2}.
#' These are required to calculate the critical value of confidence intervals for \eqn{\mu_1}
#' and \eqn{\mu_2}.
#'
#' @param z quantiles.
#' @param n number of observations.
#' @param a0 value of \eqn{\alpha_0}.
#' @param a1 value of \eqn{\alpha_1}.
#' @param a2 value of \eqn{\alpha_2}.
#' @param pct probabilities.
#' @param select Allows to select the function for different location parameters. a single model term to be selected for printing.
#'  e.g. if you just want the function for \eqn{\mu_1} set \eqn{select=1} (default).
#'
#' @return  return a function.
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' uniroot(pctl.fun, interval=c(0,2), n = 500, a0 = 2.0, a1 = 0.4, a2 = 0.5,
#'  pct = 0.025, tol = 0.0001)[[1]]
#'
#' @export pctl.fun
pctl.fun <-
  function(z, n, a0, a1, a2, pct, select=1)
{
  if (select==1){
    ((a0+a1+a2)/(a1+a2))*(1+z)^(-a0-a1)-(a0/(a1+a2))*(1+z)^(-a0-a1-a2) - (pct)^(1/n)
  }else if(select==2){
    ((a0+a1+a2)/(a1+a2))*(1+z)^(-a0-a2)-(a0/(a1+a2))*(1+z)^(-a0-a1-a2) - (pct)^(1/n)
  }else
    stop("'select' must be 1 or 2")
}
