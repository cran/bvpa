#' Bivariate Pareto Distribution
#'
#' @description
#' Implements the EM algorithm with one-step Gradient Descent method to estimate the
#' parameters of the Block-Basu bivariate Pareto distribution with location and scale.
#' We also found parametric bootstrap and asymptotic confidence intervals based on
#' the observed Fisher information of scale and shape parameters, and exact confidence
#' intervals for location parameters.
#'
#' @name bvpa-package
#' @aliases bvpa-package bvpa
#' @docType package
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @references Biplab Paul and Arabin Kumar Dey (2023). An EM algorithm for absolutely
#' continuous Marshall-Olkin bivariate Pareto distribution with location
#' and scale, Preprint.
#' @references E L Lehmann and George Casella (1998). Theory of Point Estimation,
#' Springer, New York, doi.org/10.1007/b98854.
#' @references Bradley Efron and R J Tibshirani (1994). An Introduction to the Bootstrap, CRC press,
#' New York, doi.org/10.1201/9780429246593.
#' @references A P Dempster, N M Laird and D B Rubin (1977). Maximum Likelihood from
#' Incomplete Data via the EM Algorithm, Journal of the royal statistical society:
#' series B (methodological),
#' www.jstor.org/stable/2984875.
#'
#' @keywords package
#'
#' @importFrom numDeriv grad hessian
#' @importFrom stats qnorm quantile runif uniroot

NULL



