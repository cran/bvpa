#' Parametric bootstrap confidence intervals of parameters of 3-parameter Block-Basu
#' Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Parametric bootstrap confidence interval of parameters of 3-parameter BBBVPA distribution.
#'
#' @param data bivariate observations.
#' @param a0.int initial choice of \eqn{\alpha_0}.
#' @param a1.int initial choice of \eqn{\alpha_1}.
#' @param a2.int initial choice of \eqn{\alpha_2}.
#' @param conf.lev confidence level, defult 0.95.
#' @param no.paboot number of bootstrap samples, \code{100} (default).
#' @param ... further arguments to pass to \code{estimates3}.
#'
#' @return  A matrix of lower and upper confidence interval limits (in the first and second column respectively).
#' The matrix rows are labeled by the parameter names (if any) and columns by the corresponding
#' distribution quantiles.
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' dat <- rbb.bvpa(500, 0, 0, 1.0, 1.0, 2.0, 0.4, 0.5)
#' param.boot3(dat, 2.4, 0.3, 0.6)
#'
#' @importFrom stats quantile
#'
#' @export param.boot3
param.boot3 <- function(data, a0.int, a1.int, a2.int, conf.lev=0.95, no.paboot=100, ...){
  no.obs <- nrow(data)
  res <- estimates3(data, a0.int, a1.int, a2.int, ...)
  alpha0 <- res$alpha0
  alpha1 <- res$alpha1
  alpha2 <- res$alpha2

  l.pctl <- (1-conf.lev)/2
  u.pctl <- (1+conf.lev)/2

  ve1 <- ve2 <- ve3 <- rep(0, no.paboot)

  for (i in 1:no.paboot) {
    sim.dat <- rbb.bvpa(no.obs, 0, 0, 1, 1, alpha0, alpha1, alpha2)
    res1 <- estimates3(sim.dat, a0.int, a1.int, a2.int, ...)

    ve1[i] <- res1$alpha0
    ve2[i] <- res1$alpha1
    ve3[i] <- res1$alpha2
    cat('\n', i, '\n')
  }

  CI_alp0 = c(quantile(ve1, l.pctl), quantile(ve1, u.pctl))
  CI_alp1 = c(quantile(ve2, l.pctl), quantile(ve2, u.pctl))
  CI_alp2 = c(quantile(ve3, l.pctl), quantile(ve3, u.pctl))

  rslt <- rbind(CI_alp0, CI_alp1, CI_alp2)
  rownames(rslt) <- c("alpha0", "alpha1", "alpha2")
  colnames(rslt) <- c("lower.intv", "upper.intv")
  return(rslt)
}

