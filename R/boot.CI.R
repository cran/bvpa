#' Parametric bootstrap confidence intervals of parameters of Block-Basu Bivariate
#' Pareto (BBBVPA) distribution
#'
#' @description
#' Parametric bootstrap confidence interval of parameters of BBBVPA distribution.
#'
#' @param data bivariate observations.
#' @param s1.int initial choice of \eqn{\sigma_1}.
#' @param s2.int initial choice of \eqn{\sigma_2}.
#' @param a0.int initial choice of \eqn{\alpha_0}.
#' @param a1.int initial choice of \eqn{\alpha_1}.
#' @param a2.int initial choice of \eqn{\alpha_2}.
#' @param conf.lev confidence level, defult \eqn{0.95}.
#' @param intv.m1 interval related to confidence interval of \eqn{\mu_1}, \code{c(0,2)} (default).
#' @param intv.m2 interval related to confidence interval of \eqn{\mu_1}, \code{c(0,2)} (default).
#' @param no.paboot number of bootstrap samples, \code{100} (default).
#' @param tol convergence tolerance for confidence interval of \eqn{\mu_1}.
#' and \eqn{\mu_2}, \code{0.0001} (default).
#' @param ... further arguments to pass to \code{estimates}.
#'
#' @return  A matrix of lower and upper confidence interval limits (in the first and second column respectively).
#' The matrix rows are labeled by the parameter names (if any) and columns by the corresponding
#' distribution quantiles.
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' # see the example of estimation
#'
#' @importFrom stats uniroot
#' @importFrom stats quantile
#'
#' @export param.boot
param.boot <-
  function(data, s1.int, s2.int, a0.int, a1.int, a2.int, conf.lev=0.95, intv.m1=c(0,2),
                       intv.m2=c(0,2), no.paboot=100, tol=0.0001, ...){
  no.obs <- nrow(data)
  res <- estimates(data, s1.int, s2.int, a0.int, a1.int, a2.int, ...)
  mu1 <- res$mu1
  mu2 <- res$mu2
  sigma1 <- res$sigma1
  sigma2 <- res$sigma2
  alpha0 <- res$alpha0
  alpha1 <- res$alpha1
  alpha2 <- res$alpha2

  l.pctl <- (1-conf.lev)/2
  u.pctl <- (1+conf.lev)/2

  ci.mu1L <- mu1 - sigma1*uniroot(pctl.fun, interval=intv.m1, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 1, pct = l.pctl, tol = tol)[[1]]
  ci.mu1U <- mu1 - sigma1*uniroot(pctl.fun, interval=intv.m1, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 1, pct = u.pctl, tol = tol)[[1]]

  ci.mu2L <- mu2 - sigma2*uniroot(pctl.fun, interval=intv.m2, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 2, pct = l.pctl, tol = tol)[[1]]
  ci.mu2U <- mu2 - sigma2*uniroot(pctl.fun, interval=intv.m2, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 2, pct = u.pctl, tol = tol)[[1]]

  ve1 <- ve2 <- ve3 <- ve4 <- ve5 <- ve6 <- ve7 <- rep(0, no.paboot)

  for (i in 1:no.paboot) {
    sim.dat <- rbb.bvpa(no.obs, mu1, mu2, sigma1, sigma2, alpha0, alpha1, alpha2)
    res1 <- estimates(sim.dat, s1.int, s2.int, a0.int, a1.int, a2.int, ...)

    ve1[i] <- res1$mu1
    ve2[i] <- res1$mu2
    ve3[i] <- res1$sigma1
    ve4[i] <- res1$sigma2
    ve5[i] <- res1$alpha0
    ve6[i] <- res1$alpha1
    ve7[i] <- res1$alpha2
    cat('\n', i, '\n')
  }
  CI_mu1 <- c(ci.mu1L, ci.mu1U)
  CI_mu2 <- c(ci.mu2L, ci.mu2U)
  CI_sig1 = c(quantile(ve3, l.pctl), quantile(ve3, u.pctl))
  CI_sig2 = c(quantile(ve4, l.pctl), quantile(ve4, u.pctl))
  CI_alp0 = c(quantile(ve5, l.pctl), quantile(ve5, u.pctl))
  CI_alp1 = c(quantile(ve6, l.pctl), quantile(ve6, u.pctl))
  CI_alp2 = c(quantile(ve7, l.pctl), quantile(ve7, u.pctl))

  rslt <- rbind(CI_mu1, CI_mu2, CI_sig1, CI_sig2, CI_alp0, CI_alp1, CI_alp2)
  rownames(rslt) <- c("mu1", "mu2", "sigma1", "sigma2", "alpha0", "alpha1", "alpha2")
  colnames(rslt) <- c("lower.intv", "upper.intv")
  return(rslt)
}
