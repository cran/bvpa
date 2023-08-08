#' Observed Fisher information based confidence interval of Block-Basu Bivariate
#' Pareto (BBBVPA) distribution
#'
#' @description
#' Observed Fisher information based confidence interval of Bivariate BBBVPA
#' distribution.
#'
#' @param object \code{"bbbvpa"} class object.
#' @param conf.lev confidence level, \eqn{0.95} (default).
#' @param tol convergence tolerance for confidence intervals, \code{0.0001} (default).
#' @param intv.m1 interval related to confidence interval of \eqn{\mu_1}, \code{c(0,2)} (default).
#' @param intv.m2 interval related to confidence interval of \eqn{\mu_1}, \code{c(0,2)} (default).
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
#' @importFrom stats qnorm
#' @importFrom numDeriv hessian
#'
#'
#' @export conf.intv
conf.intv <-
  function(object, conf.lev=0.95, tol=0.0001, intv.m1=c(0,2), intv.m2=c(0,2)){
  I <- object$data
  mu1 <- object$mu1
  mu2 <- object$mu2
  sigma1 <- object$sigma1
  sigma2 <- object$sigma2
  alpha0 <- object$alpha0
  alpha1 <- object$alpha1
  alpha2 <- object$alpha2

  no.obs <- nrow(I)
  l.pctl <- (1-conf.lev)/2
  u.pctl <- (1+conf.lev)/2
  c.v <- qnorm(u.pctl) # critical value/ cut-off value

  ci.mu1L <- mu1 - sigma1*uniroot(pctl.fun, interval=intv.m1, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 1, pct = l.pctl, tol = tol)[[1]]
  ci.mu1U <- mu1 - sigma1*uniroot(pctl.fun, interval=intv.m1, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 1, pct = u.pctl, tol = tol)[[1]]

  ci.mu2L <- mu2 - sigma2*uniroot(pctl.fun, interval=intv.m2, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 2, pct = l.pctl, tol = tol)[[1]]
  ci.mu2U <- mu2 - sigma2*uniroot(pctl.fun, interval=intv.m2, n = no.obs, a0 = alpha0, a1 = alpha1, a2 = alpha2, select = 2, pct = u.pctl, tol = tol)[[1]]

  OFI1 <- -hessian(mLf1, x=sigma1, mu1=mu1, a0=alpha0, a1=alpha1, a2=alpha2, I = I)
  OFI2 <- -hessian(mLf2, x=sigma2, mu2=mu2, a0=alpha0, a1=alpha1, a2=alpha2, I = I)

  ci.sig1L <- sigma1 - c.v*sqrt(1/OFI1)
  ci.sig1U <- sigma1 + c.v*sqrt(1/OFI1)

  ci.sig2L <- sigma2 - c.v*sqrt(1/OFI2)
  ci.sig2U <- sigma2 + c.v*sqrt(1/OFI2)

  Iv1 <- (I[, 1] - mu1)/sigma1
  Iv2 <- (I[, 2] - mu2)/sigma2
  I1 <- cbind(Iv1, Iv2)

  I11 <- cbind(I1[I1[, 1] < I1[,2],1], I1[I1[, 1] < I1[,2],2])
  I12 <- cbind(I1[I1[, 1] > I1[,2],1], I1[I1[, 1] > I1[,2],2])

  # n1 <- length(I11[,1])
  # n2 <- length(I12[,1])

  e.n1 <- length(I11[,1])
  e.n2 <- length(I12[,1])

  n1 <- (e.n1+e.n2)*alpha1/(alpha1+alpha2)
  n2 <- (e.n1+e.n2)*alpha2/(alpha1+alpha2)

  u1 <- alpha0/(alpha0 + alpha2)
  u2 <- alpha2/(alpha0 + alpha2)
  w1 <- alpha0/(alpha0 + alpha1)
  w2 <- alpha1/(alpha0 + alpha1)

  n0 <- (n1 + n2)*(alpha0/(alpha1 + alpha2))
  a <- 1/(alpha0 + alpha1 + alpha2)

  B <- matrix(c((n0 + u1*n1 + w1*n2)/alpha0^2, 0, 0, 0, (n1 + w2*n2)/alpha1^2, 0, 0, 0, (n2 + u2*n1)/alpha2^2), nrow=3, ncol=3, byrow = TRUE)

  S <- matrix(c(-(n0*a + sum(log(1 + I12[,1])) + sum(log(1 + I11[,2]))) + (n0 + u1*n1 + w1*n2)/alpha0,
                -(n0*a + sum(log(1 + I11[,1])) + sum(log(1 + I12[,1]))) + (n1 + w2*n2)/alpha1,
                -(n0*a + sum(log(1 + I11[,2])) + sum(log(1 + I12[,2]))) + (n2 + u2*n1)/alpha2),
              nrow=3, ncol=1, byrow = TRUE)

  OFI <- B-S%*%t(S)
  diag.elem <- diag(solve(OFI))

  ci_alp0L <- alpha0 - c.v*sqrt(diag.elem[1])
  ci_alp0U <- alpha0 + c.v*sqrt(diag.elem[1])

  ci_alp1L <- alpha1 - c.v*sqrt(diag.elem[2])
  ci_alp1U <- alpha1 + c.v*sqrt(diag.elem[2])

  ci_alp2L <- alpha2 - c.v*sqrt(diag.elem[3])
  ci_alp2U <- alpha2 + c.v*sqrt(diag.elem[3])

  res <- cbind(c(ci.mu1L, ci.mu2L, ci.sig1L, ci.sig2L, ci_alp0L, ci_alp1L, ci_alp2L), c(ci.mu1U, ci.mu2U, ci.sig1U, ci.sig2U, ci_alp0U, ci_alp1U, ci_alp2U))
  rownames(res) <- c("mu1", "mu2", "sigma1", "sigma2", "alpha0", "alpha1", "alpha2")
  colnames(res) <- c("lower.intv", "upper.intv")
  return(res)
}
