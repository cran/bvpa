#' Observed Fisher information based confidence interval of 3-parameter Block-Basu
#' Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Observed Fisher information based confidence interval of 3-parameter BBBVPA distribution.
#'
#' @param object \code{"bbbvpa"} class object.
#' @param conf.lev confidence level, \eqn{0.95} (default).
#' @param tol convergence tolerance for confidence intervals, \code{0.0001} (default).
#' @return  A matrix of lower and upper confidence interval limits (in the first and second column respectively).
#' The matrix rows are labeled by the parameter names (if any) and columns by the corresponding
#' distribution quantiles.
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' dat <- rbb.bvpa(500, 0, 0, 1.0, 1.0, 2.0, 0.4, 0.5)
#' conf.intv3(estimates3(dat, 2.4, 0.3, 0.6))
#'
#' @importFrom stats qnorm
#'
#' @export conf.intv3
conf.intv3 <-
  function(object, conf.lev=0.95, tol=0.0001){
  I <- object$data
  alpha0 <- object$alpha0
  alpha1 <- object$alpha1
  alpha2 <- object$alpha2

  no.obs <- nrow(I)
  u.pctl <- (1+conf.lev)/2
  c.v <- qnorm(u.pctl) # critical value/ cut-off value

  I1 <- I

  I11 <- cbind(I1[I1[, 1] < I1[,2],1], I1[I1[, 1] < I1[,2],2])
  I12 <- cbind(I1[I1[, 1] > I1[,2],1], I1[I1[, 1] > I1[,2],2])

  n1 <- length(I11[,1])
  n2 <- length(I12[,1])

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

  res <- cbind(c(ci_alp0L, ci_alp1L, ci_alp2L), c(ci_alp0U, ci_alp1U, ci_alp2U))
  rownames(res) <- c("alpha0", "alpha1", "alpha2")
  colnames(res) <- c("lower.intv", "upper.intv")
  return(res)
}
