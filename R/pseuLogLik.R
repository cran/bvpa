#' Pseudo log-likelihood function of Block-Basu Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Return the pseudo log likelihood value.
#'
#' @param I baivariate observations.
#' @param mu1 value of \eqn{\mu_1}.
#' @param mu2 value of \eqn{\mu_2}.
#' @param s1 value of \eqn{\sigma_1}.
#' @param s2 value of \eqn{\sigma_2}.
#' @param a0 value of \eqn{\alpha_0}.
#' @param a1 value of \eqn{\alpha_1}.
#' @param a2 value of \eqn{\alpha_2}.
#'
#' @return A scalar numeric, pseudo log likelihood of the model.
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' dat <- rbb.bvpa(500, 0.1, 0.1, 0.8, 0.8, 2.0, 0.4, 0.5)
#' pseu.logL(dat, 0.1, 0.1, 0.8, 0.8, 2.0, 0.4, 0.5)
#'
#' @export pseu.logL
pseu.logL <-
  function(I, mu1, mu2, s1, s2, a0, a1, a2){
  Iv1 <- (I[, 1] - mu1)/s1
  Iv2 <- (I[, 2] - mu2)/s2
  I1 <- cbind(Iv1, Iv2)

  I11 <- cbind(I1[I1[, 1] < I1[,2],1], I1[I1[, 1] < I1[,2],2])
  I12 <- cbind(I1[I1[, 1] > I1[,2],1], I1[I1[, 1] > I1[,2],2])

  # n1 <- length(I11[,1])
  # n2 <- length(I12[,1])

  e.n1 <- length(I11[,1])
  e.n2 <- length(I12[,1])
  # n <- (e.n1+e.n2)

  n1 <- (e.n1+e.n2)*a1/(a1+a2)
  n2 <- (e.n1+e.n2)*a2/(a1+a2)

  u1 <- a0/(a0 + a2)
  u2 <- a2/(a0 + a2)
  w1 <- a0/(a0 + a1)
  w2 <- a1/(a0 + a1)

  n0 <- (n1 + n2)*(a0/(a1 + a2))
  a <- 1/(a0 + a1 + a2)
  Q <- -a0*(n0*a + sum(log(1 + I12[,1])) + sum(log(1 + I11[,2]))) + (n0 + u1*n1 + w1*n2)*log(a0) - a1*(n0*a + sum(log(1 + I11[,1])) + sum(log(1 + I12[,1]))) + (n1 + w2*n2)*log(a1) - a2*(n0*a + sum(log(1 + I11[,2])) + sum(log(1 + I12[,2]))) + (n2 + u2*n1)*log(a2)
  return(Q)
}

