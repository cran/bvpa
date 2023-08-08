#'  Estimation of 3-parameter Block-Basu Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Parameters estimation of 3-parameter BBBVPA distribution.
#'
#' @param I bivariate observations.
#' @param a0.int initial choice of \eqn{\alpha_0}.
#' @param a1.int initial choice of \eqn{\alpha_1}.
#' @param a2.int initial choice of \eqn{\alpha_2}.
#' @param tol.est convergence tolerance, \code{0.0001} (default).
#' @param MxIter.no maximum number of iterations, \code{2000} (default).
#' @param condition convergence criterion, \code{"log.L"} (default) and \code{"p.logL"}.
#'
#' @return  Object of class "\code{bbbvpa3}", a list consisting of
#' \item{alpha0, alpha1, alpha2, iter.no}{estimates  of parameters and number of iteration.}
#' \item{data }{the supplied data \code{I}.}
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' dat <- rbb.bvpa(500, 0, 0, 1.0, 1.0, 2.0, 0.4, 0.5)
#' estimates3(dat, 2.4, 0.3, 0.6)[-5]
#'
#' @export estimates3
estimates3 <-
  function(I, a0.int, a1.int, a2.int, tol.est=0.00001, MxIter.no=2000, condition="log.L"){
  n <- nrow(I)
  if(condition=="p.logL"){
    Qold <- pseu.logL(I, 0, 0, 1, 1, a0.int, a1.int, a2.int)
  } else if(condition=="log.L"){
    Qold <- logL(I, 0, 0, 1, 1, a0.int, a1.int, a2.int)$logLik
  }else stop("'condition' must be 'p.logL' or 'log.L'")

  t = 0
  j = 0
  while(j < 1)
  {
    I1 <- I

    I11 <- cbind(I1[I1[, 1] < I1[,2],1], I1[I1[, 1] < I1[,2],2])
    I12 <- cbind(I1[I1[, 1] > I1[,2],1], I1[I1[, 1] > I1[,2],2])

    n1 <- length(I11[,1])
    n2 <- length(I12[,1])

    u1 <- a0.int/(a0.int + a2.int)
    u2 <- a2.int/(a0.int + a2.int)
    w1 <- a0.int/(a0.int + a1.int)
    w2 <- a1.int/(a0.int + a1.int)

    n0 <- (n1 + n2)*(a0.int/(a1.int + a2.int))
    a <- 1/(a0.int + a1.int + a2.int)

    alpfin0 <- (n0 + u1*n1 + w1*n2)/(n0*a + sum(log(1 + I12[,1])) + sum(log(1 + I11[,2])))
    alpfin1 <- (n1 + w2*n2)/(n0*a + sum(log(1 + I11[,1])) + sum(log(1 + I12[,1])))
    alpfin2 <- (n2 + u2*n1)/(n0*a + sum(log(1 + I11[,2])) + sum(log(1 + I12[,2])))

    if(condition=="p.logL"){
      QQ <- pseu.logL(I, 0, 0, 1, 1, alpfin0, alpfin1, alpfin2)
    } else if(condition=="log.L"){
      QQ <- logL(I, 0, 0, 1, 1, alpfin0, alpfin1, alpfin2)$logLik
    }else stop("'condition' must be 'p.logL' or 'log.L'")

    delQ <- QQ - Qold
    t <- t+1
    if((abs(delQ/Qold) < tol.est)||(t >= MxIter.no)){
      break
    }
    a0.int <- alpfin0
    a1.int <- alpfin1
    a2.int <- alpfin2
    Qold <- QQ
  }
  res <- list(alpha0=alpfin0, alpha1=alpfin1, alpha2=alpfin2, iter.no=t)
  res$data <- I
  class(res) <- c("bbbvpa3")
  return(res)
}
