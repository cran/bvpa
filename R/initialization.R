#' Initialization of Block-Basu Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Return initial choice parameters of BBBVPA distribution.
#'
#' @param data bivariate observations.
#' @param ini.run number of random initializations.
#' @param tol.ini convergence tolerance, \code{0.001} (default)..
#' @param proc different procedures, \code{"ML"} (default) and \code{"S.EM"}.
#' @param intv.s1 interval for random initialization of \eqn{\sigma_1}.
#' @param intv.s2 interval for random initialization of \eqn{\sigma_2}.
#' @param intv.a0 interval for random initialization of \eqn{\alpha_0}.
#' @param intv.a1 interval for random initialization of \eqn{\alpha_1}.
#' @param intv.a2 interval for random initialization of \eqn{\alpha_2}.
#' @param ... further arguments to pass to \code{estimates}.
#' @return numeric vector.
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' # see the example of estimation
#'
#' @export intliz
intliz <-
  function(data, ini.run=100, tol.ini=0.001, proc="ML", intv.s1=c(0,5), intv.s2=c(0,5),
           intv.a0=c(0,5), intv.a1=c(0,5), intv.a2=c(0,5), ...){
  if(!(proc  %in% c("ML", "S.EM"))) stop("'proc' must be 'ML' or 'S.EM'")
  ini.mat <- matrix(0, ini.run, 6)
  ini.est <- matrix(0, ini.run, 5)
  Nt <- 0
  for(i in 1:ini.run) {
     ini.mat[i, 1] <- runif(1, intv.s1[1], intv.s1[2])
     ini.mat[i, 2] <- runif(1, intv.s2[1], intv.s2[2])

     ini.mat[i, 3] <- runif(1, intv.a0[1], intv.a0[2])
     ini.mat[i, 4] <- runif(1, intv.a1[1], intv.a1[2])
     ini.mat[i, 5] <- runif(1, intv.a2[1], intv.a2[2])

     res <- estimates(data, ini.mat[i, 1], ini.mat[i, 2], ini.mat[i, 3], ini.mat[i, 4],
                      ini.mat[i, 5], tol.est = tol.ini, ...)
     mu1 <- res$mu1
     mu2 <- res$mu2
     ini.est[i,1] <- res$sigma1
     ini.est[i,2] <- res$sigma2
     ini.est[i,3] <- res$alpha0
     ini.est[i,4] <- res$alpha1
     ini.est[i,5] <- res$alpha2
     res1 <- logL(data, mu1, mu2, ini.est[i,1], ini.est[i,2], ini.est[i,3], ini.est[i,4],
                  ini.est[i,5])
     ini.mat[i, 6] <- res1$logLik
     Nt <- Nt + 1
     cat('\n', ini.mat[i, 6], res1$n1, res1$n2, ini.est[i,1], ini.est[i,2], ini.est[i,3], ini.est[i,4], ini.est[i,5], Nt,'\n')
  }
  if(proc=="ML"){
    return(apply(ini.est, 2, min))
  }else if(proc=="S.EM"){
    return(ini.mat[which.max(ini.mat[,6]),1:5])
  }else stop("'proc' must be 'ML' or 'S.EM'")
}
