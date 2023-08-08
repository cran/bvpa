#' Initialization of 3-parameter Block-Basu Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Return initial choice parameters of 3-parameter BBBVPA distribution.
#'
#' @param data bivariate observations.
#' @param ini.run number of random initializations.
#' @param tol.ini convergence tolerance, \code{0.001} (default)..
#' @param proc different procedures, \code{"ML"} (default) and \code{"S.EM"}.
#' @param intv.a0 interval for random initialization of \eqn{\alpha_0}.
#' @param intv.a1 interval for random initialization of \eqn{\alpha_1}.
#' @param intv.a2 interval for random initialization of \eqn{\alpha_2}.
#' @param ... further arguments to pass to \code{estimates3}.
#' @return  numeric vector.
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' dat <- rbb.bvpa(500, 0, 0, 1.0, 1.0, 2.0, 0.4, 0.5)
#' intliz3(dat)
#'
#' @export intliz3
intliz3 <-
  function(data, ini.run=100, tol.ini=0.001, proc="ML", intv.a0=c(0,5), intv.a1=c(0,5),
           intv.a2=c(0,5), ...){
    if(!(proc  %in% c("ML", "S.EM"))) stop("'proc' must be 'ML' or 'S.EM'")
    ini.mat <- matrix(0, ini.run, 4)
    ini.est <- matrix(0, ini.run, 3)

    Nt <- 0
    for(i in 1:ini.run) {
      ini.mat[i, 1] <- runif(1, intv.a0[1], intv.a0[2])
      ini.mat[i, 2] <- runif(1, intv.a1[1], intv.a1[2])
      ini.mat[i, 3] <- runif(1, intv.a2[1], intv.a2[2])

      res <- estimates3(data, ini.mat[i, 1], ini.mat[i, 2], ini.mat[i, 3], tol.est = tol.ini, ...)
      ini.est[i,1] <- res$alpha0
      ini.est[i,2] <- res$alpha1
      ini.est[i,3] <- res$alpha2
      res1 <- logL(data, 0, 0, 1, 1, ini.est[i,1], ini.est[i,2], ini.est[i,3])
      ini.mat[i, 4] <- res1$logLik

      Nt <- Nt + 1
      cat('\n', ini.mat[i, 4], res1$n1, res1$n2, ini.est[i,1], ini.est[i,2], ini.est[i,3], Nt,'\n')
    }
    if(proc=="ML"){
      return(apply(ini.est, 2, min))
    }else if(proc=="S.EM"){
      return(ini.mat[which.max(ini.mat[,4]), 1:3])
    }else stop("'proc' must be 'ML' or 'S.EM'")
  }
