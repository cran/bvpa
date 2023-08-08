#' Simulate from a Block-Basu Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Produces one or more samples from the specified BBBVPA distribution.
#'
#' @param n number of observations.
#' @param mu1 value of \eqn{\mu_1}
#' @param mu2 value of \eqn{\mu_2}
#' @param sig1 value of \eqn{\sigma_1}
#' @param sig2 value of \eqn{\sigma_2}
#' @param alp0 value of \eqn{\alpha_0}
#' @param alp1 value of \eqn{\alpha_1}
#' @param alp2 value of \eqn{\alpha_2}
#' @return  numeric matrix.
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' cor(rbb.bvpa(500, 0.1, 0.1, 0.8, 0.8, 2.0, 0.4, 0.5))
#'
#' @importFrom stats runif
#'
#' @export rbb.bvpa
rbb.bvpa <-
  function(n, mu1, mu2, sig1, sig2, alp0, alp1, alp2){
  g <- matrix(0,n,2)
  t <- 0
  i <- 1
  while(t < 1){
    u0 <- (1 - runif(1,0,1))^(-1/alp0) - 1
    u1 <- mu1 + ((1 - runif(1,0,1))^(-1/alp1) - 1)*sig1
    u2 <- mu2 + ((1 - runif(1,0,1))^(-1/alp2) - 1)*sig2

    if(((mu1 + sig1*u0) > u1) || ((mu2 + sig2*u0) > u2)){
      g[i,1] <- min((mu1 + sig1*u0),u1)
      g[i,2] <- min((mu2 + sig2*u0),u2)
      i <- i + 1
      if(i == (n+1)){
        break
      }
    }
  }
  return(g)
}


