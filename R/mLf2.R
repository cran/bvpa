#' Marginal log-liklihood function of variable X2
#'
#' @description
#' Return the marginal log-liklihood value of variable \eqn{X_2}.
#'
#' @param I baivariate observations.
#' @param mu2 value of \eqn{\mu_2}.
#' @param s2 value of \eqn{\sigma_2}.
#' @param a0 value of \eqn{\alpha_0}.
#' @param a1 value of \eqn{\alpha_1}.
#' @param a2 value of \eqn{\alpha_2}.
#'
#' @return A scalar numeric, the marginal log-liklihood value of variable \eqn{X_2}.
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' dat <- rbb.bvpa(500, 0.1, 0.1, 0.8, 0.8, 2.0, 0.4, 0.5)
#' mLf2(dat, 0.1, 0.8, 2.0, 0.4, 0.5)
#'
#' @export mLf2
mLf2 <-
  function(I, mu2, s2, a0, a1, a2){
    fn2 = ((a0 + a1 + a2)/(a1 + a2))*((a0 + a2)/s2)*(1 + (I[,2] - mu2)/s2)^(- a0 - a2 - 1) - (a0/(a1 + a2))*((a0 + a1 + a2)/s2)*(1 + (I[,2] - mu2)/s2)^(-a0 - a1 - a2 - 1)
    ff2 <- sum(log(fn2))
    return(ff2)
  }
