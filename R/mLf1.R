#' Marginal log-liklihood function of variable X1
#'
#' @description
#' Return the marginal log-liklihood value of variable \eqn{X_1}.
#'
#' @param I baivariate observations.
#' @param mu1 value of \eqn{\mu_1}.
#' @param s1 value of \eqn{\sigma_1}.
#' @param a0 value of \eqn{\alpha_0}.
#' @param a1 value of \eqn{\alpha_1}.
#' @param a2 value of \eqn{\alpha_2}.
#'
#' @return A scalar numeric, the marginal log-liklihood value of variable \eqn{X_1}.
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' dat <- rbb.bvpa(500, 0.1, 0.1, 0.8, 0.8, 2.0, 0.4, 0.5)
#' mLf1(dat, 0.1, 0.8, 2.0, 0.4, 0.5)
#'
#' @export mLf1
mLf1 <-
  function(I, mu1, s1, a0, a1, a2){
    fn1 = ((a0 + a1 + a2)/(a1 + a2))*((a0 + a1)/s1)*(1 + (I[,1] - mu1)/s1)^(- a0 - a1 - 1) - (a0/(a1 + a2))*((a0 + a1 + a2)/s1)*(1 + (I[,1] - mu1)/s1)^(-a0 - a1 - a2 - 1)
    ff1 <- sum(log(fn1))
    return(ff1)
  }
