#' Estimation of Block-Basu Bivariate Pareto (BBBVPA) distribution
#'
#' @description
#' Parameters estimation of BBBVPA distribution.
#'
#' @param I bivariate observations.
#' @param s1.int initial choice of \eqn{\sigma_1}.
#' @param s2.int initial choice of \eqn{\sigma_2}.
#' @param a0.int initial choice of \eqn{\alpha_0}.
#' @param a1.int initial choice of \eqn{\alpha_1}.
#' @param a2.int initial choice of \eqn{\alpha_2}.
#' @param tol.est convergence tolerance, \code{0.00001} (default).
#' @param MxIter.no maximum number of iterations, \code{2000} (default).
#' @param rate step size or learning rate for gradient descent, \code{0.0001} (default).
#' @param condition convergence criterion, \code{"log.L"} (default) and \code{"p.logL"}.
#'
#' @return  object of class "\code{bbbvpa}", a list consisting of
#' \item{mu1, mu2, sigma1, sigma2, alpha0, alpha1, alpha2, iter.no}{estimates  of parameters and number of iteration.}
#' \item{data }{the supplied data \code{I}.}
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Arabin Kumar Dey <arabin@iitg.ac.in>
#'
#' @examples
#' \donttest{
#' # Read data
#' data(precipitation)
#' data <- as.vector(precipitation[,2])
#' data[is.na(data)]<-0
#' n <- length(data)
#' # Construct the three-dimensional data set
#' data3d <- function(data){
#'  u <- 12
#'  Y <- c()
#'  indx <- indx1 <- indx2 <- indx3 <- 0
#'  r <- 5
#'  i <- 2
#'  while(i < n){
#'    i <- i + 1
#'    if(data[i] > u || sum(data[(i-1):i]) > u || sum(data[(i-2):i]) > u){
#'      if(data[i] > u){imax <- i}
#'      if(sum(data[(i-1):i]) > u){imax <- i - 3 + which(data[(i-1):i] == max(data[(i-1):i]))[1]}
#'      if(sum(data[(i-2):i]) > u){imax <- i - 3 + which(data[(i-2):i] == max(data[(i-2):i]))[1]}
#'      if(max(indx) > (imax-r)){
#'        cluster <- data[(max(indx)+3):(imax+r)]
#'      } else{
#'        cluster <- data[(imax-r):(imax+r)]
#'      }
#'      cluster2 <- sapply(c(1:(length(cluster)-1)), function(j) sum(cluster[j:(j+1)]))
#'      cluster3 <- sapply(c(1:(length(cluster)-2)), function(j) sum(cluster[j:(j+2)]))
#'      indx1 <- append(indx1,imax-r-1+which(cluster==max(cluster))[1])
#'      indx2 <- append(indx2,imax-r-1+which(cluster2==max(cluster2)))
#'      indx3 <- append(indx3,imax-r-1+which(cluster3==max(cluster3)))
#'      Y <- rbind(Y, c(max(cluster),max(cluster2),max(cluster3)))
#'      indx <- append(indx,imax)
#'      i <- i + r
#'    }
#'  }
#'  return(Y)
#' }
#' I <- data3d(data)[,c(1,3)]
#' iniz <- intliz(I)
#' iniz
#' est <- estimates(I, iniz[1], iniz[2], iniz[3], iniz[4], iniz[5])
#' est[-9]
#' param.boot(I, iniz[1], iniz[2], iniz[3], iniz[4], iniz[5])
#' conf.intv(est)
#' }
#'
#' @importFrom numDeriv grad
#'
#' @export estimates
estimates <-
  function(I, s1.int, s2.int, a0.int, a1.int, a2.int, tol.est=0.00001, MxIter.no=2000,
                      rate=0.0001, condition="log.L"){
  n <- nrow(I)
  mufin1 <- min(I[,1])
  mufin2 <- min(I[,2])

  if(condition=="p.logL"){
    Qold <- pseu.logL(I, mufin1, mufin2, s1.int, s2.int, a0.int, a1.int, a2.int)
  } else if(condition=="log.L"){
    Qold <- logL(I, mufin1, mufin2, s1.int, s2.int, a0.int, a1.int, a2.int)$logLik
  }else stop("'condition' must be 'p.logL' or 'log.L'")
  t = 0
  j = 0
  while(j < 1)
  {
    sigfin1 <- s1.int + rate*grad(mLf1, s1.int, mu1=mufin1, a0=a0.int, a1=a1.int, a2=a2.int, I = I)
    sigfin2 <- s2.int + rate*grad(mLf2, s2.int, mu2=mufin2, a0=a0.int, a1=a1.int, a2=a2.int, I = I)

    Iv1 <- (I[, 1] - mufin1)/sigfin1
    Iv2 <- (I[, 2] - mufin2)/sigfin2
    I1 <- cbind(Iv1, Iv2)

    I11 <- cbind(I1[I1[, 1] < I1[,2],1], I1[I1[, 1] < I1[,2],2])
    I12 <- cbind(I1[I1[, 1] > I1[,2],1], I1[I1[, 1] > I1[,2],2])

    # n1 <- length(I11[,1])
    # n2 <- length(I12[,1])

    e.n1 <- length(I11[,1])
    e.n2 <- length(I12[,1])
    # n <- (e.n1+e.n2)

    n1 <- (e.n1+e.n2)*a1.int/(a1.int+a2.int)
    n2 <- (e.n1+e.n2)*a2.int/(a1.int+a2.int)

    u1 <- a0.int/(a0.int + a2.int)
    u2 <- a2.int/(a0.int + a2.int)
    w1 <- a0.int/(a0.int + a1.int)
    w2 <- a1.int/(a0.int + a1.int)

    n0 <- (n1 + n2)*(a0.int/(a1.int+a2.int))
    a <- 1/(a0.int + a1.int + a2.int)

    alpfin0 <- (n0 + u1*n1 + w1*n2)/(n0*a + sum(log(1 + I12[,1])) + sum(log(1 + I11[,2])))
    alpfin1 <- (n1 + w2*n2)/(n0*a + sum(log(1 + I11[,1])) + sum(log(1 + I12[,1])))
    alpfin2 <- (n2 + u2*n1)/(n0*a + sum(log(1 + I11[,2])) + sum(log(1 + I12[,2])))

    if(condition=="p.logL"){
      QQ <- pseu.logL(I, mufin1, mufin2, sigfin1, sigfin2, alpfin0, alpfin1, alpfin2)
    } else if(condition=="log.L"){
      QQ <- logL(I, mufin1, mufin2, sigfin1, sigfin2, alpfin0, alpfin1, alpfin2)$logLik
    }else stop("'condition' must be 'p.logL' or 'log.L'")

    delQ <- QQ - Qold
    t <- t+1
    if((abs(delQ/Qold) < tol.est)||(t >= MxIter.no)){
      break
    }
    s1.int <- sigfin1
    s2.int <- sigfin2
    a0.int <- alpfin0
    a1.int <- alpfin1
    a2.int <- alpfin2
    Qold <- QQ
  }
  res <- list(mu1=mufin1, mu2=mufin2, sigma1=sigfin1, sigma2=sigfin2,
              alpha0=alpfin0, alpha1=alpfin1, alpha2=alpfin2, iter.no=t)
  res$data <- I
  class(res) <- c("bbbvpa")
  return(res)
}
