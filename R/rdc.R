#' Randomized Dependence Coefficient (RDC)
#'
#' This function calculates the Random Dependence Coefficient (RDC).
#' Calculate the dependence between the random samples as the highest canonical
#' correlation between the k random non-linear projections of their copula
#' transformations.
#'
#' It also provides the p-value for the independence hypothesis, assuming the
#' normality of the data through the Bartlett's approximation.
#'
#' @param x A vector, matrix or numeric data frame
#' @param y A vector, matrix or numeric data frame
#' @param k Number of non-linear projections of the copula, by default k = 20
#' @param s Variance to draw i.i.d. projection coefficients in N ~(0, sI),
#' by defect is 1/6
#' @param f Function that is used for the generation of random non-linear
#' projections, if it is not indicated it uses the sinusoidal projections (sin)
#'
#' @importFrom stats cancor
#' @export
#' @examples
#' #Lineal dependence
#' x<-rnorm(200, sd=2.5)
#' y<-2*x+rnorm(200, sd=2.5)
#' rdc(x,y)
#'
#' #No dependency
#' x<-rnorm(200, mean=2, sd =1.8)
#' y<-rnorm(200, mean=0, sd =2.5)
#' rdc(x,y)


rdc <- function(x,y,k=20,s=1/6,f=sin) {
  x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
  y <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)),1)
  x<- s/ncol(x)*x%*%matrix(rnorm(ncol(x)*k),ncol(x))
  y <- s/ncol(y)*y%*%matrix(rnorm(ncol(y)*k),ncol(y))
  can<-cancor(cbind(f(x),1),cbind(f(y),1))$cor

  k<-length(can)
  chi<-(((2*k+3)/2)-nrow(x))*log(prod(1-can^2))
  rdc<-as.data.frame(cbind(can[1],pchisq(chi,k^2, lower.tail = FALSE)))
  colnames(rdc)<-c("rdc","p-value")
  return(rdc)
}

