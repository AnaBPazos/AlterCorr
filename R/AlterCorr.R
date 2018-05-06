#' Alternatives to the Pearson correlation with p-values
#'
#' This function returns the value of the correlation between two variables and
#' its corresponding p-value for four types of coefficients; Pearson correlation
#' (r), Maximum information coefficient (MIC); Random dependency coefficient (RDC)
#' and Correlation of distances (dCor).
#'
#' @param x A numeric vector
#' @param y A numeric vector of the same dimension as \code{x}
#' @param type Specifies the type of correlations to compute. "Pearson" for the
#' Pearson Correlation, "MIC" for the Maximum Information Coefficient, "RDC" for
#' the Random Dependency Coefficient and finally "dCor" for the Correlation of
#' Distances.
#' @param R Number of permutations to be made for MIC and dCor to calculate
#' the p-values
#'
#' @keywords correlation
#' @import Hmisc
#' @import minerva
#' @import energy
#' @export
#' @examples
#' #Lineal dependence
#' x<-rnorm(200, sd=2.5)
#' y<-2*x+rnorm(200, sd=2.5)
#' AlterCorr(x, y, type="pearson")
#' AlterCorr(x, y, type="MIC",R=10)
#' AlterCorr(x, y, type="RDC")
#' AlterCorr(x, y, type="dCor",R=100)
#'
#' #No dependency
#' x<-rnorm(200, mean=2, sd =1.8)
#' y<-rnorm(200, mean=0, sd =2.5)
#' AlterCorr(x, y, type="pearson")
#' AlterCorr(x, y, type="MIC",R=10)
#' AlterCorr(x, y, type="RDC")
#' AlterCorr(x, y, type="dCor",R=100)

AlterCorr<-function(x,y,type=c("pearson","MIC","RDC","dCor"),R=100){
  type <- match.arg(type)
  if (! is.null(R)) {
    R <- floor(R)
    if (R < 1) R <- 100
  } else {
    R <- 100
  }
  if (type=="pearson") {
    corrP <- rcorr(x,y, type="pearson")
    result<-as.data.frame(cbind(corrP$r[1,2], corrP$P[1,2]))
  }
  if (type=="MIC") {
    result <- MIC(x, y, R=R)
  }
  if (type=="RDC") {
    result <- rdc(x, y, k=20, s=1/6, f=sin)
  }
  if (type=="dCor") {
    dcor <- dcor.test(x, y, index=1.0, R=R)
    result<-as.data.frame(cbind(dcor$statistic, dcor$p.value))
  }
  colnames(result)<-c("Correlation","pvalue")
  return(result)
}

