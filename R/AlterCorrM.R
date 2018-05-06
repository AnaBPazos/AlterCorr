#' Alternatives to the Pearson correlation with p-values for Matrix
#'
#' This function returns the value of the correlation between two matrix and
#' its corresponding p-value for four types of coefficients; Pearson correlation
#' (r), Maximum information coefficient (MIC); Random dependency coefficient (RDC)
#' and Correlation of distances (dCor).
#'
#' Depending on the selected comparison option, the correlation between each
#' value of the matrix will be calculated with respect to the other matrix or
#' it will be calculated in pairs.
#'
#'
#' @param x A numeric vector or matrix or dataframe
#' @param y A numeric matrix or dataframe, with the same number of rows as \code{x}
#' @param type Specifies the type of correlations to compute. "Pearson" for the
#' Pearson Correlation, "MIC" for the Maximum Information Coefficient, "RDC" for
#' the Random Dependency Coefficient and finally "dCor" for the Correlation of
#' Distances.
#' @param comparison Method with which the correlations will be calculated.
#' The "pairs" option is to compare two matrices that share rows and columns but
#' the values are of different variables and we want to make pairwise comparisons.
#' The case of "alls" is for when you have two matrices with a common dimension
#' (observations, individuals, samples), and you want to calculate the correlations
#' of all the variables with all.
#' @param R Number of permutations to be made for MIC and dCor to calculate
#' the p-values
#'
#' @note It is not recommended to use the MIC coefficient for large matrix or
#' with a large number of repetitions, since it requires a high computational
#' time.
#'
#' @keywords correlation
#' @import Hmisc
#' @import minerva
#' @import energy
#' @export
#' @examples
#' #Correlations using the "pair" comparison
#' x<-replicate(20,rnorm(20, sd=2.5))
#' y<-2*x+replicate(20,rnorm(20, sd=2.5))
#' AlterCorrM(x,y,type="pearson",comparison="pairs")
#' AlterCorrM(x,y,type="MIC",comparison="pairs",R=50)
#' AlterCorrM(x,y,type="RDC",comparison="pairs")
#' AlterCorrM(x,y,type="dCor",comparison="pairs",R=50)
#'
#' #Correlations using the "all" comparison
#' x<-replicate(5,rnorm(10, sd=2.5))
#' colnames(x)<-letters[1:ncol(x)]
#' y<-replicate(7,rnorm(10, mean=1.5, sd =0.85))
#' colnames(y)<-letters[1:ncol(y)]
#'
#' AlterCorrM(x,y,type="pearson",comparison="all")
#' AlterCorrM(x,y,type="MIC",comparison="all",R=50)
#' AlterCorrM(x,y,type="RDC",comparison="all")
#' AlterCorrM(x,y,type="dCor",comparison="all",R=50)
#'

AlterCorrM<-function(X,Y,type=c("pearson","MIC","RDC","dCor"),
                     comparison=c("all","pairs"),R=100){
  type <- match.arg(type)
  comp <- match.arg(comparison)

  if (! is.null(R)) {
    R <- floor(R)
    if (R < 1) R <- 100
  } else {
    R <- 100
  }

  if ((nrow(X)!=nrow(Y)) && type=="all") stop('The number of observations does not match')
  if ((nrow(X)!=nrow(Y) ||(ncol(X)!=ncol(Y))) && type=="pairs") {
    stop('Matrices dimensions do not match, select comparation="alls"')
  }

  if (comp=="all"){
    Corrs<-matrix(NA,nrow=ncol(X),ncol=ncol(Y))
    Pvalue<-matrix(NA,nrow=ncol(X),ncol=ncol(Y))
    for (i in 1:ncol(X)){
      for (j in 1:ncol(Y)){
        z<-AlterCorr(X[,i], Y[,j], type=type,R=R)
        Corrs[i,j] <- z$Correlation
        Pvalue[i,j] <- z$pvalue
      }
    }
    colnames(Corrs)<-colnames(Y)
    colnames(Pvalue)<-colnames(Y)
    rownames(Corrs)<-colnames(X)
    rownames(Pvalue)<-colnames(X)
  }

  if (comp=="pairs"){
    Corrs<-matrix(NA,nrow=ncol(X),ncol=1)
    Pvalue<-matrix(NA,nrow=ncol(X),ncol=1)
    for (i in 1:ncol(X)){
      z<-AlterCorr(X[,i], Y[,i], type=type,R=R)
      Corrs[i,1] <- z$Correlation
      Pvalue[i,1] <- z$pvalue
    }
    rownames(Corrs)<-colnames(X)
    rownames(Pvalue)<-colnames(X)
  }

  res<-(list(Correlation=Corrs, pvalue=Pvalue))
  return(res)
}

