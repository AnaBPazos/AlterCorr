# AlterCorr"

#' Alternatives to the Pearson correlation with p-values for Matrix
#'
#' This package calculate the correlations between two variables or matrices 
#' or dataframes, for four independence coefficients; Pearson correlation, 
#' Maximum Information Coefficient (MIC), Random dependency coefficient (RDC) 
#' and Correlation of distances (dCor). These coefficients differ from the 
#' Pearson correlation in that they detect non-linear dependencies and are not 
#' subject to the normality of the data. Together with the correlations, 
#' calculate the p-values and the adjusted p-values for each type of coefficient. #' When the function is used for matrices and dataframes, the difference between
#' the calculation of paired data (must have the same dimensions) or independent 
#' (they must only coincide in the same number of columns).
#' 
#' 
#' It consists of four functions:
#' <B>AlterCorr</B> that does the calculations for vectors of variables.
#' <B>AlterCorrM</B> that performs calculations for value matrices using the
#' AlterCorr function.
#' The <B>MIC</B> function calculates the maximum information coefficient and its
#' p-values with a permutations test.
#' The <B>RDC</B> function calculates the coefficient of random dependence, 
#' together with the significance test using the Bartlett approximation.
