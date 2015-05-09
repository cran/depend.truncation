\name{depend.truncation-package}
\alias{depend.truncation-package}
\alias{depend.truncation}
\docType{package}
\title{Statistical Inference for Parametric and Semiparametric Models Based on Dependently Truncated Data}
\description{
Suppose that one can observe bivariate random variables (X, Y) only when X<=Y holds.
Data (Xj, Yj), subject to Xj<=Yj, for all j=1,...,n, are called truncated data
 (or truncation data).
 This package performs parametric and semiparametric 
inference methods for models on (X, Y), especially when X and Y are dependent (correlated).
Semi-parametric approaches are based on Chaieb et al. (2006), Emura et al. (2011), Emura and Wang (2012) and Emura and Murotani (2015). Parametric approaches are based on Emura and Konno (2012).
 A regression approach is based on Emura and Wang (2015).
 This package also includes truncated data on the number of deaths at each year (1963-1980) for Japanese male centenarians (Emura and Murotani 2015). 
}

\details{
\tabular{ll}{
Package: \tab depend.truncation\cr
Type: \tab Package\cr
Version: \tab 2.3\cr
Date: \tab 2015-05-09\cr
License: \tab GPL-2\cr
}
}
\author{
Takeshi Emura <emura@stat.ncu.edu.tw>
}
\references{
Chaieb L, Rivest LP, Abdous B (2006), Estimating Survival Under a Dependent Truncation, Biometrika 93: 655-669.

Emura T, Wang W, Hung HN (2011), Semiparametric Inference for Copula Models for Truncated Data, Stat Sinica 21: 349-367.

Emura T, Konno Y (2012), Multivariate Normal Distribution Approaches for Dependently Truncated Data.
Statistical Papers 53 (No.1), 133-149.

Emura T, Wang, W (2012), Nonparametric Maximum Likelihood Estimation for Dependent Truncation Data Based on Copulas, Journal of Multivariate Analysis 110, 171-88

Emura T, Murotani K (2015), An Algorithm for Estimating Survival Under a Copula-based Dependent Truncation Model, TEST. DOI: 10.1007/s11749-015-0432-8

Emura T, Wang W (2015), Semiparametric Inference for an Accelerated Failure Time Model with Dependent Truncation, to appear, Annals of the Institute of Statistical Mathematics.
}