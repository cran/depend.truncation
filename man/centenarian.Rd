\name{centenarian}
\alias{centenarian}
\docType{data}
\title{
Japanese Centenarians Data
}
\description{
The number of deaths at each year (1963-1980) for Japanese male centenarians (Table 2 of Emura and Murotani (2015)). 
See also the original reference Sibuya & Hanayama (2004).
}
\usage{data("centenarian")}
\format{
  A data frame with 21 observations on the following 19 variables.
  \describe{
    \item{\code{X}}{:the age at death}
    \item{\code{X1963}}{:the number of deaths between 1963 and 1964}
    \item{\code{X1964}}{:the number of deaths between 1964 and 1965}
    \item{\code{X1965}}{:the number of deaths between 1965 and 1966}
    \item{\code{X1966}}{:the number of deaths between 1966 and 1967}
    \item{\code{X1967}}{:the number of deaths between 1967 and 1968}
    \item{\code{X1968}}{:the number of deaths between 1968 and 1969}
    \item{\code{X1969}}{:the number of deaths between 1969 and 1970}
    \item{\code{X1970}}{:the number of deaths between 1970 and 1971}
    \item{\code{X1971}}{:the number of deaths between 1971 and 1972}
    \item{\code{X1972}}{:the number of deaths between 1972 and 1973}
    \item{\code{X1973}}{:the number of deaths between 1973 and 1974}
    \item{\code{X1974}}{:the number of deaths between 1974 and 1975}
    \item{\code{X1975}}{:the number of deaths between 1975 and 1976}
    \item{\code{X1976}}{:the number of deaths between 1976 and 1977}
    \item{\code{X1977}}{:the number of deaths between 1977 and 1978}
    \item{\code{X1978}}{:the number of deaths between 1978 and 1979}
    \item{\code{X1979}}{:the number of deaths between 1979 and 1980}
    \item{\code{X1980}}{:the number of deaths between 1980 and 1981}
  }
}

\references{
Sibuya M, Hanayama N (2004), Estimation of Human Longevity Distribution Based on Tabulated Statistics. Proceeding of ISM 52: 117-34

Emura T, Murotani K (2015), An Algorithm for Estimating Survival Under a Copula-based Dependent Truncation Model, TEST 24 (No.4): 734-751.
}
\examples{
## Below is the centenarians data analysis of Emura & Murotani (2015) ##
data(centenarian)
Death=centenarian[,1]
Year=1963:1980
data.mat=centenarian[,-1]

X=T=NULL
for(i in 1:length(Death)){
  for(j in 1:length(Year)){
    X=c( X,rep(Death[i],data.mat[i,j]) )
    T=c( T,rep(Year[j]-i+1,data.mat[i,j]) )  ### T= Year at age 100.5 ###    
  }
}

x.trunc=X
z.trunc=max(Year)+0.5-T+100

m=length(x.trunc)
d=rep(1,m)

set.seed(1)
x.trunc=x.trunc+runif(length(x.trunc),min=-0.01,max=0)
z.trunc=z.trunc+runif(length(z.trunc),min=0,max=0.01)

### Copula-based estimator ####
## CHAIEB.Frank(x.trunc, z.trunc, d, a = 1/10)
}
\keyword{datasets}
