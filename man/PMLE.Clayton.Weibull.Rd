\name{PMLE.Clayton.Weibull}
\alias{PMLE.Clayton.Weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric Inference for Models with Dependent Truncation Data}
\description{Maximum likelihood estimation (MLE) for dependent truncation data under 
the Clayton copula with Weibull margins for a bivariate lifetimes (L, X). The truncated data (L_j, X_j), subject to L_j<=X_j for all j=1, ..., n,
 are used to obtain the MLE for the population parameters of (L, X).}
\usage{
PMLE.Clayton.Weibull(l.trunc, x.trunc, GOF = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l.trunc}{vector of truncation variables satisfying l.trunc<=x.trunc}
  \item{x.trunc}{vector of variables satisfying l.trunc<=x.trunc}
  \item{GOF}{if TRUE, a goodness-of-fit test statistics is computed}
}

\details{Relevant paper is submitted for review}
\value{
 \item{n }{sample size}
 \item{alpha }{alpha and its standard error}
 \item{lambda_L }{hazard of L and its standard error}
 \item{lambda_X }{hazard of X and its standard error}
 \item{nu_L }{shape parameter of L and its standard error}
 \item{nu_X }{shape parameter of X and its standard error}
 \item{c }{ inclusion probability, defined by c=Pr(L<=X)}
 \item{C }{ Cramer-von Mises goodness-of-fit test statistics}
 \item{K }{ Kolmogorov-Smirnov goodness-of-fit test statistics}
}
\references{
Emura T, Pan CH (2016), Parametric maximum likelihood inference for copula models for dependently left-truncated data, submitted for publication.
}
\author{Takeshi Emura}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
########### Brake pad data from Lawless (2003) ##############
l.trunc = c( 22.2 , 23.0 , 24.0 , 28.6 , 21.8 , 17.0 , 26.0 , 23.2 , 18.9 , 21.9 , 27.3 ,
       13.8 , 24.0 , 20.1 , 15.7 , 26.8 , 27.9 , 15.3 , 28.8 , 16.0 , 23.6 , 53.8 ,
       21.7 , 28.8 , 17.0 , 16.5 , 15.7 , 28.0 , 13.3 , 16.5 , 24.2 , 17.6 , 27.8 , 
       18.3 , 17.7 , 20.0 , 13.2 , 16.9 , 14.9 , 15.5 ,  7.0 , 15.8 , 15.0 , 38.3 , 
       11.2 , 38.2 , 26.7 , 17.1 , 29.0 , 18.3 , 18.4 , 18.2 , 15.9 , 16.4 , 23.6 , 
       19.2 , 23.3 , 20.4 , 20.9 , 28.5 , 23.2 , 17.9 , 46.1 , 39.3 , 11.8 , 17.7 , 
       30.9 , 22.4 , 45.0 , 18.2 , 30.2 , 21.8 , 18.2 , 23.0 , 27.2 , 10.9 , 25.5 , 
       12.4 , 39.9 , 17.7 , 26.3 , 14.1 , 21.0 , 11.2 , 10.8 , 25.7 , 32.4 , 13.6 , 
       19.1 , 16.1 , 53.3 , 57.3 , 36.5 , 19.7 , 20.8 , 30.8 , 20.0 , 39.6 )

x.trunc = c( 38.7 , 49.2 , 42.4 , 73.8 , 46.7 , 44.1 , 61.9 , 39.3 , 49.8 , 46.3 , 56.2 , 
       50.5 , 54.9 , 54.0 , 49.2 , 44.8 , 72.2 , 107.8 , 81.6 , 45.2 , 124.6 , 64.0 , 
       83.0 , 143.6 , 43.4 , 69.6 , 74.8 , 32.9 , 51.5 , 31.8 , 77.6 , 63.7 , 83.0 , 
       24.8 , 68.8 , 68.8 , 89.1 , 65.0 , 65.1 , 59.3 , 53.9 , 79.4 , 47.4 , 61.4 , 
       72.8 , 54.0 , 37.2 , 44.2 , 50.8 , 65.5 , 86.7 , 43.8 , 100.6 , 67.6 , 89.5 , 
       60.3 , 103.6 , 82.6 , 88.0 , 42.4 , 68.9 , 95.7 , 78.1 , 83.6 , 18.6 , 92.6 , 
       42.4 , 34.3 , 105.6 , 20.8 , 52.0 , 77.2 , 68.9 , 78.7 , 165.5 , 79.5 , 55.0 , 
       46.8 , 124.5 , 92.5 , 110.0 , 101.2 , 59.4 , 27.8 , 33.6 , 69.0 , 75.2 , 58.4 , 
       105.6 , 56.2 , 55.9 , 83.8 , 123.5 , 69.0 , 101.9 , 87.6 , 38.8 ,74.7 )

PMLE.Clayton.Weibull(l.trunc,x.trunc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Parametric method }% __ONLY ONE__ keyword per line
