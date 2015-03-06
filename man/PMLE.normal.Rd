\name{PMLE.Normal}
\alias{PMLE.Normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric Inference for Models with Dependent Truncation Data}
\description{Maximum likelihood estimation (MLE) for dependent truncation data under 
the bivariate normal distribution. A bivariate normal distribution is assumed for bivariate
 random variables (L, X). The truncated data (L_j, X_j), subject to L_j<=X_j for all j=1, ..., n,
 are used to obtain the MLE for the population parameters of (L, X).}
\usage{
PMLE.Normal(l.trunc, x.trunc, testimator = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l.trunc}{vector of truncation variables satisfying l.trunc<=x.trunc}
  \item{x.trunc}{vector of variables satisfying l.trunc<=x.trunc}
  \item{testimator}{if TRUE, testimator is computed instead of MLE}
}
\details{
PMLE.Normal performs the maximum likelihood estimation for dependently left-truncated data under 
the bivariate normal distribution. "PMLE.Normal" implements the methodologies developed in 
Emura T. & Konno Y. (2012, Statistical Papers 53, 133-149)and can produce the maximum likelihood estimates
 and their standard errors. Furthermore, "PMLE.Normal" tests the independence assumption between truncation variable
and variable of interest via likelihood ratio test. The MLE is obtained by minimizing -logL using "nlm", where L is the log-likelihood.
}
\value{
 \item{mu_L }{mean of L and its standard error}
 \item{mu_X }{mean of X and its standard error}
 \item{var_L }{variance of L and its standard error}
 \item{var_X }{variance of X and its standard error}
 \item{cov_LX }{covariance between L and X and its standard error}
 \item{prob }{ inclusion probability, defined by Pr(L<=X), and its standard error}
 \item{test }{ Likelihood ratio statistic and p-value}
}
\references{
Emura T, Konno Y (2012), Multivariate Normal Distribution Approaches for Dependently Truncated Data.
Statistical Papers 53 (No.1), 133-149.

Emura T, Konno Y (2014), Erratum to: Multivariate Normal Distribution Approaches for Dependently Truncated Data, Statistical Papers 55 (No.4): 1233-36
}
\author{Takeshi EMURA}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
l.trunc=c(1,2,3,4,5,6,7,8,8)
x.trunc=c(2,4,4,5,5,7,7,9,10)
PMLE.Normal(l.trunc,x.trunc,testimator=FALSE)
}

\keyword{ Parametric method }% __ONLY ONE__ keyword per line
