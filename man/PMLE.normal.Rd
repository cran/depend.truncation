\name{PMLE.normal}
\alias{PMLE.normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric Inference for models with dependent truncation data}
\description{Maximum likelihood estimation for dependent truncation data under 
the bivariate normal distribution}
\usage{
PMLE.normal(l.trunc, x.trunc, testimator = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l.trunc}{vector of truncation variables}
  \item{x.trunc}{vector of variables that precede "l.trunc"}
  \item{testimator}{if TRUE, testimator is computed instead of MLE}
}
\details{
PMLE.normal performs the maximum likelihood estimation for dependently left-truncated data under 
the bivariate normal distribution. "PMLE.normal" implements the methodologies developed in 
Emura T. & Konno Y. (2012, Statistical Papers 53, 133-149)and can produce the maximum likelihood estimates
 and their standard errors. Furthermore, "PMLE.normal" tests the independence assumption between truncation variable
and variable of interest via likelihood ratio test. The MLE is obtained by minimizing -logL using "nlm", where L is the log-likelihood.
}
\value{
 \item{mu_L }{mean of L and its standard error}
 \item{mu_X }{mean of X and its standard error}
 \item{var_L }{variance of L and its standard error}
 \item{var_X }{variance of X and its standard error}
 \item{cov_LX }{covariance between L and X and its standard error}
 \item{prob }{ Pr(L<=X) and its standard error}
 \item{test }{ Likelihood ratio statistic and p-value}
}
\references{
Emura, T. and Konno, Y. 2012, Multivariate Normal Distribution Approaches for Dependently Truncated Data.
Statistical Papers 53 (No.1), 133-149.
}
\author{Takeshi EMURA}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
l.trunc=c(1,2,3,4,5,6,7,8,8)
x.trunc=c(2,4,4,5,5,7,7,9,10)
PMLE.normal(l.trunc,x.trunc,testimator=FALSE)
}

\keyword{ Parametric method }% __ONLY ONE__ keyword per line
