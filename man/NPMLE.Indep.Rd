\name{NPMLE.Indep}
\alias{NPMLE.Indep}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Semiparametric Inference under the Independence Copula with Dependent Truncation}

\description{Nonparametric maximum likelihood estimation for dependent truncation data under 
the independence copula models (Emura & Wang, 2012). The forms of the marginal distributions for X and Y are completely unspecified, but
the parametric form of copula is specified as the independence copula. 
}

\usage{
NPMLE.Indep(x.trunc, y.trunc,
x.fix = median(x.trunc), y.fix = median(y.trunc), plotX = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.trunc}{vector of variables satisfying x.trunc<=y.trunc}
  \item{y.trunc}{vector of variables satisfying x.trunc<=y.trunc}
  \item{x.fix}{vector of fixed points at which marginal distribution function of X is calculated}
  \item{y.fix}{vector of fixed points at which marginal survival function of Y are calculated}
  \item{plotX}{if TRUE, the culumative distribution function for X is plotted}
}

\details{
The function produces the nonparametric maximum likelihood estimate (NPMLE) for the marginal distributions
 and the estimate of the association parameter under the independence copula model. The method follows Emura & Wang (2012).
The maximization of the likelihood (L) is conducted by minimizing -logL by "nlm". Standard errors are calculated
from the diagonal elements in the observed Fisher information matrix, which is obtained from the output of the "nlm".  
The method can handle ties but cannot handle right-censoring.
}
\value{
\item{Hx }{estimate of the marginal cumulative reverse-hazard function for X}
\item{Hx_se }{standard error of the estimate of the marginal cumulative reverse-hazard function for X}
\item{Ay }{estimate of the marginal cumulative hazard function for Y}
\item{Ay_se }{standard error of the estimate of the marginal cumulative hazard function for Y}
\item{Fx }{estimate of the marginal distribution function for X}
\item{Fx_se }{standard error of the estimate of the marginal distribution function for X}
\item{Sy }{estimate of the marginal survival function for Y}
\item{Sy_se }{standard error of the estimate of the marginal survival function for Y}
\item{conv }{an integer indicating why the optimization process terminated in "nlm": If conv=1,
the maximization of the likelihood function is properly done. Please refer R function "nlm"
 for more details.}
\item{Grad}{L_2 norm for the gradient vector at the solution (close to zero if the solution is proper)}
\item{MinEigen}{Minimum eigenvalue of the Hessian matrix at the solution (positive if the solution is proper)}
}



\references{
Emura T, Wang W (2012)  Nonparametric maximum likelihood estimation for dependent truncation data based on copulas, Journal of Multivariate Analysis 110, 171-88
}
\author{Takeshi EMURA}


\examples{
##### Generate independent truncated data with unit exponential marginals #####
set.seed(1)
m=35 #### sample size ####
x.trunc=y.trunc=numeric(m)
l=1
while(l<=m){
  x=rexp(1,rate=1)
  y=rexp(1,rate=1)
  if((x<y)){
  x.trunc[l]=x;y.trunc[l]=y
  l=l+1
  }
}

NPMLE.Indep(x.trunc, y.trunc, x.fix =1, y.fix = 1, plotX = TRUE)

#### compare the above estimate with the true value #####
c(Fx_true=1-exp(-1),Sy_true=exp(-1))
}

\keyword{ Copula }
\keyword{ Semiparametric method }% __ONLY ONE__ keyword per line
