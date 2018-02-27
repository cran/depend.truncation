\name{Logrank.stat}
\alias{Logrank.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The weighted log-rank statistics for testing quasi-independence (without ties in data)
}
\description{
The three log-rank statistics (L_0, L_1, and L_log) corresponding to 3 different weights.
}
\usage{
Logrank.stat(x.trunc, z.trunc, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.trunc}{vector of variables satisfying x.trunc<=z.trunc}
  \item{z.trunc}{vector of variables satisfying x.trunc<=z.trunc}
  \item{d}{censoring indicator(0=censoring,1=failure) for z.trunc}
}
\details{
If there is no tie in the data, the function "Logrank.stat.tie" and "Logrank.stat" give identical results.
However, "Logrank.stat" is computationally more efficient. The simulations of Emura & Wang (2010) are 
based on "Logrank.stat" since simulated data are generated from continuous distributions. The real data analyses 
of Emura & Wang (2010) are based on "Logrank.stat.tie" since there are many ties in the data.
}
\value{
  \item{L0}{Logrank statistics (most powerfull to detect the Clayton copula type dependence)}
  \item{L1}{Logrank statistics (most powerfull to detect the Frank copula type dependence)}
  \item{Llog}{Logrank statistics (most powerfull to detect the Gumbel copula type dependence)}
}
\references{
Emura T, Wang W (2010) Testing quasi-independence for truncation data. Journal of Multivariate Analysis 101, 223-239
}
\author{Takeshi Emura}

\examples{
x.trunc=c(10,5,7,1,3,9)
z.trunc=c(12,11,8,6,4,13)
d=c(1,1,1,1,0,1)
Logrank.stat(x.trunc,z.trunc,d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Copula }
\keyword{ Quasi-independence test }% __ONLY ONE__ keyword per line
