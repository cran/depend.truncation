\name{PMLE.Clayton.Weibull}
\alias{PMLE.Clayton.Weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric Inference for Bivariate Weibull Models with Dependent Truncation}
\description{Maximum likelihood estimation (MLE) for dependent truncation data under 
the Clayton copula with Weibull margins for a bivariate lifetimes (L, X). The truncated data (L_j, X_j), subject to L_j<=X_j for all j=1, ..., n,
 are used to obtain the MLE for the population parameters of (L, X).}
\usage{
PMLE.Clayton.Weibull(l.trunc, x.trunc, GOF = TRUE,
 Err=2,alpha_max=20,alpha_min=10^-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l.trunc}{vector of truncation variables satisfying l.trunc<=x.trunc}
  \item{x.trunc}{vector of variables satisfying l.trunc<=x.trunc}
  \item{GOF}{if TRUE, a goodness-of-fit test statistics is computed}
  \item{Err}{tuning parameter in the NR algorithm}
  \item{alpha_max}{upper bound for the copula parameter}
  \item{alpha_min}{lower bound for the copula parameter}
}

\details{Relevant paper is submitted for review}
\value{
 \item{n }{sample size}
 \item{alpha }{dependence parameter}
 \item{lambda_L }{scale parameter of L}
 \item{lambda_X }{scale parameter of X}
 \item{nu_L }{shape parameter of L}
 \item{nu_X }{shape parameter of X}
 \item{mean_X }{Mean lifetime of X, defined as E[X]}
 \item{logL }{Maximized log-likelihood}
 \item{c }{ inclusion probability, defined by c=Pr(L<=X)}
 \item{C }{ Cramer-von Mises goodness-of-fit test statistics}
 \item{K }{ Kolmogorov-Smirnov goodness-of-fit test statistics}
}
\references{
Emura T, Pan CH (2017), Parametric likelihood inference and goodness-of-fit for dependently left-truncated data, a copula-based approach, Statistical Papers, doi:10.1007/s00362-017-0947-z.
}
\author{Takeshi Emura}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
l.trunc=c(22.207,23.002,23.982,28.551,21.789,17.042,25.997,23.220,18.854,21.857,
    27.321,13.767,23.982,20.110,15.779,26.821,27.934,15.292,28.843,15.985,
    23.580,53.770,21.731,28.844,17.046,16.506,15.696,27.959,13.272,16.482,
    24.210,17.626,27.770,
    18.264,17.694,20.014,13.152,16.886,14.894,15.531,6.951,15.841,14.974,
    38.292,11.204,38.156,26.652,17.101,28.953,18.325,18.391,18.220,15.896,
    16.447,23.642,19.170,23.257,20.428,20.947,28.462,23.210,17.900,46.134,
    39.300,11.768,17.717,
    30.863,22.350,44.976,18.169,30.164,21.822,18.201,22.895,27.189,10.915,
    25.503,12.350,39.869,17.698,26.296,14.091,21.011,11.201,10.757,25.692,
    32.372,13.592,19.102,16.112,53.281,57.298,36.450,19.651,20.755,30.788,20.0,39.62)

x.trunc = c(38.701,49.173,42.409,73.823,46.738,44.071,61.904,39.327,49.828,46.314,
    56.150,50.549,54.930,54.039,49.170,44.795,72.238,107.783,81.609,45.228,
    124.637,64.018,82.957,143.550,43.382,69.644,74.750,32.881,51.483,31.767,
    77.633,63.745,82.965,
    24.818,68.762,68.762,89.100,64.979,65.127,59.289,53.926,79.370,47.385,
    61.395,72.826,53.980,37.220,44.224,50.826,65.460,86.726,43.819,100.605,
    67.615,89.542,60.266,103.580,82.570,87.960,42.385,68.914,95.666,78.135,
    83.643,18.617,92.629, 
    42.415,34.346,106.569,20.758,52.003,77.179, 68.934,78.661,165.543,79.547,
    55.009,46.774,124.526,92.504,109.986,101.161,59.422,27.772,33.598,69.038,
    75.222,58.373,105.610,56.158,55.913,83.770,123.468,68.994,101.869,87.627,
    38.790,74.734)

u.min=10
l.trunc=l.trunc[-41]-u.min
x.trunc=x.trunc[-41]-u.min

PMLE.Clayton.Weibull(l.trunc,x.trunc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Parametric method }% __ONLY ONE__ keyword per line
