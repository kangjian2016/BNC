% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BNCmainfun.R
\name{biGaussianNull}
\alias{biGaussianNull}
\title{{biGaussianNull}}
\usage{
biGaussianNull(
  rstat,
  null.quantile = c(0.25, 0.75),
  method = c("biGaussianMean0", "biGaussianModeSplit")
)
}
\arguments{
\item{rstat}{The observed test statistics.}

\item{null.quantile}{Data quantiles, pre-fixed so that we could reduce the impact of extreme values for estimating prior null distribution. Default value is c(0.25, 0.75).}

\item{method}{A char. The methods we used to do estimation: either "biGaussian"-- EM algorithm for mixtures of two univariate normals, or "biGaussianWith0Cutoff"-- assume all negative test statistics forms one normal and all positive test statistics forms the other one normal. And proportion parameter is proportional to a number of observations each class, or "biGaussianMean0"-- null is formed by two half normals. "biGaussianModeSplit"-- split data from median value, then flip each part to the other side to estimate normal distribution.}
}
\value{
A list with element
 \item{alpha}{cutoff value (if given)}
 \item{mu}{mean positon for two normals}
 \item{sd}{sd value for two normals}
 \item{prop}{proportion for each normal}
}
\description{
Compute a prior null distribution without given any prior biological knowledge. Methods used to estimate this prior include: "biGaussian","biGaussianWith0Cutoff","biGaussianMean0","biGaussianModeSplit" etc.
}
\keyword{internal}
