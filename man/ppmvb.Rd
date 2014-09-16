\name{ppmvb}
\alias{ppmvb}
\title{Bayesian ppm-function}
\usage{
ppmvb(Q, ...)
}
\arguments{
  \item{...}{Parameters to ppm. Method will be "logi".}
}
\description{
Otherwise the same as \code{\link{ppm}} with method="logi"
but using a Bayesian fitter.
}
\details{
The variational approximation uses the package
\code{\link{vblogistic}}. 'verb' and 'eps' parameters are
renamed here 'verbosity' and 'epsilon' as ppm uses the
former internally.

The returned object is apart from $internal$glmfit the same
as that of ppm. Some of the methods (e.g. plot) wont work
at the moment.

Method \link{vbsummary} provides some extra Bayesian
information.

The $internal$glmfit is of class "vblogitfit" with print,
marginals and plot summary.

Prior arguments:

* m0: Gaussian prior mean vector.

* S0: Gaussian prior covariance matrix.

Make sure the dimensions match, e.g. try basic ppm and look
how many parameters are estimated.
}
\examples{
library(spatstat)
x <- rStrauss(100, 0.1, 0.07)
# glm:
f0 <- ppm(x, interaction=Strauss(0.06), method="logi")
# vb-logistic, with bad prior:
f1 <- ppmvb(x, interaction=Strauss(0.06), m0=log(c(100, 0.3)), S0=diag(c(10, 3)), verbose=TRUE )

print( exp( rbind(coef(f0), coef(f1))) )

summary(f1)



s <- vbsummary(f1)
# this is of class vblogitfit
print(s)
plot(s)
# in exp scale
plot(s, log=FALSE)
}

