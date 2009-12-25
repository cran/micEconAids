\name{vcov.aidsEst}
\alias{vcov.aidsEst}

\title{Covariance matrix of an Almost Ideal Demand System}

\description{
   These method returns the covariance matrix of the coefficients from
   an Almost Ideal Demand System (AIDS).
}

\usage{
   \method{vcov}{aidsEst}( object, \dots )
}

\arguments{
  \item{object}{an object of class \code{aidsEst}.}
  \item{\dots}{currently not used.}
}

\value{
The \code{vcov} method for objects of class \code{aidsEst}
returns a symmetric matrix: the covariance matrix of the coefficients.
}

\author{Arne Henningsen}

\seealso{\code{\link{aidsEst}}, \code{\link{coef.aidsEst}}}

\keyword{methods}
