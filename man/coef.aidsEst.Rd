\name{coef.aidsEst}
\alias{coef.aidsEst}
\alias{print.coef.aidsEst}

\title{Coefficients of an Almost Ideal Demand System}

\description{
   These methods return and print the coefficients from
   an Almost Ideal Demand System.
}

\usage{
   \method{coef}{aidsEst}( object, \dots )

   \method{print}{coef.aidsEst}( x, \dots )
}

\arguments{
  \item{object}{an object of class \code{aidsEst}.}
  \item{x}{an object of class \code{coef.aidsEst}.}
  \item{\dots}{further arguments for methods}
}

\value{
The \code{coef} method returns an object of class \code{coef.aidsEst}
containing following objects:
   \item{alpha0}{a scalar, coefficient alpha0
      (only for the AIDS with translog price index).}
   \item{alpha}{a vector of the alpha coefficients.}
   \item{beta}{a vector of the beta coefficients.}
   \item{gamma}{a matrix of the gamma coefficients.}
   \item{delta}{a matrix of the delta coefficients
      (only if the model was estimated with demand shifters).}
}

\author{Arne Henningsen}

\seealso{\code{\link{aidsEst}}}

\keyword{methods}
