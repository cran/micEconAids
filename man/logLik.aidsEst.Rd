\name{logLik.aidsEst}
\alias{logLik.aidsEst}
\title{Log-Likelihood value of an object of class aidsEst}

\description{
   This method extracts the log-likelihood value
   of a fitted Almost Ideal Demand System (AIDS)
}

\usage{
\method{logLik}{aidsEst}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{aidsEst}.}
   \item{\dots}{currently ignored.}
}

\value{
   A numeric scalar (the log-likelihood value) with 2 attributes:
   \code{nobs} (total number of observations in all equations) and
   \code{df} (number of free parameters, i.e. coefficients
   + elements of the residual covariance matrix).
}

\seealso{\code{\link{aidsEst}}.}

\author{Arne Henningsen}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   ## Repeating the demand analysis of Blanciforti, Green & King (1986)
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, priceIndex = "SL" )
   logLik( estResult )
}

\keyword{models}
