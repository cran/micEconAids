\name{summary.aidsEst}
\alias{summary.aidsEst}
\alias{print.summary.aidsEst}
\title{Summarizing the Estimation of an Almost Ideal Demand System}

\description{
   \code{summary.aidsEst} summarizes the estimation results
   of an Almost Ideal Demand System (AIDS).
}

\usage{
\method{summary}{aidsEst}( object, \dots )

\method{print}{summary.aidsEst}( x, \dots )
}

\arguments{
   \item{object}{an object of class \code{aidsEst}.}
   \item{x}{an object of class \code{summary.aidsEst}.}
   \item{\dots}{currently ignored.}
}

\value{
   \code{summary.aidsEst} returns a list of class \code{summary.aidsEst}
   that is currently identical to the provided \code{object}
   (except for its class).
}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsElas}}.}

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
   print( summary( estResult ) )
}

\keyword{models}
