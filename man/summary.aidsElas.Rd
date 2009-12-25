\name{summary.aidsElas}
\alias{summary.aidsElas}
\alias{print.summary.aidsElas}
\title{Summarizing the Elasticities of an Almost Ideal Demand System}

\description{
   These functions summarize and print the estimated elasticities
   of an Almost Ideal Demand System (AIDS).
}

\usage{
\method{summary}{aidsElas}( object, \dots )

\method{print}{summary.aidsElas}( x, \dots )
}

\arguments{
   \item{object}{an object of class \code{aidsElas}.}
   \item{x}{an object of class \code{summary.aidsElas}.}
   \item{\dots}{currently ignored.}
}

\value{
   \code{summary.aidsEst} returns a list of class \code{summary.aidsElas}.
   It is identical to the provided \code{object} (except for its class),
   but it contains following additional element:
\item{table}{a matrix with 4 columns: all elasticities,
   their standard errors (if available), their t-values (if available),
   and their P-values (if available).}
}

\seealso{\code{\link{aidsElas}}, \code{\link{aidsEst}}.}

\author{Arne Henningsen}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   ## Repeating the demand analysis of Blanciforti, Green & King (1986)
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "IL" )
   summary( elas( estResult ) )
}

\keyword{models}
