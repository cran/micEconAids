\name{fitted.aidsEst}
\alias{fitted.aidsEst}

\title{Fitted values of an Almost Ideal Demand System}

\description{
   This method extracts the fitted demanded quantities and expenditure shares
   from an estimated Almost Ideal Demand System.
}

\usage{
   \method{fitted}{aidsEst}( object, \dots )
}

\arguments{
  \item{object}{an object of class \code{aidsEst}.}
  \item{\dots}{currently unused.}
}

\value{
The \code{fitted} method returns a list
containing following objects:
   \item{shares}{a data frame for the fitted expenditure shares.}
   \item{quant}{a data frame for the fitted demanded quantities.}
}

\author{Arne Henningsen}

\seealso{\code{\link{aidsEst}}}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86 )
   fitted( estResult )
}

\keyword{methods}
