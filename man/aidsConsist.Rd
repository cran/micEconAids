\name{aidsConsist}
\alias{aidsConsist}
\alias{print.aidsConsist}
\alias{checkConsist.aidsEst}
\title{Check Consistency of the AIDS}
\description{
   Check whether the specification of the AIDS is consistent
   with microeconomic demand theory (i.e. utility maximisation).
}

\usage{
aidsConsist( priceNames, totExpName, coef, data,
   priceIndex = "TL", basePrices = NULL, baseShares = NULL,
   shareNames = NULL )

\method{print}{aidsConsist}( x, \dots )

\method{checkConsist}{aidsEst}( object, observedShares = FALSE, \dots )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{coef}{a list containing the coefficients alpha, beta, gamma,
      and (only for the translog price index) alpha0.}
   \item{data}{a data frame containing the data.}
   \item{priceIndex}{a character string specifying the price index
      (see \code{\link{aidsPx}})
      or a numeric vector providing the log values of the price index
      (passed to \code{\link{aidsMono}}.}
   \item{basePrices}{a vector specifying the base prices for the
      Paasche, Laspeyres, and Tornqvist price index
      (passed to \code{\link{aidsMono}}.}
   \item{baseShares}{a vector specifying the base expenditure shares for the
      Laspeyres, simplified Laspeyres, and Tornqvist index
      (passed to \code{\link{aidsMono}}.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares (passed to \code{\link{aidsConcav}}.}
   \item{x}{an object of class \code{aidsConsist}.}
   \item{object}{an object of class \code{aidsEst}.}
   \item{observedShares}{logical. Using observed shares? (see details).}
   \item{\dots}{currently not used.}
}

\details{
   The \code{checkConsist} method for objects of class \code{aidsEst}
   extracts all relevant elements from an object returned by \code{aidsEst}
   and passes them as arguments to \code{aidsConsist}.
   The optional argument \code{observedShares} determines
   whether fitted (default) or observed expenditure shares are
   used in the price index of the LA-AIDS to check monotonicity
   and in the calculation of the subsitution matrix to check concavity.
}

\value{
   \code{aidsConsist} and the \code{checkConsist} method
   return a list of class \code{aidsConsist}
   that contains following elements:
   \item{addingUp}{logical. Is the addinp-up condition fulfilled?}
   \item{homogeneity}{logical. Is the homogeneity condition fulfilled?}
   \item{symmetry}{logical. Is the symmetry condition fulfilled?}
   \item{mono}{an object returned by \code{\link{aidsMono}}.}
   \item{concav}{an object returned by \code{\link{aidsMono}}
      (only if the symmetry condition is fulfilled and
      the translog price index is used,
      i.e. argument \code{priceIndex} is \code{"TL"}.}
}
\references{
   Deaton, A.S. and J. Muellbauer (1980a)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Deaton, A.S. and J. Muellbauer (1980b)
   \emph{Economics and Consumer Behavior},
   Cambridge University Press, Cambridge.
}

\author{Arne Henningsen}

\seealso{\code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )
   aidsConsist( priceNames, "xFood", data = Blanciforti86,
      coef = estResult$coef )

   # the same can be obtained in an easier way
   checkConsist( estResult )
}

\keyword{models}
