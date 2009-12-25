\name{aidsMono}
\alias{aidsMono}
\alias{print.aidsMono}
\title{Monotonicity of the AIDS}
\description{
   \code{aidsMono} checks whether the expenditure function
   of an estimated Almost Ideal Demand System (AIDS)
   is monotonic increasing in prices,
   which implies that all demanded quantities and expenditure shares
   are non-negative.
}

\usage{
aidsMono( priceNames, totExpName, coef, data,
   priceIndex = "TL", basePrices = NULL, baseShares = NULL )

\method{print}{aidsMono}( x, header = TRUE, \dots )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{coef}{a list containing the coefficients alpha, beta, gamma,
      and (only for the translog price index) alpha0.}
   \item{data}{a data frame containing the data.}
   \item{priceIndex}{a character string specifying the price index
      (see \code{\link{aidsPx}})
      or a numeric vector providing the log values of the price index.}
   \item{basePrices}{a vector specifying the base prices for the
      Paasche, Laspeyres, and Tornqvist price index.}
   \item{baseShares}{a vector specifying the base expenditure shares for the
      Laspeyres, simplified Laspeyres, and Tornqvist index.}
   \item{x}{an object of class \code{aidsMono}.}
   \item{header}{logical. Print a header?}
   \item{\dots}{currently unused.}
}

\details{
   Internally, \code{aidsMono} passes its arguments to \code{aidsCalc}
   and then checks for each observation, wether all expenditure shares
   are non-negative.

   If argument \code{priceIndex} specifies a price index of the LA-AIDS,
   'fitted' values are used for current and lagged expenditure shares
   in these price indices
   (see \code{\link{aidsCalc}}).
   However, if argument \code{priceIndex} is a numeric vector
   containing the log values of a price index
   (e.g.\ the price index used in the estimation),
   this price index is used for the calculations.
}

\value{
   \code{aidsMono} returns a list of class \code{aidsMono}
   that contains following elements:
   \item{monotony}{a logical vector indicating whether the monotony
      condition is fulfilled at each observation.}
   \item{nValidObs}{number of (valid) observation at which monotonicity
      could be checked.}
   \item{nMonoObs}{number of observation at which the monotonicity codition
      is fulfilled.}
   \item{monoPercent}{percent of observations where the monotony condition is
      fulfilled.}
   \item{priceIndex}{a chacter string indicating the price index
      specified by argument \code{priceIndex}
      (\code{"numeric"} if the price index is specified numerically).}
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

\seealso{\code{\link{aidsEst}}, \code{\link{aidsCalc}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   ## AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )
   aidsMono( priceNames, "xFood", coef = coef( estResult ),
      data = Blanciforti86 )

   ## LA-AIDS with Tornqvist price index
   estResultLaT <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, priceIndex = "T" )
   # with fitted expenditure shares in the price index
   aidsMono( priceNames, "xFood", coef = coef( estResultLaT ),
      data = Blanciforti86, priceIndex = "T",
      basePrices = estResultLaT$basePrices,
      baseShares = estResultLaT$baseShares )
   # with observed expenditure shares in the price index
   aidsMono( priceNames, "xFood", coef = coef( estResultLaT ),
      data = Blanciforti86, priceIndex = estResultLaT$lnp )
}

\keyword{models}
