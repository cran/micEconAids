\name{aidsCalc}
\alias{aidsCalc}
\alias{predict.aidsEst}
\title{ Shares and Quantities of the Almost Ideal Demand System }

\description{
  Given prices, total expenditure and coefficients this function
  calculates the demanded quantities and expenditure shares based
  on the Almost Ideal Demand System.
}

\usage{
aidsCalc( priceNames, totExpName, coef, data, priceIndex = "TL",
   basePrices = NULL, baseShares = NULL, shifterNames = NULL )

\method{predict}{aidsEst}( object, newdata = NULL,
   observedShares = FALSE, \dots )

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
   \item{shifterNames}{a vector of strings containing the names 
      of the demand shifters.}
   \item{object}{an object of class \code{aidsEst}.}
   \item{newdata}{an optional data frame which should contain the variables
      for the prediction.
      If omitted, the data frame used for the estimation is used
      also for the prediction.}
   \item{observedShares}{logical. Using observed shares? (see details).}
   \item{\dots}{currently not used.}
}

\details{
   The \code{predict} method for objects of class \code{aidsEst}
   extracts all relevant elements from an object returned by \code{aidsEst}
   and passes them as arguments to \code{aidsCalc}.
   The optional argument \code{observedShares} determines
   whether fitted (default) or observed expenditure shares are
   used in the price index of the LA-AIDS.
}

\value{
   \code{aidsCalc} and
   the \code{predict} method for objects of class \code{aidsEst}
   return a list with following elements:
   \item{shares}{a data frame containing the calculated expenditure shares.}
   \item{quantities}{a data frame containing the calculated quantites.}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.
}

\author{Arne Henningsen}

\seealso{ \code{\link{aidsEst}}, \code{\link{aidsPx}} }

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   ## LA-AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, priceIndex = "S" )

   # using observed shares in the Stone index
   lnp <- aidsPx( "S", priceNames, Blanciforti86, shareNames )
   fitted <- aidsCalc( priceNames, "xFood", coef = coef( estResult ),
      data = Blanciforti86, priceIndex = lnp )
   fitted$shares  # equal to estResult$wFitted
   fitted$quant   # equal to estResult$qFitted
   # now the same with the predict method
   fitted2 <- predict( estResult, observedShares = TRUE )
   all.equal( fitted, fitted2 ) 

   # using fitted shares in the Stone index
   fitted <- aidsCalc( priceNames, "xFood", coef = estResult$coef,
      data = Blanciforti86, priceIndex = "S" )
   # now the same with the predict method
   fitted2 <- predict( estResult )
   all.equal( fitted, fitted2 ) 

   ## AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )

   fitted <- aidsCalc( priceNames, "xFood", coef = coef( estResult ),
      data = Blanciforti86 )

   fitted$shares  # equal to estResult$wFitted
   fitted$quant   # equal to estResult$qFitted

   fitted2 <- predict( estResult )
   all.equal( fitted, fitted2 ) 
}

\keyword{ models }
