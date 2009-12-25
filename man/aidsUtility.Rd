\name{aidsUtility}
\alias{aidsUtility}
\alias{aidsUtilityDeriv}
\title{Indirect Utility Function of the Almost Ideal Demand System}

\description{
   These functions calculate the utility level
   given prices and total expenditure
   using the indirect utility function
   of the Almost Ideal Demand System
   and the partial derivatives of this indirect utility function
   with repect to prices and total expenditure.
}

\usage{
aidsUtility( priceNames, totExpName, coef, data )

aidsUtilityDeriv( priceNames, totExpName, coef, data, rel = FALSE )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{coef}{a list containing the coefficients in elements
      alpha0 (scalar), alpha (vector), beta (vector), gamma (matrix),
      and possibly beta0 (scalar, if not given, it is assumed to be 1).}
   \item{data}{a data frame containing the data.}
   \item{rel}{logical. If \code{TRUE} the returned partial derivatives are
      given in relative terms (like elasticities),
      i.e.\ they indicate the percentage change in the utility level
      when a price or total expenditure is increased by 1\%.}
}

\value{
   \code{aidsUtility} returns a numeric vector
   that contains the utility levels;
   \code{aidsUtilityDeriv} returns a data.frame
   that contains the partial derivatives
   of the indirect utility function
   with repect to prices and total expenditure.
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.
}

\author{Arne Henningsen}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsCalc}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   ## estimate the (non-linear) AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )

   # calculate the utility levels of each year
   utility <- aidsUtility( priceNames, "xFood", coef = coef( estResult ),
      data = Blanciforti86 )

   utilityDeriv <- aidsUtilityDeriv( priceNames, "xFood",
      coef = coef( estResult ), data = Blanciforti86 )

   utilityEla <- aidsUtilityDeriv( priceNames, "xFood",
      coef = coef( estResult ), data = Blanciforti86, rel = TRUE )
}

\keyword{ models }
