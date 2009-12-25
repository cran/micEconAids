\name{aidsPx}
\alias{aidsPx}
\title{ Price Index for the AIDS }

\description{
   Calculate log price indices for the AIDS and LA-AIDS.
}

\usage{
   aidsPx( priceIndex, priceNames, data, shareNames = NULL, base = 1, coef = NULL,
      shifterNames = NULL )
}

\arguments{
   \item{priceIndex}{the price index to be used (see details).}
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{data}{a data frame containing the required variables.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{base}{the base to calculate the indices (see details).}
   \item{coef}{a list containing the coefficients alpha0, alpha, beta, and gamma
      (only needed for the translog price index).}
   \item{shifterNames}{an optional vector of strings containing the names of
      the demand shifters
      that modify the alphas of the Translog price index.}
}
\details{
   Currently 7 different price indices are implemented:
   \itemize{
   \item translog price index ('TL'),
   \item Stone price index ('S'),
   \item Stone price index with lagged shares ('SL'),
   \item loglinear analogue to the Paasche price index ('P'),
   \item loglinear analogue of the Laspeyres price index ('L'),
   \item simplified loglinear analogue of the Laspeyres price index ('Ls'), and
   \item Tornqvist price index ('T').
   }

   The argument \code{base} can be either
   \itemize{
   \item a single number: the row number of the base prices and quantities,
   \item a vector indicating several observations: The means of these
      observations are used as base prices and quantities, or
   \item a logical vector with the same length as the \code{data}: The
      means of the observations indicated as 'TRUE' are used as base prices
      and quantities.
   \item a list with elements \code{prices} and/or \code{shares}
      that are numeric vectors specifying the values of the base prices
      and/or base expenditure shares.
   }
}

\value{
   A vector containing the log price index.
   If argument \code{priceIndex} is \code{"P"}, \code{"Ls"}, or \code{"T"},
   the returned vector has the attribute(s)
   \code{basePrices} and/or \code{baseShares}
   that are numeric vectors of the base prices and/or base expenditure shares
   for calculating the price indices.
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Moschini, G. (1995)
   Units of Measurement and the Stone Index in Demand System Estimation.
   \emph{American Journal of Agricultural Economics}, 77, p. 63-68.
}

\author{Arne Henningsen}

\seealso{ \code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )

   # Stone price index
   aidsPx( "S", c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      Blanciforti86, c( "wFood1", "wFood2", "wFood3", "wFood4" ) )

   # (log-linear analogue of the) Paasche price index
   aidsPx( "P", c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      Blanciforti86, c( "wFood1", "wFood2", "wFood3", "wFood4" ),
      base = row.names(Blanciforti86) == "1970" )

   # Tornqvist price index
   aidsPx( "T", c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      Blanciforti86, c( "wFood1", "wFood2", "wFood3", "wFood4" ),
      base = list( prices = rep( 100, 4 ), shares = rep( 0.25, 4 ) ) )
}

\keyword{ models }
