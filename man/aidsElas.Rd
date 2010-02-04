\name{aidsElas}
\alias{aidsElas}
\alias{elas.aidsEst}
\alias{print.aidsElas}
\title{ Elasticities of the AIDS model }

\description{
   These functions calculate and print
   the demand elasticities of an AIDS model.
}

\usage{
   aidsElas( coef, prices = NULL, shares = NULL, totExp = NULL,
      method = "AIDS", priceIndex = "TL", basePrices = NULL, baseShares = NULL,
      quantNames = NULL, priceNames = NULL, coefCov = NULL, df = NULL )

   \method{elas}{aidsEst}( object, method = NULL, observedShares = FALSE, \dots )

   \method{print}{aidsElas}( x, \dots )
}

\arguments{
   \item{coef}{a list containing the coefficients alpha, beta and gamma.}
   \item{prices}{a vector of the prices at which the elasticities should be calculated.}
   \item{shares}{a vector of the shares at which the elasticities should be calculated.}
   \item{totExp}{total expenditure at which the elasticities should be calculated.}
   \item{method}{the elasticity formula to be used (see details).}
   \item{priceIndex}{the price index (see details).}
   \item{basePrices}{a vector specifying the base prices for the
      Paasche, Laspeyres, and Tornqvist price index.}
   \item{baseShares}{a vector specifying the base expenditure shares for the
      Laspeyres, simplified Laspeyres, and Tornqvist index.}
   \item{quantNames}{an optional vector of strings containing the names of
      the quantities to label elasticities.}
   \item{priceNames}{an optional vector of strings containing the names of
      the prices to label elasticities.}
   \item{coefCov}{variance covariance matrix of the coefficients (optional).}
   \item{df}{degrees of freedom to calculate P-values of the elasticities 
      (optional).}
   \item{object}{an object of class \code{aidsEst}.}
   \item{observedShares}{logical. Using observed shares for calculating the
      demand elasticities?}
   \item{x}{an object of class \code{aidsElas}.}
   \item{\dots}{additional arguments of \code{elas.aidsEst}
      are passed to \code{aidsEla};
      additional arguments of \code{print.aidsElas}
      are currently ignored.}
}

\details{
   Currently, \code{aidsElas} and \code{elas.aidsEst} can calculate
   elasticities only for models without demand shifters.
   However, the user can calculate elasticies for models with demand shifters
   by removing the coefficients of the demand shifters
   (\eqn{\delta_i}{delta_i}, \code{coef$delta}),
   adjusting the coefficients \eqn{\alpha_i}{alpha_i} (\code{coef$alpha})
   \sQuote{by hand},
   and then calling \code{aidsElas}.

   Argument \code{priceIndex} has two effects:
   first it determines the price index that is used
   for calculating (fitted) expenditure shares,
   if argument \code{shares} is not provided (see \code{\link{aidsCalc}});
   second it determines which version of the formulas for calculating
   demand elasticities of the LA-AIDS are used,
   because formulas \code{B1}/\code{LA}, \code{B2}, and \code{Go}/\code{Ch}
   have different versions depending on the price index.

   \code{elas.aidsEst} is a wrapper function to \code{aidsElas}
   that extracts the
   estimated coefficients (\code{coef}),
   mean expenditure shares (\code{wMeans}),
   mean prices (\code{pMeans}),
   names of the prices (\code{priceNames}),
   estimated coefficient variance covariance matrix (\code{coef$allcov}), and
   degrees of freedom (\code{est$df})
   from the object of class \code{aidsEst}
   and passes them to \code{aidsElas}.
   If argument \code{method} in \code{elas.aidsEst} is not specified,
   the default value depends on the estimation method.
   If the demand system was estimated by the linear approximation (LA),
   the default method is 'Ch'.
   If the demand system was estimated by the iterative linear least squares
   estimator (ILLE),
   the default method is 'AIDS'.

   At the moment the elasticity formulas of the orginal AIDS (\code{AIDS}),
   the formula of Goddard (1983) or Chalfant (1987) (\code{Go} or \code{Ch}),
   the formula of Eales and Unnevehr (1988) (\code{EU}),
   the formula of Green and Alston (1990) or the first of Buse (1994)
   (\code{GA} or \code{B1}) and
   the second formula of Buse (1994) (\code{B2})
   are implemented.

   The variance covariance matrices of the elasticities are calculated using
   the formula of Klein (1953, p. 258) (also known as the delta method).
   At the moment this is implemented only for the elasticity formulas of the 
   orginal AIDS.
}

\value{
  a list of class \code{aidsElas} containing following elements:
  \item{method}{the elasticity formula used to calculate these elasticities.}
  \item{priceIndex}{the price index used (see details).}
  \item{df}{degrees of freedom to calculate P-values of the elasticities
    (only if argument \code{df} is provided).}
  \item{exp}{vector of expenditure elasticities.}
  \item{hicks}{matrix of Hicksian (compensated) price elasticities.}
  \item{marshall}{matrix of Marshallian (uncompensated) price elasticities.}
  \item{allVcov}{variance covariance matrix of all elasticities.}
  \item{expVcov}{variance covariance matrix of the expenditure elasticities.}
  \item{hicksVcov}{variance covariance matrix of the Hicksian (compensated) 
    price elasticities.}
  \item{marshallVcov}{variance covariance matrix of the Marshallian 
    (uncompensated) price elasticities.}
  \item{expStEr}{standard errors of the expenditure elasticities.}
  \item{hicksStEr}{standard errors of the Hicksian (compensated) price 
    elasticities.}
  \item{marshallStEr}{standard errors of the Marshallian (uncompensated) 
    price elasticities.}
  \item{expTval}{t-values of the expenditure elasticities.}
  \item{hicksTval}{t-values of the Hicksian (compensated) price elasticities.}
  \item{marshallTval}{t-values of the Marshallian (uncompensated) price 
    elasticities.}
  \item{expPval}{P-values of the expenditure elasticities.}
  \item{hicksPval}{P-values of the Hicksian (compensated) price elasticities.}
  \item{marshallPval}{P-values of the Marshallian (uncompensated) price 
    elasticities.}
}

\references{
   Chalfant, J.A. (1987)
   A Globally Flexible, Almost Ideal Demand System.
   \emph{Journal of Business and Economic Statistics}, 5, p. 233-242.

   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Eales J.S. and L.J. Unnevehr (1988)
   Demand for beef and chicken products: separability and structural change.
   \emph{American Journal of Agricultural Economics}, 70, p. 521-532.
   
   Klein L.R. (1953)
   \emph{A Textbook of Econometrics}. Row, Petersen and Co., New York.
}

\author{Arne Henningsen}

\seealso{\code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86 )
   wMeans <- colMeans( Blanciforti86[ , c( "wFood1", "wFood2",
      "wFood3", "wFood4" ) ] )
   aidsElas( estResult$coef, shares = wMeans, method = "Ch",
      priceIndex = "S" )

   ## Repeating the evaluation of different elasticity formulas of
   ## Green & Alston (1990)
   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   # AIDS estimation and elasticities
   estResultA <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86[ -1, ],
      method = "IL", maxiter = 100 )
   diag( elas( estResultA, method = "AIDS" )$marshall )

   # LA-AIDS estimation
   estResultLA <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, priceIndex = "SL", maxiter = 100 )

   # LA-AIDS + formula of AIDS
   diag( elas( estResultLA, method = "AIDS" )$marshall )

   # LA-AIDS + formula of Eales + Unnevehr
   diag( elas( estResultLA, method = "EU" )$marshall )

   # LA-AIDS + formula of Goddard or Chalfant:
   diag( elas( estResultLA, method = "Go" )$marshall )
   diag( elas( estResultLA, method = "Ch" )$marshall )

   # LA-AIDS + formula of Green + Alston (= 1st of Buse):
   diag( elas( estResultLA, method = "GA" )$marshall )
}

\keyword{ models }
