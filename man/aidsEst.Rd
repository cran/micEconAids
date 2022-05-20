\name{aidsEst}
\alias{aidsEst}
\alias{print.aidsEst}
\title{Estimating the Almost Ideal Demand System (AIDS) }

\description{
   \code{aidsEst} does an econometric estimation
   of the Almost Ideal Demand System (AIDS)
}

\usage{
aidsEst( priceNames, shareNames, totExpName, data,
      method = "LA", priceIndex = "Ls", pxBase = 1,
      hom = TRUE, sym = TRUE,
      shifterNames = NULL, instNames = NULL,
      estMethod = ifelse( is.null( instNames ), "SUR", "3SLS" ),
      ILmaxiter = 50, ILtol = 1e-5, alpha0 = 0, restrict.regMat = FALSE, ... )

\method{print}{aidsEst}( x, ... )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{data}{a data frame containing all required variables.}
   \item{method}{character string specifying the method to estimate the AIDS:
      either 'LA' or 'IL' (see deatils).}
   \item{priceIndex}{character string specifying the price index
      for the 'Linear Approximation':
      either 'S', 'SL', 'P', 'L', 'Ls', or 'T' (see details).}
   \item{pxBase}{The base to calculate the LA-AIDS price indices
      (see \code{\link{aidsPx}}).}
   \item{hom}{logical. Should the homogeneity condition be imposed?}
   \item{sym}{logical. Should the symmetry condition be imposed?}
   \item{shifterNames}{an optional vector of strings containing the names of
      the demand shifters.}
   \item{instNames}{a vector of strings containing the names of instrumental
      variables.}
   \item{estMethod}{estimation method (e.g. 'SUR' or '3SLS',
      see \code{\link[systemfit]{systemfit}}).}
   \item{ILmaxiter}{maximum number of iterations of the
      'Iterated Linear Least Squares Estimation'.}
   \item{ILtol}{tolerance level of the 'Iterated Linear Least Squares
      Estimation'.}
   \item{alpha0}{the intercept of the translog price index (\eqn{\alpha_0}).}
   \item{restrict.regMat}{logical. Method to impose homogeneity and symmetry restrictions:
      either via restrict.matrix (default) or via restrict.regMat
      (see \code{\link[systemfit]{systemfit}}).}
   \item{x}{An object of class \code{aidsEst}.}
   \item{...}{additional arguments of \code{aidsEst} are passed to
      \code{\link[systemfit]{systemfit}};
      additional arguments of \code{print.aidsEst} are currently ignored.}
}

\details{
   Argument \code{method} can specify two different estimation methods:
   The 'Linear Approximate AIDS' (LA) and the 'Iterative Linear Least Squares
   Estimator' (IL) proposed by Blundell and Robin (1999).\cr
   Argument \code{priceIndex} can specify six different price indices
   for the LA-AIDS:
   \itemize{
   \item Stone price index ('S'),
   \item Stone price index with lagged shares ('SL'),
   \item loglinear analogue to the Paasche price index ('P'),
   \item loglinear analogue of the Laspeyres price index ('L'),
   \item simplified loglinear analogue of the Laspeyres price index ('Ls'), and
   \item Tornqvist price index ('T').
   }

   The 'Iterative Linear Least Squares Estimator' (IL) needs starting
   values for the (translog) price index.
   Starting values are taken from an initial estimation
   of the 'Linear Approximate AIDS' (LA) with the price index
   specified by argument \code{priceIndex}.
}

\value{
   a list of class \code{aidsEst} containing following objects:
   \item{coef}{a list containing the vectors/matrix of the estimated
      coefficients (alpha, beta, and gamma).}
   \item{r2}{\eqn{R^2}-values of all share equations.}
   \item{r2q}{\eqn{R^2}-values of the estimated quantities.}
   \item{wFitted}{fitted expenditure shares.}
   \item{wResid}{residuals of the expenditure shares.}
   \item{qObs}{observed quantities / quantitiy indices.}
   \item{qFitted}{fitted quantities / quantitiy indices.}
   \item{qResid}{residuals of the estimated quantities.}
   \item{est}{estimation result, i.e. the object returned
      by \code{\link[systemfit]{systemfit}}.}
   \item{iter}{iterations of SUR/3SLS estimation(s).
      If the AIDS is estimated by the 'Iterated Linear Least Squares
      Estimator' (ILLE):
      a vector containing the SUR/3SLS iterations at each iteration.}
   \item{ILiter}{number of iterations of the 'Iterated Linear Least Squares
      Estimation'.}
   \item{method}{the method used to estimate the aids (see details).}
   \item{priceIndex}{the name of the price index (see details).}
   \item{lnp}{log of the price index used for estimation.}
   \item{pMeans}{means of the prices.}
   \item{wMeans}{means of the expenditure shares.}
   \item{xtMean}{mean of total expenditure.}
   \item{sMeans}{means of shifter variables (only if shifter variables are used).}
   \item{call}{the call of \code{aidsEst}.}
   \item{priceNames}{names of the prices.}
   \item{shareNames}{names of the expenditure shares.}
   \item{totExpName}{name of the variable for total expenditure.}
   \item{basePrices}{the base prices of the Paasche, Laspeyres,
      or Tornqvist price index.}
   \item{baseShares}{the base shares of the Laspeyres, simplified Laspeyres,
      or Tornqvist price index.}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Blundell, R. and J.M. Robin (1999)
   Estimationin Large and Disaggregated Demand Systems:
   An Estimator for Conditionally Linear Systems.
   \emph{Journal of Applied Econometrics}, 14, p. 209-232.
}

\seealso{\code{\link{summary.aidsEst}}, \code{\link{aidsElas}},
   \code{\link{aidsCalc}}.}

\author{Arne Henningsen}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   ## Repeating the demand analysis of Blanciforti, Green & King (1986)
   ## Note: Blanciforti, Green & King (1986) use scaled data,
   ##       which leads to slightly different results
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, priceIndex = "SL", maxiter = 100 )
   print( estResult )
   elas( estResult )

   ## Estimations with a demand shifter: linear trend
   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
   Blanciforti86$trend <- c( 0:( nrow( Blanciforti86 ) - 1 ) )
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, shifterNames = "trend" )
   print( estResult )

   # Estimations with two demand shifters: linear + quadratic trend
   Blanciforti86$trend2 <- c( 0:( nrow( Blanciforti86 ) - 1 ) )^2
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, shifterNames = c( "trend", "trend2" ) )
   print( estResult )
}

\keyword{models}
