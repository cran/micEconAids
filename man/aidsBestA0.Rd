\name{aidsBestA0}
\alias{aidsBestA0}
\title{Find 'best' Value for alpha 0 in the AIDS}

\description{
   Search for the intercept of the translog price index
   (\eqn{\alpha_0}) that gives the best fit to the
   Almost Ideal Demand System (see Michalek and Keyzer, 1992)
}

\usage{
   aidsBestA0( priceNames, shareNames, totExpName,
      a0min = -50, a0max = 50, stoprange = 3, stopiter = 10,
      verbose = FALSE, \dots )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{a0min}{lower limit of the range for searching for \eqn{\alpha_0}.}
   \item{a0max}{upper limit of the range for searching for \eqn{\alpha_0}.}
   \item{stoprange}{stop searching when the search interval is smaller than
      or equal to \code{stoprange}.}
   \item{stopiter}{maximal number of iterations.}
   \item{verbose}{print each determinant of the residual covariance matrix
      immediately after its calculation.}
   \item{\dots}{arguments passed to \code{\link{aidsEst}}.}
}

\details{
   The demand system is estimated using the Iterative Linear Least Squares
   Estimator (ILLE) suggested by Blundell and Robin (1999).
   This iterative procedure is equivalent to the method proposed by
   Michalek and Keyzer (1992). However, the latter do not correct
   the coefficient covariance matrix.

   The fit of the model is measured in terms of the likelihood value.
   Since the determinant of the residual covariance matrix is monotonically
   decreasing with the likelihood value, we search for the smallest
   determinant of the residual covariance matrix.

   Since each call of \code{\link{aidsEst}} generally takes a long time,
   the search algorithm is constructed to minimize the
   calls of the function \code{\link{aidsEst}}.
}

\value{
   a list containing following objects:
  \item{alpha0}{\eqn{\alpha_0} that gives the best fit.}
  \item{allValues}{all \eqn{\alpha_0} values that have been tested and
      the determinants of the corresponding residual covariance matrices.}
  \item{iter}{number of iterations.}
}

\references{
   Blundell, R. and J.M. Robin (1999)
   Estimationin Large and Disaggregated Demand Systems:
   An Estimator for Conditionally Linear Systems.
   \emph{Journal of Applied Econometrics}, 14, p. 209-232.

   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Michalek, J. and M. A. Keyzer (1992)
   Estimation of a two-stage LES-AIDS consumer demand system
   for eight EC countries.
   \emph{European Review of Agricultural Economics}, 19 (2), p. 137-163.
}

\author{Arne Henningsen}

\seealso{ \code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   bestA0 <- aidsBestA0( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, useMatrix = FALSE )
   # may take some time (argument 'useMatrix = FALSE' decreases
   # the computation time only if the model and data set are small)
   print( bestA0$alpha0 )
   plot( bestA0$allValues ) # this should be convex
}

\keyword{ models }
