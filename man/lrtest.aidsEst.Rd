\name{lrtest.aidsEst}
\alias{lrtest.aidsEst}

\title{Likelihood Ratio test for Almost Ideal Demand Systems}

\description{
   Testing hypothesis in Almost Ideal Demand Systems
   by a Likelihood Ratio test.
}

\usage{
   \method{lrtest}{aidsEst}( object, \dots )
}

\arguments{
   \item{object}{a fitted model object of class \code{aidsEst}.}
   \item{\dots}{further fitted model objects of class \code{aidsEst}.}
}

\details{
   \code{lrtest.aidsEst} consecutively compares
   the fitted model object \code{object}
   with the models passed in \code{\dots}.
}

\value{
   An object of class \code{anova},
   which contains the log-likelihood value,
   degrees of freedom, the difference in degrees of freedom,
   likelihood ratio Chi-squared statistic and corresponding p value.
   See documentation of \code{\link[lmtest]{lrtest}}
   in package "lmtest".
}

\author{Arne Henningsen}

\seealso{\code{\link{aidsEst}}, \code{\link[lmtest]{lrtest}}
   (package "lmtest"),
   \code{\link{logLik.aidsEst}}}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   # names of prices and expenditure shares
   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   # estimation with symmetry and homogeneity imposed
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86 )

   # estimation with only homogeneity imposed
   estResultHom <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, sym = FALSE )

   # unrestricted estimation
   estResultUnr <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, sym = FALSE, hom = FALSE )

   # LR tests
   lrtest( estResult, estResultHom, estResultUnr, estResult )

   ## Estimation with a demand shifter: linear trend
   Blanciforti86$trend <- c( 0:( nrow( Blanciforti86 ) - 1 ) )
   estResultTrend <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, shifterNames = "trend" )

   # LR tests
   lrtest( estResult, estResultTrend )
}

\keyword{models}

