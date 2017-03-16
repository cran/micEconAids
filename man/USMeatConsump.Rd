\name{USMeatConsump}
\docType{data}
\alias{USMeatConsump}

\title{U.S. Meat Consumption Data}

\description{
   The \code{USMeatConsump} data set contains quarterly retail prices
   and consumption quantities for four meat product categories:
   beef, pork, chicken, and turkey.
   The data period ranges from the first quarter of 1975
   to the third quarter of 1999.
   Hence, there are 99 observations.
}

\usage{data(USMeatConsump)}

\format{
   This data frame contains the following columns:
   \describe{
      \item{year}{Year.}
      \item{qtr}{Quarter of the year.}
      \item{t}{Time trend.}
      \item{pop}{Population [million].}
      \item{cpi}{Consumer price index.}
      \item{total_exp}{Total per capita expenditure.}
      \item{meat_exp}{Per capita expenditure on meat.}
      \item{beef_q}{Per capita consumption of beef [pound].}
      \item{pork_q}{Per capita consumption of pork[pound].}
      \item{chick_q}{Per capita consumption of chicken [pound].}
      \item{turkey_q}{Per capita consumption of turkey [pound].}
      \item{beef_p}{Retail price of beef [cents / pound].}
      \item{pork_p}{Retail price of pork [cents / pound].}
      \item{chick_p}{Retail price of chicken [cents / pound].}
      \item{turkey_p}{Retail price of turkey [cents / pound].}
      \item{beef_w}{Expenditure share of beef (in meat).}
      \item{pork_w}{Expenditure share of pork (in meat).}
      \item{chick_w}{Expenditure share of chicken (in meat).}
      \item{turkey_w}{Expenditure share of turkey (in meat).}
   }
}

\source{
   SAS,
   \emph{SAS/ETS Examples: Estimating an Almost Ideal Demand System Model},
   \url{https://support.sas.com/rnd/app/ets/examples/aids/index.htm}.
}

\examples{
   ## replicating the LA-AIDS estimation of the SAS example
   # loading data set
   data( USMeatConsump )

   # adding shifter variables for modeling seasonal effects
   USMeatConsump$co1 <- cos( 1 / 2 * 3.14159 * USMeatConsump$t )
   USMeatConsump$si1 <- sin( 1 / 2 * 3.14159 * USMeatConsump$t )

   # Scaling prices by their means
   USMeatConsump$beef_pm <- USMeatConsump$beef_p / mean( USMeatConsump$beef_p )
   USMeatConsump$pork_pm <- USMeatConsump$pork_p / mean( USMeatConsump$pork_p )
   USMeatConsump$chick_pm <- USMeatConsump$chick_p / mean( USMeatConsump$chick_p )
   USMeatConsump$turkey_pm <- USMeatConsump$turkey_p / mean( USMeatConsump$turkey_p )

   # Estimation of the model
   meatModel <- aidsEst( c( "beef_pm", "pork_pm", "chick_pm", "turkey_pm" ),
      c( "beef_w", "pork_w", "chick_w", "turkey_w" ),
      "meat_exp", shifterNames = c( "co1", "si1", "t" ),
      priceIndex ="S", data = USMeatConsump, maxiter=1000 )
   summary( meatModel )
}


\keyword{datasets}
