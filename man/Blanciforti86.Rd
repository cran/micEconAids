\name{Blanciforti86}
\docType{data}
\alias{Blanciforti86}

\title{U.S. consumption data}

\description{
   The \code{Blanciforti86} data frame contains annual U.S. consumption data
   from 1947 to 1978 or 1981.
   These data include eleven aggregate commodity groups:
   (1) food, (2) alcohol plus tobacco, (3) clothing, (4) housing,
   (5) utilities, (6) transportation, (7) medical care, (8) durable goods,
   (9) other nondurable goods, (10) other services, and
   (11) other miscellaneous goods;
   four food groups:
   (1) meats, (2) fruits and vegetables, (3) cereal and bakery products, and
   (4) miscellaneous foods;
   and four meat groups:
   (1) beef and veal, (2) pork, (3) fish, and (4) poultry.
   For each of these groups the consumption expenditures, price indices,
   and expenditure shares are available.

   This data frame contains the following variables / columns:
   \describe{
      \item{year}{The year.}
      \item{xAggX}{Expenditure on the aggregate commodity group X
         (in Millions of US-Dollars).}
      \item{xAgg}{Total expenditure on all eleven aggregate commodity groups
         (in Millions of US-Dollars).}
      \item{xcAggX}{Deflated expenditure on the aggregate commodity group X
         (in Million of 1972 US-Dollars).}
      \item{xcAgg}{Total deflated expenditure on all eleven aggregate commodity
         groups (in Million of 1972 US-Dollars).}
      \item{pAggX}{Price index for the aggregate commodity group X (1972 = 100).}
      \item{wAggX}{Expenditure share of the aggregate commodity group X.}
      \item{xFoodX}{Per capita expenditure for food group X (in US-Dollars).}
      \item{xFood}{Total per capita expenditure for all four food groups
         (in US-Dollars).}
      \item{xcFoodX}{Deflated per capita expenditure for food group X
         (in 1972 US-Dollars).}
      \item{xcFood}{Total deflated per capita expenditure for all four food
         groups (in 1972 US-Dollars).}
      \item{pFoodX}{Price index for food group X (1972 = 100).}
      \item{wFoodX}{Expenditure share of food group X.}
      \item{xFoodUsdoc}{Expenditure on food published by the US Department of
         Commerce (USDOC) (in Millions of US-Dollars).}
      \item{xFoodUsda}{Expenditure on food published by the US Department of
         Agriculture (USDA) (in Millions of US-Dollars).}
      \item{xFoodNew}{Expenditure on food (new estimates)
         (in Millions of US-Dollars).}
      \item{xMeatX}{Per capita expenditure for meat group X (in US-Dollars).}
      \item{xMeat}{Total per capita expenditure for all four meat groups
         (in US-Dollars).}
      \item{pMeatX}{Price index for meat group X (1972 = 100).}
      \item{cMeatX}{Per capita consumption of meat group X (in pounds).}
      \item{pMeat}{Price index for the meat group (1972 = 100).}
      \item{pMeatNew}{Price index for the meat group (new estimate)
         (1972 = 100).}
      \item{population3}{Population as presented in Table 5.A.3 (in  Millions).}
      \item{population12}{Population as presented in Table 5.A.12 (in  Millions).}
   }
}

\usage{data(Blanciforti86)}


\source{
   Blanciforti, Laura. A., Richard D. Green and Gordon A. King (1986)
   U.S. Consumer Behavior Over the Postwar Period:
   An Almost Ideal Demand System Analysis.
   Giannini Foundation Monograph Number 40, August 1986.
}

\keyword{datasets}
