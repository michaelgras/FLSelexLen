% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLselex.R
\name{selex.backtest}
\alias{selex.backtest}
\title{function to do nyears backtest of selex pattern in FLStocks}
\usage{
selex.backtest(
  sel,
  stock,
  sr = NULL,
  Fref = NULL,
  byears = 10,
  nyears = 3,
  plim = 0.975
)
}
\arguments{
\item{sel}{list of selex parameters of class of FLPars()}

\item{stock}{stock object of class FLStock}

\item{sr}{optional spawner-recruitment function FLSR}

\item{Fref}{option to input current F value, specify "catch" or  refpts = c("F0","Fmsy","F0.1","Fspr30","Fsq")}

\item{byears}{number of backtest years}

\item{nyears}{number of years for referencNULe conditions}

\item{plim}{set fbar for ages with Selectivy >= plim (default 0.975)}
}
\value{
FLStocks object
}
\description{
function to do nyears backtest of selex pattern in FLStocks
}
