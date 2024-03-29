% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mk.R
\name{mk}
\alias{mk}
\title{Equal rates Markov process model of discrete character evolution}
\usage{
mk(x, phy)
}
\arguments{
\item{x}{A named vector of character state data having class \code{factor}.}

\item{phy}{An object of class \code{tree}.}
}
\value{
A function having class \code{mk} that returns the likelihood of 
different input rates.
}
\description{
Creates a likelihood function for fitting a fully symmetric Markov process  
model of character evolution to categorical phenotype data.
}
\examples{
data(squamatatree)
data(squamatareprod)
phy = read.newick(text=squamatatree)
lik = mk(squamatareprod, phy)
optimize(lik, c(0, 1), maximum=TRUE)
}
