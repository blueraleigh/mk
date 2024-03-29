% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mk.R
\name{mk.asr}
\alias{mk.asr}
\title{Ancestral character state reconstruction}
\usage{
mk.asr(lik)
}
\arguments{
\item{lik}{A function having class \code{mk}.}
}
\value{
A function that returns marginal ancestral state probabilities for
different input rates.
}
\description{
Calculate marginal ancestral state probabilities under a fully symmetric 
Markov process model of character evolution.
}
\examples{
data(squamatatree)
data(squamatareprod)
phy = read.newick(text=squamatatree)
lik = mk(squamatareprod, phy)
fit = optimize(lik, c(0, 1), maximum=TRUE)
mk.asr(lik)(fit$maximum)
}
