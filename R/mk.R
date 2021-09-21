#' Equal rates Markov process model of discrete character evolution
#'
#' Creates a likelihood function for fitting a fully symmetric Markov process  
#' model of character evolution to categorical phenotype data.
#'
#' @param x A named vector of character state data having class \code{factor}.
#' @param phy An object of class \code{tree}.
#' @return A function having class \code{mk} that returns the likelihood of 
#' different input rates.
#' @examples
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' lik = mk(squamatareprod, phy)
#' optimize(lik, c(0, 1), maximum=TRUE)
mk = function(x, phy) {
    stopifnot(!is.null(names(x)))
    stopifnot(is.factor(x))
    
    x = x[tiplabels(phy)]

    g = matrix(0, Ntip(phy), nlevels(x))

    idx = (match(names(x), tiplabels(phy))-1) + (match(x, levels(x))-1)*Ntip(phy)

    g[idx+1] = 1
    g[which(!(tiplabels(phy) %in% names(x))), ] = 1

    model = .Call(C_mk_model, t(g), phy)

    structure(function(rate) {
        .Call(C_mk_loglikelihood, as.numeric(rate), model)
    }, class="mk")
}


#' Ancestral character state reconstruction
#'
#' Calculate marginal ancestral state probabilities under a fully symmetric 
#' Markov process model of character evolution.
#'
#' @param lik A function having class \code{mk}.
#' @return A function that returns marginal ancestral state probabilities for
#' different input rates.
#' @examples
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' lik = mk(squamatareprod, phy)
#' fit = optimize(lik, c(0, 1), maximum=TRUE)
#' mk.asr(lik)(fit$maximum)
mk.asr = function(lik) {
    stopifnot(inherits(lik, "mk"))
    model = environment(lik)$model
    function(rate) {
        t(.Call(C_mk_marginal_asr, as.numeric(rate), model))
    }
}
