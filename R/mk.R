#' @example
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' lik = mk(squamatareprod, phy)
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


#' @example
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' lik = mk(squamatareprod, phy)
#' mk.asr(lik)(0.0015)
mk.asr = function(lik) {
    stopifnot(inherits(lik, "mk"))
    model = environment(lik)$model
    function(rate) {
        t(.Call(C_mk_marginal_asr, as.numeric(rate), model))
    }
}
