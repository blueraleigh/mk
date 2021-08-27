#' @example
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' lik = mk2(squamatareprod, phy)
mk2 = function(x, phy) {
    stopifnot(!is.null(names(x)))
    stopifnot(is.factor(x))
    stopifnot(nlevels(x) == 2)
    
    x = x[tiplabels(phy)]

    g = matrix(0, Ntip(phy), 2L)

    idx = (match(names(x), tiplabels(phy))-1) + (match(x, levels(x))-1)*Ntip(phy)

    g[idx+1] = 1
    g[which(!(tiplabels(phy) %in% names(x))), ] = 1

    model = .Call(C_mk2_model, t(g), phy)

    structure(function(rate) {
        if (length(rate) < 2)
            rate = c(rate, 1)
        .Call(C_mk2_loglikelihood, as.numeric(rate[1]), as.numeric(rate[2]), model)
    }, class="mk2")
}


#' @example
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' lik = mk2(squamatareprod, phy)
#' mk2.asr(lik)(0.0015)
mk2.asr = function(lik) {
    stopifnot(inherits(lik, "mk2"))
    model = environment(lik)$model
    function(rate) {
        if (length(rate) < 2)
            rate = c(rate, 1)
        t(.Call(C_mk2_marginal_asr, as.numeric(rate[1]), as.numeric(rate[2]), model))
    }
}
