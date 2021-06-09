#' @example
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' x = mptp.mk(squamatareprod, phy)
#' r8t = findInterval(x$avg.rates, seq(0, max(x$avg.rates), length.out=33))
#' edge.color = colorRampPalette(
#'    c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))(33)[r8t]
#' plot(phy, edge.color=edge.color)
mk.shift = function(x, phy) {
    stopifnot(!is.null(names(x)))
    stopifnot(is.factor(x))
    x = x[tiplabels(phy)]
    obj = .Call(C_mk_shift, x, phy)

    list(
        avg.rates = obj[[1L]]
        , aic.weights = structure(obj[[2L]], names=0:(Nnode(phy)-1))
        , rate = function(index) {
            .Call(C_mk_shift_backtrack, as.integer(index), phy)
        }
    )
}
