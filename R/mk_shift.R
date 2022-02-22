#' Heterogeneous rate models of discrete character evolution
#'
#' Fits multi-rate fully symmetric Markov process models of character 
#' evolution to categorical phenotype data.
#'
#' @param x A named vector of character state data having class \code{factor}.
#' @param phy An object of class \code{tree}.
#' @return A list with three components:
#' \describe{
#' \item{avg.rates}{Model-averaged branch rates for each branch in \code{phy}.
#' Rates are calculated by averaging over the set of evaluated rate-shift
#' configurations. The i-th rate in this vector corresponds to the branch 
#' subtending the node having index i.}
#' \item{aic.weights}{Akaike weights for each rate-shift configuration. The
#' rates in \code{avg.rates} are computed using these weights. Note that only
#' even-valued indices correspond to evaluated rate-shift configurations.}
#' \item{rate}{A function to return branch rates associated with a specific
#' rate-shift configuration. It takes a single argument -- the index of a valid
#' rate-shift configuration -- and returns a vector of estimated branch rates.}
#' }
#' @examples
#' \dontrun{
#' # For a more complete example see
#' demo(mk.shift, package="mk", ask=TRUE)
#' 
#' data(squamatatree)
#' data(squamatareprod)
#' phy = read.newick(text=squamatatree)
#' fit = mk.shift(squamatareprod, phy)
#' rate.bin = findInterval(fit$avg.rates, 
#'     seq(0, max(fit$avg.rates), length.out=33))
#' edge.color = colorRampPalette(
#'     c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))(33)[rate.bin]
#' plot(phy, edge.color=edge.color, lwd=0.5)
#' }
mk.shift = function(x, phy) {
    stopifnot(!is.null(names(x)))
    stopifnot(is.factor(x))
    
    if (length(setdiff(tiplabels(phy), names(x))))
        stop("Some terminal nodes are missing data")

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
