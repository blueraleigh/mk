### Fitting multi-rate fully symmetric Markov process models of character 
### to empirical data using mk::mk.shift

library(mk)

# Data on squamate reproductive modes from
# Pyron, R.A. and F.T. Burbrink. 2014. Ecology Letters 17(1): 13-21
data(squamatatree)
data(squamatareprod)

phy = phylo::read.newick(text=squamatatree)

# A named vector of class ‘factor’ with reproductive mode
# classifications for 3962 squamates.
head(squamatareprod)

# As a first exploration of the data we'll simply plot the distribution
# of parsimony-inferred changes
mpr = parsimony::mpr.fitch(phy, squamatareprod)

# The minimum number of changes required to explain the distribution of
# reproductive modes observed at the tips of the tree
mpr$score

# The total number of MPR histories of character state change
mpr$mpr_count

# A function to draw MPR histories of character state change at random
mpr$simulate

# We'll just do one. This returns a matrix. Each column is a single
# history and rows represent nodes in the phylogeny. The cell value is
# the character state assigned to the node whose index is the same as the
# row index. For 'tree' class phylogenies used in the mk package, node indices
# are identical to the node indices of an ape 'phylo' object in 'cladewise'
# node order
h = mpr$simulate(1)

# We'll use the results stored in h to color branches according to their
# inferred reproductive mode. Oviparous (egg-laying) lineages will receive
# a blue color and viviparous ("live bearing") lineages will receive a
# red color.
edge.color = ifelse(h[,1] == 1L, 4L, 2L)


# Let's also use h to identify all instances of character state change
# To do so, we'll traverse the phylogeny and identify where a node's 
# inferred character state differs from that of its ancestor
state.change = logical(phylo::Nnode(phy))
for (node in phylo::descendants(root(phy), phy))
{
    if (h[node,] != h[phylo::parent(phy, node),])
        state.change[node] = TRUE
}

# Confirm that this equals the number of MPR inferred state changes from before
sum(state.change) == mpr$score

# Note that the plot function invisibly returns the plotting coordinates,
# which we'll store for later annotations
L = plot(phy, edge.color=edge.color)
title("Maximum parsimony inferred character states")

# Now highlight the nodes that underwent a state change
points(L[[1]][state.change, 1], L[[1]][state.change, 3], pch=21, 
    bg=edge.color[state.change], cex=0.8)


# We can observe that state changes are definitely clustered in certain
# regions of the phylogeny and that changes away from oviparity are much
# more common. 

# The conspicuous phylogenetic clustering of changes is certainly
# suggestive of rate heterogeneity, so let's fit a multi-rate model
# of reproductive mode evolution and see how it aligns with the
# distribution of inferred character state changes
fit = mk::mk.shift(squamatareprod, phy)

# This returns a list with three components.
#
# The first component is a vector of lineage-specific rates of character
# evolution that have been average over multiple possible rate-shift
# configurations. The second component is a set of weights assigned to
# different rate-shift configurations. These weights were used to compute
# the average rates in the first component. The final component is a function
# that can be used to retrieve the lineage-specific rates of character
# evolution associated with a particular rate-shift configuration.
str(fit)


# Let's first look at the model-averaged rates. Note that each value in
# this vector is the rate assigned to the edge leading to the node whose index 
# is the same as the vector index.
edge.rate = fit$avg.rates

# Let's discretize the rates and then map them to a set of colors, where
# warmer colors mean faster rates
rate.bin = findInterval(edge.rate, seq(0, max(edge.rate), length.out=33))
edge.color = colorRampPalette(
     c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))(33)[rate.bin]

plot(phy, edge.color=edge.color, lwd=0.5)
title("Model-averaged rates")
# Now highlight the nodes that underwent a state change
points(L[[1]][state.change, 1], L[[1]][state.change, 3], pch=21, bg=8, cex=0.8)

# And what you should see is that all the regions of phylogeny that experienced 
# changes are inferred to have elevated rates of evolution compared to regions
# that underwent little to no change.


# We can repeat the exact same process using individual rate-shift configurations
# as well, not just the model-averaged rates. Let's look at the individual
# configuration that has the highest score (highest AIC weight)
best_index = which.max(fit$aic.weights)

# Note that we have to be a little careful here due to different indexing
# conventions between R and C. Rather than pass best_index directly, we need
# to pass best_index - 1. This is why the names associated with the
# fit$aic.weights vector are always one less than the vector position.
edge.rate = fit$rate(best_index-1L)

rate.bin = findInterval(edge.rate, seq(0, max(edge.rate), length.out=33))
edge.color = colorRampPalette(
     c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))(33)[rate.bin]

# The overall pattern here is quite similar to the model-averaged pattern
plot(phy, edge.color=edge.color, lwd=0.5)
points(L[[1]][state.change, 1], L[[1]][state.change, 3], pch=21, bg=8, cex=0.8)
title("Best rate-shift configuration rates")
