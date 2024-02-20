##### 2023-10-08

### SUBSET TREES FUNCTION

### Subsetting trees - randomly samples n tips from a given topology
subset.trees <- function(n, tree, ntaxa) {
	#n - number of subsets to be made
	#tree - tree (object of class phylo) to be subsetted
	#ntaxa - number of tips to retain in the subsetted trees

	trees <- replicate(n, drop.tip(tree,sample(tree$tip.label)[1:(length(tree$tip.label)-ntaxa)]), simplify=FALSE)
	return(trees)
	}
