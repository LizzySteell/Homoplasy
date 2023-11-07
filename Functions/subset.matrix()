# 2023-10-08

### Function for subsetting matrices based on a list of subset trees

subset.matrix <- function(trees, matrix) {
#trees - list of subset trees (object class phylo) from the 'subset.trees()' function
#matrix - character taxon matrix to be subsetted
	n <- length(trees)
	tips <- TipLabels(trees)
	datasets <- list()
		for (i in 1:n) {
		dataset <- matrix[c(unlist(tips[i])),]
		datasets[[i]] <- dataset
		}
	return(datasets)
	}
#Returns a list of subset matrices for each subsampled tree supplied.
