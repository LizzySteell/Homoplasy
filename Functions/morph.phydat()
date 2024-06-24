### 2023-11-08

### morph.phydat: PhyDat object wrapper function V1

# To make phyDat objects (from Phangorn package) for morphological datasets. Includes changing gaps to ambiguities, changing uncertainties to polymorphisms, and creating a contrast matrix.
# Dependencies: packages 'phangorn'

# Polymorphisms and uncertainties are treated as the same here.

# Uses the contrast.matrix() and uncert.to.poly() functions - see functions documentation.


morph.phydat <- function(matrix){
	
	# Matrix of class matrix array which a phylogenetic dataset; taxa as rows, characters as columns.
	# Polymorphisms or uncertainties indicated with round brackets only
	
	# Create contrast matrix
	contrast <- contrast.matrix(matrix)
	
	# Remove gaps
	matrix[matrix == "-"] <- "?"
	
	# Change uncertainties into polymorphisms, if any
	
	matrix <- uncert.to.poly(matrix)
	
	# Create phyDat object
	
	phydat <- phyDat(matrix, type="USER", contrast=contrast)
	
	return(phydat)
	
}
