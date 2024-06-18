### 2024-02-08

### poly.to.amb(): A function to change polymorphisms to ambiguities '?' in a matrix

poly.to.amb <- function(matrix){
	
	# matrix: Phylogenetic character matrix of class matrix array
	
	# First change uncertainties (curly brackets) to polymorphisms
	
	matrix <- uncert.to.poly(matrix)
	
#Make a subfunction to extract polymorphisms from parentheses - round brackets
	 	paren_extract <- function(x) {
    # Find characters within round brackets
     extracted <- gsub("[()]", "", x)
    return(extracted)}
		
		scores <- rownames(table(matrix))
		
		# Extract polies
		
		poly <- unlist(lapply(scores, paren_extract))

			
			if(length(poly) != 0){ #If polies present:
				
				## Puts the polies into a vector with brackets
				round_bracket <- parenthesize(poly)
				
				#Change polies to ambiguities
				for(i in 1:length(round_bracket)){
					matrix[matrix == round_bracket[i]] <- "?"
				}
				
				
				
			}
			return(matrix)
				}
