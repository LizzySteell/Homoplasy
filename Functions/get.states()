### 2024-02-08

### Function 'get.states()'

### To get charater states from a phylogenetic matrix. Ignores polymorphisms (round brackets only), gaps and '?'


get.states <- function(data) {
	
	#data - matrix array of phylogenetic dataset
		
		#Change gaps to ?
		data[data == "-"] <- "?"
		
	scores <- row.names(table(data)) #Makes a vector of all the possible scores, including ambiguities and polymorphisms
	 	
	 	#Make a subfunction to extract polymorphisms from parentheses - round brackets
	 	paren_extract <- function(x) {
   # Find characters within round brackets
     extracted <- gsub("[()]", "", x)
    return(extracted)}
		
		#Make vectors of polymorphisms - with and without round brackets
			 #This has no brackets
			polys <- unlist(lapply(scores, paren_extract))
			
			if(length(polys) != 0){
			
			 #This has brackets
			polys_brackets <- unlist(lapply(polys, parenthesize))
			
		states <- scores[!scores %in% c("?", polys_brackets)] #Makes a vector of the possible states - '0', '1', '2' etc. Excludes ambiguities and polymorphisms.
		
			} else { #If no polymorphisms/uncertainties
			
			states <- scores[!scores %in% c("?")] #Makes a vector of the possible states - '0', '1', '2' etc. Excludes ambiguities.
		}
		
		return(states) }
