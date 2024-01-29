### 2024-01-29

# Cost matrix function

### This function builds a standard cost matrix for ordered morphological characters whereby one character state change between two character states that are adjacent in value 'costs' one step, but a change between states that are sequentially higher or lower in value cost an increasing number of steps.

cost.matrix <- function(data)

# data - Phylogenetic matrix (in simple matrix array format) where row names are taxa and columns are characters. Gaps ('-') are treated as ambiguities ('?')

	{
		#Change all '-' values to '?'
	 	data[data == '-'] <- '?'
	 	
	 	#get character states
	 	
	 	states <- get.states(data)
	 	
	 	matrix <- matrix(nrow=length(states), ncol=length(states))
	 	
	 	dimnames(matrix) <- list(states, states)
	 	
	 	#fill in matrix
	 	
	 	for(i in 1:length(states)){
	 		
	 		for(j in 1:length(matrix[i,])){
	 			
	 			matrix[i,][j] <- j-i
	 		}
	 	}
	 	
	 	#Make values all positive
	 	
	 	matrix <- abs(matrix)
	 	
	}
