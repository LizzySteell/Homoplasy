### 2024-02-08

### uncert.to.poly: A function to change uncertainties, denoted with curly brackets, to polymorphisms in a phylogenetic matrix

uncert.to.poly <- function(matrix) {
  # Define curly_extract function
  curly_extract <- function(x) {
    # Find characters within curly braces
    extracted <- gsub("[{}]", "", x)
    return(extracted)
  }
  
  scores <- rownames(table(matrix))
  
  # Extract uncertainties
  curly <- unlist(lapply(scores, curly_extract))
  
  if (length(curly) != 0) { # If uncertainties present
    for (i in curly) {
      round_value <- paste0("(", i, ")")
      matrix[matrix == paste0("{", i, "}")] <- round_value
    }
  }
  
  return(matrix)
}
