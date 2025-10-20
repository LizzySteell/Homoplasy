### Adapt RHI.char function to accommodate ordered characters - using copilot's help

####### 2025-10-20

#V1.4

# Updated RHI.char that uses the package's existing cost.matrix()
# - Removed local cost.matrix definition and now calls the repository's cost.matrix()
# - Keeps ordered argument optional (NULL), and uses Sankoff parsimony with cost.matrix() per ordered character
# - Uses the package-local get.states() and uncert.to.poly() helpers

RHI.char <- function(data, tree, n = 100, ordered = NULL){

  if(!requireNamespace("phangorn", quietly = TRUE)) stop("phangorn package required.")
  if(!requireNamespace("ape", quietly = TRUE)) stop("ape package required.")

  # Ensure package-local helpers exist
  if(!exists("get.states", mode = "function")) stop("get.states() not found — make sure Functions/get.states is sourced or the package is loaded.")
  if(!exists("uncert.to.poly", mode = "function")) stop("uncert.to.poly() not found — make sure Functions/uncert.to.poly is sourced or the package is loaded.")
  if(!exists("cost.matrix", mode = "function")) stop("cost.matrix() not found — make sure Functions/cost.matrix is sourced or the package is loaded.")

  # Convert to a simple matrix for state extraction (polymorphisms in round brackets)
  data_matrix <- uncert.to.poly(as.matrix(data))

  # Number of characters
  nchar <- dim(as.matrix(data))[2]

  # Validate ordered argument and convert names -> indices if necessary
  if(!is.null(ordered)){
    if(is.character(ordered)){
      coln <- colnames(as.matrix(data))
      if(is.null(coln)) stop("data has no column names; supply ordered as numeric indices.")
      idx <- match(ordered, coln)
      if(any(is.na(idx))) stop("Some ordered character names not found in data column names.")
      ordered_idx <- sort(unique(idx))
    } else if(is.numeric(ordered)){
      if(any(ordered < 1 | ordered > nchar)) stop("ordered indices out of range.")
      ordered_idx <- sort(unique(as.integer(ordered)))
    } else {
      stop("ordered must be NULL, a numeric vector of indices, or a character vector of column names.")
    }
  } else {
    ordered_idx <- integer(0)
  }

  # Put each character (single-site phyDat) into a list
  char_list <- vector("list", nchar)
  for(i in seq_len(nchar)){
    char_list[[i]] <- data[, i, drop = FALSE]
  }

  # lmin for individual characters: total number states - 1
  lmin_chars <- integer(nchar)
  for(k in seq_len(nchar)){
    states_k <- get.states(data_matrix[, k])
    lmin_chars[k] <- length(states_k) - 1
  }

  # Subfunction to randomise tips and keep topology
  rtips <- function(tree){
    new_tips <- sample(tree$tip.label)
    Rtree <- tree
    Rtree$tip.label <- new_tips
    return(Rtree)
  }

  # Replicates tip randomising n times and creates a list of trees
  random_trees <- replicate(n, rtips(tree), simplify = FALSE)

  # Calculate Lnull (median) per character, accounting for ordered characters using Sankoff + cost.matrix()
  null_tlengths_char <- vector("list", nchar)
  for(m in seq_len(nchar)){
    states_m <- get.states(data_matrix[, m])
    nstates <- length(states_m)

    if(m %in% ordered_idx && nstates > 1){
      cost <- cost.matrix(states_m)
      null_tlengths_char[[m]] <- unlist(lapply(random_trees, function(rt){
        phangorn::parsimony(rt, char_list[[m]], method = "sankoff", cost = cost)
      }))
    } else {
      null_tlengths_char[[m]] <- unlist(lapply(random_trees, function(rt){
        phangorn::parsimony(rt, char_list[[m]])
      }))
    }
  }

  # Calculate median Lnull per character
  Lnull_char_med <- as.numeric(sapply(null_tlengths_char, function(x) stats::quantile(x, prob = 0.5)))

  # Calculate observed character length for each character (use Sankoff for ordered)
  char_lengths <- numeric(nchar)
  for(b in seq_len(nchar)){
    states_b <- get.states(data_matrix[, b])
    nstates_b <- length(states_b)

    if(b %in% ordered_idx && nstates_b > 1){
      cost_b <- cost.matrix(states_b)
      char_lengths[b] <- phangorn::parsimony(tree, char_list[[b]], method = "sankoff", cost = cost_b)
    } else {
      char_lengths[b] <- phangorn::parsimony(tree, char_list[[b]])
    }
  }

  # Subtract lmin
  char_lengths_minus_lmin <- char_lengths - lmin_chars
  Lnull_minus_lmin <- Lnull_char_med - lmin_chars

  # Compute RHI per character
  RHI_char <- char_lengths_minus_lmin / Lnull_minus_lmin

  # Handle divide-by-zero and NaN/Inf:
  zero_denom <- (Lnull_minus_lmin == 0)
  RHI_char[zero_denom & (char_lengths_minus_lmin == 0)] <- 0
  RHI_char[zero_denom & (char_lengths_minus_lmin != 0)] <- NA

  RHI_char[is.nan(RHI_char)] <- 0

  # Name the vector with column names if present, else with indices
  coln <- colnames(as.matrix(data))
  if(!is.null(coln)) names(RHI_char) <- coln else names(RHI_char) <- paste0("char", seq_len(nchar))

  return(RHI_char)
}
