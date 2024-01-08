# balanced trees
test_balanced_tree <- function(tree_depth){
  tree_depth <- as.integer(tree_depth)
  indices <- sample(1L:9L, 2**tree_depth-1, replace=TRUE)
  thresholds <- trunc(runif(2**tree_depth-1, min=0, max=10)*10) / 10

  indices[(2**(tree_depth-1)):(length(indices))] <- -1L

  print(indices)
  print(thresholds)
  .Call("test_bfs_q2tree", indices, thresholds, length(indices))
}

test_unbalanced_tree <- function(tree_depth, subtree_depth){
  tree_depth <- as.integer(tree_depth)
  subtree_depth <- vapply(as.integer(subtree_depth), \(x) min(x, tree_depth-2L), integer(1L))
  indices <- sample(1L:9L, 2**tree_depth-1, replace=TRUE)
  thresholds <- trunc(runif(2**tree_depth-1, min=1, max=10)*10) / 10

  indices[(2**(tree_depth-1)):(length(indices))] <- -1L
  print(indices)
  l <- vector('list', length=length(subtree_depth))
  for(i in seq_along(subtree_depth)){
    v <- subtree_depth[i]
    trim_tree <- sample(1:(2**v), 1L)
    trim_tree <- 2**(v)-1 + trim_tree
    num_added <- 0
    while(num_added+v+1 < tree_depth){
      nv <- 2*trim_tree[(length(trim_tree)-(2**num_added)+1):length(trim_tree)]
      trim_tree <- c(trim_tree, nv, nv+1)
      num_added <- num_added+1
    }
    l[[i]] <- trim_tree
  }
  to_flip <- integer(length(l))
  to_remove <- integer(0L)
  for(i in seq_along(l)){
    to_flip[i] <- l[[i]][1L]
    to_remove <- c(to_remove, l[[i]][-1L])
  }
  to_flip <- unique(to_flip)
  to_remove <- unique(to_remove)
  print(to_flip)
  print(to_remove)
  indices[to_flip] <- -1L
  indices <- indices[-to_remove]
  thresholds <- thresholds[-to_remove]
  print(indices)
  print(thresholds)
  .Call("test_bfs_q2tree", indices, thresholds, length(indices))
}
