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

test_grow_tree <- function(nvars, nclasses, nsamp, max_depth, to_samp=as.integer(floor(sqrt(nvars)))){
  # R CMD INSTALL . && Rscript -e "library('machineRy'); test_grow_tree(5L, 3L, 100L, max_depth=5L)"
  set.seed(123L)
  cat("Generating test dataset...\n")
  d <- matrix(NA_real_, nrow=nsamp, ncol=nvars)
  cls <- integer(nsamp)
  nclasses <- as.integer(nclasses)

  for(i in seq_len(nvars)){
    m = runif(1, min=0, max=5)
    s = runif(1, min=0, max=2)
    d[,i] <- abs(rnorm(nsamp, mean=m, sd=s))
  }

  eqs <- lapply(seq_len(nclasses), \(i){
    runif(nvars, min=0, max=5)
  })

  for(i in seq_len(nrow(d))){
    p <- vapply(eqs, \(x) sum(abs(x * d[i,])), numeric(1L))
    cls[i] <- sample(seq_len(nclasses), 1,  prob = p)
  }

  print(cls)
  print(head(d))
  cat("Growing tree...\n")
  .Call("R_learn_tree_classif",
        d, nrow(d), ncol(d),
        cls, nclasses, to_samp, max_depth)
}
