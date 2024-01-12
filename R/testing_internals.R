# balanced trees
test_balanced_tree <- function(tree_depth){
  tree_depth <- as.integer(tree_depth)
  indices <- sample(1L:9L, 2**tree_depth-1, replace=TRUE)
  thresholds <- trunc(runif(2**tree_depth-1, min=0, max=10)*10) / 10
  gini_gain <- runif(2**tree_depth-1)

  indices[(2**(tree_depth-1)):(length(indices))] <- -1L

  print(indices)
  print(thresholds)
  .Call("test_bfs_q2tree", indices, thresholds, gini_gain, length(indices))
}

test_unbalanced_tree <- function(tree_depth, subtree_depth){
  tree_depth <- as.integer(tree_depth)
  subtree_depth <- vapply(as.integer(subtree_depth), \(x) min(x, tree_depth-2L), integer(1L))
  indices <- sample(1L:9L, 2**tree_depth-1, replace=TRUE)
  thresholds <- trunc(runif(2**tree_depth-1, min=1, max=10)*10) / 10
  gini_gain <- runif(2**tree_depth-1)

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
  gini_gain <- gini_gain[-to_remove]
  print(indices)
  print(thresholds)
  .Call("test_bfs_q2tree", indices, thresholds, gini_gain, length(indices))
}

test_grow_tree <- function(nvars, nclasses, nsamp,
                           max_depth, to_samp=as.integer(floor(sqrt(nvars))),
                           min_nodesize=1L){
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
  r <- .Call("R_learn_tree_classif",
        d, nrow(d), ncol(d),
        cls, nclasses, to_samp,
        max_depth, min_nodesize)
  names(r) <- c("indices", "thresholds", "gini_gain")

  exptr <- .Call("R_get_treeptr", NULL, r[[1]], r[[2]], r[[3]], length(r[[1]]))

  ## get the original values back with inverse.rle
  res <- structure(exptr, Indices=rle(r[[1]]), Thresholds=rle(r[[2]]),
                   Gini=rle(r[[3]]), Size=length(r[[1]]),
                   nvars=nvars, nclasses=nclasses)
  return(res)
}

test_run_rf <- function(nentries){
  df <- data.frame(v1=runif(nentries), v2=runif(nentries),
                   v3=rnorm(nentries), v4=sample(c(T,F), nentries, r=T),
                   y=as.factor(sample(1:3, nentries, r=TRUE)))

  rf <- RandForest(y~., data=df, ntree=100L, nodesize=1)
  subsamp <- sample(seq_len(nentries), 100L)
  res <- predict(rf, df[subsamp,])
  res <- cbind(res, as.integer(df$y[subsamp]))

  pred_corr <- vapply(seq_len(nrow(res)), \(i){
    vals <- which(res[i,1:3] == max(res[i,1:3]))
    res[i,4]%in%vals
  }, logical(1L))
  cat("We got ", sum(pred_corr), "% correct\n", sep='')

  require('randomForest')
  pd <- randomForest(y~., data=df, ntree=100L, nodesize=1L)
  pred <- predict(pd, df[subsamp,])
  res <- cbind(pred, df$y[subsamp])
  pred_corr <- vapply(seq_len(nrow(res)), \(i){
    res[i,1] == res[i,2]
  }, logical(1L))
  cat("randomForest got ", sum(pred_corr), "% correct\n", sep='')
  return(NULL)
  #return(list(model=rf, results=res, data=df))
}
