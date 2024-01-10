
.safecheck_numeric <- function(v, argname, mustBePositive=FALSE){
  if(!is.numeric(v) || is.na(v) || is.null(v) || (mustBePositive && v < 0))
    stop("invalid value for ", argname)
  if(length(v) > 1){
    warning("discarding values for ", argname, " after the first entry")
    v <- v[1]
  }
  as.integer(v)
}

# These are going to change eventually
# for now, i'll just use placeholder functions
RandForest <- function(formula, data, subset, verbose=TRUE,
                   weights, na.action,
                   method='rf.fit', contrasts=NULL, ...){
  ## copying a lot of this from glm()
  if(missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if(identical(method, "model.frame"))
    return(mf)
  mt <- attr(mf, 'terms')
  y <- model.response(mf, "any")
  if(length(dim(y)) == 1L){
    nm <- rownames(y)
    dim(y) <- NULL
    if(!is.null(nm))
      names(y) <- nm
  }
  if(!is.empty.model(mt))
    x <- model.matrix(mt, mf, contrasts)
  else
    x <- matrix(NA_real_,nrow(y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  if(!is.numeric(x))
    stop("values supplied must be coercable to numeric")

  r <- RandForest.fit(x, y, weights=weights, verbose, ...)
  attr(r, 'formula') <- formula
  attr(r, 'contrasts') <- contrasts
  r
}

show.DecisionTree <- function(object){
  all_indices <- inverse.rle(object$indices)
  num_leaves <- sum(all_indices == -1)
  num_internal <- sum(all_indices != -1)
  cat(paste("A DecisionTree object with", num_leaves,
            "leaves and", num_internal, "internal nodes\n"))
}

show.RandForest <- function(object){
  l <- length(object)
  cat("A RandForest predictor containing", l, "trees:\n")
  if(l < 5){
    to_show <- seq_along(object)
  } else {
    to_show <- c(1,2,3,-1,l)
  }
  nc <- nchar(l)
  for(i in to_show){
    if(i==-1){
      ch <- paste0('\t', rep(' ', nc-1L), '.\n', collapse='')
      cat(rep(ch, 3))
    } else {
      firstchar <- paste0('\t', rep(' ', nc-nchar(i)), collapse='')
      cat(firstchar, i, '. ', sep='')
      show(object[[i]])
    }
  }
}

print.RandForest <- function(x, ...){
  show.RandForest(x)
}

print.DecisionTree <- function(x, ...){
  show.DecisionTree(x)
}

initDTStructure <- function(l){
  structure(list(pointer=l[[1]],
            indices=rle(l[[2]]),
            thresholds=rle(l[[3]]),
            ginis=rle(l[[4]])),
            class="DecisionTree")
}

RandForest.fit <- function(x, y=NULL, verbose=TRUE, ntree=10,
                           mtry=floor(sqrt(ncol(x))),
                           weights=NULL, replace=TRUE,
                           sampsize=if(replace) nrow(x) else ceiling(0.632*nrow(x)),
                           nodesize=1L, max_depth=NULL, ...){
  if(is.null(max_depth))
    max_depth <- -1L
  max_depth <- .safecheck_numeric(max_depth, 'max_depth', FALSE)
  ntree <- .safecheck_numeric(ntree, 'ntree')
  mtry <- .safecheck_numeric(mtry, 'mtry')
  nodesize <- .safecheck_numeric(nodesize, 'nodesize')
  sampsize <- .safecheck_numeric(sampsize, 'sampsize')

  classresponse <- as.integer(y)
  classnames <- levels(y)
  nclasses <- length(classnames)

  l <- vector('list', length(ntree))
  nr <- nrow(x)
  if(verbose){
    startt <- Sys.time()
    pb <- txtProgressBar(max=ntree, style=3)
  }
  for(i in seq_len(ntree)){
    subsamp <- sample(seq_len(nr), sampsize, replace=replace)
    r <- .Call("R_learn_tree_classif",
               x[subsamp,], length(subsamp), ncol(x),
               classresponse[subsamp], nclasses, mtry,
               max_depth, nodesize)
    l[[i]] <- initDTStructure(r)
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose){
    cat('\n  ')
    print(round(difftime(Sys.time(), startt), digits=2L))
  }

  attr(l, "class_levels") <- classnames
  attr(l, "num_vars") <- ncol(x)
  class(l) <- 'RandForest'
  l
}

predict.RandForest <- function(rf, newdata=NULL, na.action=na.pass){
  tt <- terms(attr(rf, 'formula'), data=newdata)
  noData <- (missing(newdata) || is.null(newdata))
  if(noData){
    x <- model.matrix(rf)
    mmDone <- TRUE
    return()
  } else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action)
    x <- model.matrix(Terms, m, contrasts.arg=attr(rf, 'contrasts'))
    mmDone <- FALSE
  }

  nentries <- nrow(x)
  nc <- ncol(x)
  results <- matrix(0.0, nrow=nentries, ncol=length(attr(rf, "class_levels")))
  colnames(results) <- attr(rf, "class_levels")
  for(i in seq_along(rf)){
    treeobj <- rf[[i]]
    for(i in seq(2,4))
      treeobj[[i]] <- inverse.rle(treeobj[[i]])
    p <- as.integer(.Call("R_rfpredict", treeobj, t(x), nc, nentries))
    idxs <- cbind(seq_len(nentries), p)
    results[idxs] <- results[idxs] + 1.0
  }

  results <- results / length(rf)
  results
}
