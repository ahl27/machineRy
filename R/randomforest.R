## TODOs:
## 1. handle discrete/categorical data
## 2. handle regression
## 3. optimize performance


.safecheck_numeric <- function(v, argname, mustBePositive=FALSE){
  if(!is.numeric(v) || anyNA(v) || any(is.null(v)) || (mustBePositive && v < 0))
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
                   weights, na.action, method='rf.fit',
                   rf.mode=c('auto', 'classification', 'regression'), contrasts=NULL, ...){
  rf.mode <- match.arg(rf.mode)
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
  if(rf.mode == 'auto'){
    if(is.factor(y)){
      rf.mode <- "classification"
    } else if (is.numeric(y)){
      rf.mode <- "regression"
    } else {
      stop("response variable must be of type 'numeric' or 'factor'")
    }
  }

  isClassif <- rf.mode == 'classification'
  if(isClassif && !is.factor(y)){
    stop("response must be of type 'factor' for rf.mode='classification'")
  }
  if(!isClassif && !is.numeric(y)){
    stop("response must be of type 'numeric' for rf.mode='regression'")
  }

  ## TODO: check if this is robust -- what if we have multiple predictors?
  ## the response guaranteed to be the first outcome?
  all_data_types <- table(attr(mt, 'dataClasses'))
  if(rf.mode == "classification")
    all_data_types['factor'] <- all_data_types['factor'] - 1L
  all_data_types <- all_data_types[all_data_types!=0]

  if(length(dim(y)) == 1L){
    # I don't know what this is intended to catch
    nm <- rownames(y)
    dim(y) <- NULL
    if(!is.null(nm))
      names(y) <- nm
  }
  if(!isClassif){
    y <- unname(as.numeric(y))
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

  ## TODO: track names of variables -- `y~.` doesn't store names in `.`
  r <- RandForest.fit(x, y, method=rf.mode, weights=weights, verbose, ...)
  attr(r, 'formula') <- formula
  attr(r, 'contrasts') <- contrasts
  attr(r, 'mode') <- rf.mode
  r
}

show.DecisionTree <- function(object){
  #all_indices <- inverse.rle(object$indices)
  all_indices <- object$indices
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

initDTStructure <- function(l, predType, classnames){
  # structure(list(pointer=l[[1]],
  #           indices=rle(l[[2]]),
  #           thresholds=rle(l[[3]]),
  #           ginis=rle(l[[4]])),
  #           class="DecisionTree")
  structure(list(pointer=l[[1]],
            indices=l[[2]],
            thresholds=l[[3]],
            ginis=l[[4]],
            type=predType,
            class_levels=classnames),
            class="DecisionTree")
}

RandForest.fit <- function(x, y=NULL, verbose=TRUE, ntree=10,
                           mtry=floor(sqrt(ncol(x))),
                           weights=NULL, replace=TRUE,
                           sampsize=if(replace) nrow(x) else ceiling(0.632*nrow(x)),
                           nodesize=1L, max_depth=NULL,
                           method=NULL, ...){
  ## method will always be filled in by the main call
  useClassification <- method=="classification"
  if(is.null(max_depth))
    max_depth <- -1L
  max_depth <- .safecheck_numeric(max_depth, 'max_depth', FALSE)
  ntree <- .safecheck_numeric(ntree, 'ntree')
  mtry <- .safecheck_numeric(mtry, 'mtry')
  nodesize <- .safecheck_numeric(nodesize, 'nodesize')
  sampsize <- .safecheck_numeric(sampsize, 'sampsize')
  l <- vector('list', length(ntree))
  nr <- nrow(x)

  if(useClassification){
    response <- as.integer(y)
    classnames <- levels(y)
    nclasses <- length(classnames)
  } else {
    response <- as.numeric(y)
    classnames <- NULL
    if(length(unique(response)) <= 5 && nr > 5){
      warning("response variable has 5 or less unique values; are you sure you want to do regression?")
    }
    nclasses <- 0L
  }
  if(all(y==y[1])){
    stop("response variable has only one unique value!")
  }

  if(verbose){
    startt <- Sys.time()
    pb <- txtProgressBar(max=ntree, style=3)
  }
  for(i in seq_len(ntree)){
    subsamp <- sample(seq_len(nr), sampsize, replace=replace)
    r <- .Call("R_learn_tree",
               x[subsamp,], length(subsamp), ncol(x),
               response[subsamp], nclasses, mtry,
               max_depth, nodesize,
               useClassification)
    l[[i]] <- initDTStructure(r, method, classnames)
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose){
    cat('\n  ')
    print(round(difftime(Sys.time(), startt), digits=2L))
  }

  if(useClassification)
    attr(l, "class_levels") <- classnames
  attr(l, "num_vars") <- ncol(x)
  class(l) <- 'RandForest'
  l
}

predict.RandForest <- function(object, newdata=NULL, na.action=na.pass, ...){
  tt <- terms(attr(object, 'formula'), data=newdata)
  predmode <- attr(object, "mode")
  predmode <- match.arg(attr(object, 'mode'), c("classification", "regression"))
  useClassification <- predmode == 'classification'
  noData <- (missing(newdata) || is.null(newdata))
  if(noData){
    x <- model.matrix(object)
    mmDone <- TRUE
    return()
  } else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action)
    x <- model.matrix(Terms, m, contrasts.arg=attr(object, 'contrasts'))
    mmDone <- FALSE
  }

  nentries <- nrow(x)
  nc <- ncol(x)
  out_nc <- ifelse(useClassification, length(attr(object, "class_levels")), 1L)
  results <- matrix(0.0, nrow=nentries, ncol=out_nc)
  if(useClassification)
    colnames(results) <- attr(object, "class_levels")
  for(i in seq_along(object)){
    treeobj <- object[[i]]
    #for(i in seq(2,4))
    #  treeobj[[i]] <- inverse.rle(treeobj[[i]])
    p <- .Call("R_rfpredict", treeobj, t(x), nc, nentries)
    if(useClassification){
      p <- as.integer(p)
      idxs <- cbind(seq_len(nentries), p)
      results[idxs] <- results[idxs] + 1.0
    } else {
      results[,1] <- results[,1] + p
    }
  }

  results <- results / length(object)
  results
}

as.dendrogram.DecisionTree <- function(object){
  inds <- object$indices
  ths <- object$thresholds
  gains <- object$ginis
  isClassif <- object$type == "classification"
  levs <- attr(d, "class_levels")

  ## set up the root node
  if(length(inds) == 1L){
    d <- 1L
    attr(d, "leaf") <- TRUE
    attr(d, "label") <- ifelse(isClassif, levs[as.integer(ths[1])], sprintf("%0.2f", ths[1]))
  } else {
    d <- vector('list', length=2L)
  }
  attr(d, "ind") <- inds[1]
  attr(d, "thresh") <- ths[1]
  attr(d, "gain") <- gains[1]

  ## dend attributes to fix later
  attr(d, "members") <- 0L
  attr(d, "height") <- 0
  attr(d, "midpoint") <- 0
  class(d) <- "dendrogram"

  ## set up a queue
  cur_q <- list(1,2)
  ctr <- 2L
  leafctr <- 1L
  while(ctr <= length(inds)){
    acc <- cur_q[[1L]]
    i <- inds[ctr]
    th <- ths[ctr]
    g <- gains[ctr]
    cur_q <- cur_q[-1L]
    if(i != -1){
      ## append the next two entries to list
      node <- vector("list", 2L)
      cur_q <- c(cur_q, list(c(acc, 1L), c(acc,2)))
      #attr(node, "midpoint") <- 0
    } else {
      node <- leafctr
      leafctr <- leafctr + 1L
      attr(node, 'leaf') <- TRUE
      attr(node, "label") <- ifelse(isClassif, levs[as.integer(ths[ctr])], sprintf("%0.2f", ths[ctr]))
    }
    attr(node, "ind") <- inds[ctr]
    attr(node, "thresh") <- ths[ctr]
    attr(node, "gain") <- gains[ctr]
    attr(node, "members") <- 0L
    attr(node, "height") <- 0
    class(node) <- "dendrogram"

    d[[acc]] <- node
    ctr <- ctr + 1L
  }

  ## now we have to fix all the internal values
  .midpoint <- \(node) attr(node, 'midpoint')
  .height <- \(node) attr(node, 'height')
  all_v <- unlist(d)
  d <- dendrapply(d, \(y){
    if(!is.leaf(y)){
      is_leaf <- c(is.leaf(y[[1]]), is.leaf(y[[2]]))
      if(all(is_leaf)){
        mp <- 0.5
      } else if(is_leaf[1]){
        mp <- ((.midpoint(y[[2]]) + 1) / 2)
      } else if(is_leaf[2]){
        mp <- ((.midpoint(y[[1]]) + 1) / 2) + .midpoint(y[[1]])
      } else {
        mp <- (attr(y[[1]], 'members') + .midpoint(y[[1]]) + .midpoint(y[[2]])) / 2
      }
      attr(y, 'members') <- length(unlist(y))
      attr(y, 'height') <- max(.height(y[[1]]), .height(y[[2]])) + 1L
      attr(y, 'midpoint') <- mp
    } else {
      tmp <- as.integer(y)
      y[] <- match(tmp, all_v)
      attr(y, 'members') <- 1L
    }
    y
  }, how='post.order')
  return(d)
}

plot.DecisionTree <- function(x, y, text.adj=c(0,-0.05), text.cex=0.5, plotGain=FALSE, ...){
  ## TODO: figure out how to parse extended `text.` args into `text`
  x <- as.dendrogram(x)
  ## extract values from the tree
  ## builds a data.frame using postorder traversal
  r <- dendrapply(x, \(x){
    if(is.leaf(x)){
      data.frame(x=numeric(0L), y=numeric(0L), label=character(0L))
    } else {
      yv <- attr(x, "edgeplotpar")
      cur_v <- data.frame(x=attr(x, "xpos"), y=yv$h, label=yv$lab)
      rbind(x[[1]], x[[2]], cur_v)
    }
  }, how='post.order')

  ## plot dendrogram
  plot(x, ..., axes=FALSE)
  x <- dendrapply(x, \(dend){
    if(!is.leaf(dend)){
      xr <- unlist(dend)
      xr <- c(min(xr), max(xr))
      xv <- attr(dend, 'midpoint') + xr[1]
      yv <- attr(dend, 'height')
      vname <- paste0("V", attr(dend, 'ind'))
      thresh <- sprintf("%0.2f", attr(dend, 'thresh'))
      if(plotGain){
        gain <- sprintf("%+0.1f", attr(dend, 'gain'))
        text(xv+text.adj[1], yv+text.adj[2], bquote(.(vname) <= .(thresh) ~~ ( .(gain) )), cex=text.cex)
      } else {
        text(xv+text.adj[1], yv+text.adj[2], bquote(.(vname) <= .(thresh)), cex=text.cex)
      }
    }
    dend
  }, how='post.order')

  invisible(NULL)
}