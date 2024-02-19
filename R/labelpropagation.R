# 3340 total

# optimization priorities:
# - run_label_prop
# - tapply() call in consensus_clustering

logistic_transform <- function(x, slope=1, scale=0.5){
  1 / (1+exp(-slope*(x-scale)))
}

sigTransform_network <- function(l, slope=1, scale=0.5, cutoff=0.1){
  for(i in seq_along(l)){
    v <- l[[i]][[2]]
    v <- v / max(v)
    v <- logistic_transform(v, slope, scale)
    v <- v / max(v)
    p <- which(v >= cutoff)
    l[[i]][[1]] <- l[[i]][[1]][p]
    l[[i]][[2]] <- v[p]
  }
  l
}

consensus_clustering <- function(cl){
  nclust <- ncol(cl)
  l <- vector("list", nrow(cl))
  for(i in seq_along(l)){
    # huge time sink, big memory sink
    r <- table(unlist(apply(cl, 2, \(x) which(x==x[i])))) # 630
    l[[i]][[1]] <- as.integer(names(r))
    names(r) <- NULL
    r <- r / nclust
    l[[i]][[2]] <- r
  }
  # functionally infinite iterations,
  # 1e9 guards against infinite loops
  run_label_prop(l, Inf) # 620
}

run_clabel_prop <- function(network, max_iter=Inf){
  n <- length(network)
  #startqueue <- seq_len(n)
  startqueue <- sample(n)
  if(is.infinite(max_iter) || is.na(max_iter) || is.null(max_iter))
    max_iter <- -1L
  if(!is.integer(max_iter))
    max_iter <- as.integer(max_iter)
  partition <- .Call("R_fastLP", network, max_iter, n, startqueue)
  partition <- match(partition, unique(partition))
  nn <- names(network)
  if(is.null(nn)) nn <- seq_len(n)
  names(partition) <- nn
  return(partition)
}

run_label_prop <- function(network, max_iter=Inf){
  ## Expects a list
  ## Each entry is length 2:
  ##    network[[x]][[1]]: vertices that vertex x is connected to
  ##    network[[x]][[2]]: edge weights for edges in network[[x]][[1]]
  n <- length(network)
  partition <- seq_len(n)
  queue <- sample(n)
  #queue <- rev(seq_len(n))
  outQueue <- logical(n)

  pos <- n
  end <- 1L
  iter <- 1L

  while(iter < max_iter){
    i <- queue[pos]

    ## This is the section to optimize
    groups <- network[[i]][[1L]]
    if(length(groups) == 0){
      outQueue[i] <- TRUE
      pos <- pos - 1L
      if(pos == 0L){
        pos <- n
        iter <- iter + 1L
      }
      next
    }
    tabs <- tapply(network[[i]][[2L]], partition[groups], sum) # 100
    group <- names(tabs)[which.max(tabs)]

    curp <- as.integer(group)

    if(curp == partition[i]){
      outQueue[i] <- TRUE # remove from queue
    } else {
      # Update grouping
      partition[i] <- curp

      # Add all neighbors in different partitions to queue
      w <- groups[partition[groups] != curp]
      w <- w[outQueue[w]]
      outQueue[w] <- FALSE

      # Move to next element of queue
      end <- ((end - 2L) %% n) + 1L
      queue[end] <- i
      for (j in seq_along(w)) {
        end <- ((end - 2L) %% n) + 1L
        queue[end] <- w[j]
      }
    }
    if(i == 0){
      print(tabs)
      print(partition)
      print(w)
    }

    pos <- pos - 1L
    if(pos == 0L){
      pos <- n
      iter <- iter + 1L
    }

    if(pos == end) break
  }

  partition <- match(partition, unique(partition))
  nn <- names(network)
  if(is.null(nn)) nn <- seq_len(n)
  names(partition) <- nn
  return(partition)
}

fast_convert_igraph <- function(g, add_self_loop=FALSE){
  # TODO: handle named vertices
  df <- as_data_frame(g, 'both')
  prelim_g <- .Call("R_convertgraph",
                    as.integer(df$edges[[1]])-1L, as.integer(df$edges[[2]])-1L,
                    df$edges[[3]], nrow(df$edges), nrow(df$vertices))
  for(i in seq_along(prelim_g)){
    prelim_g[[i]][[1]] <- prelim_g[[i]][[1]] + 1L
    prelim_g[[i]][[2]] <- prelim_g[[i]][[2]] / max(prelim_g[[i]][[2]])
    if(add_self_loop){
      prelim_g[[i]][[1]] <- c(prelim_g[[i]][[1]], i)
      prelim_g[[i]][[2]] <- c(prelim_g[[i]][[2]], 0.5)
    }
  }

  prelim_g
}

LP_igraph <- function(g, max_iterations, add_self_loop=FALSE, consensus=FALSE, useC=FALSE){
  el <- fast_convert_igraph(g, add_self_loop)

  if(consensus){
    ncores <- max(parallel::detectCores()-1L, 1)
    ncores <- 1
    to_iter <- seq(0, 2, by=0.1)
    if(ncores > 1){
      a <- mclapply(to_iter, \(i){
        run_label_prop(sigTransform_network(el, slope=i, cutoff=0.2), max_iterations)
      }, mc.cores=ncores)
      assignments <- do.call(cbind, a)
    } else {
      assignments <- matrix(nrow=length(el), ncol=length(to_iter))
      for(i in seq_along(to_iter)){
        l <- sigTransform_network(el, slope=to_iter[i], cutoff=0.2) # 40
        assignments[,i] <- run_label_prop(l, max_iterations) # 1920
      }
    }
    res <- consensus_clustering(assignments) # 1230
  } else {
    if(useC){
      res <- run_clabel_prop(el, max_iterations)
    } else {
      res <- run_label_prop(el, max_iterations)
    }
  }
  res
}

deparse_communities <- function(commPred){
  n <- sum(lengths(commPred))
  res <- integer(n)
  names(res) <- seq_len(n)
  for(i in seq_along(commPred)){
    for(j in commPred[[i]]){
      res[j] <- i
    }
  }

  return(res)
}
