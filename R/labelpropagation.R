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
  nvert <- nrow(cl)
  l <- vector("list", nvert)

  for(i in seq_along(l)){
    cnts <- .Call("R_fastcount", cl, nvert, i) / nclust
    to_pull <- cnts!=0
    p <- (seq_len(nvert))[to_pull]
    cnts <- cnts[to_pull]
    l[[i]] <- list(p, cnts)
  }

  run_clabel_prop(l, nvert)
}

run_clabel_prop <- function(network, max_iter=Inf){
  n <- length(network)
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

fast_convert_igraph <- function(g, add_self_loop=FALSE){
  df <- as_data_frame(g, 'both')
  vertex_names <- df$vertices[[1]]
  vertex_map <- match(vertex_names, (vertex_names <- unique(vertex_names)))
  names(vertex_map) <- vertex_names
  e1 <- vertex_map[(df$edges[[1]])]
  e2 <- vertex_map[(df$edges[[2]])]
  w <- df$edges[[3]]
  prelim_g <- .Call("R_convertgraph",
                    e1-1L, e2-1L, w,
                    nrow(df$edges), nrow(df$vertices))
  for(i in seq_along(prelim_g)){
    prelim_g[[i]][[1]] <- prelim_g[[i]][[1]] + 1L
    prelim_g[[i]][[2]] <- prelim_g[[i]][[2]] / max(prelim_g[[i]][[2]])
    if(add_self_loop){
      prelim_g[[i]][[1]] <- c(prelim_g[[i]][[1]], i)
      prelim_g[[i]][[2]] <- c(prelim_g[[i]][[2]], 0.5)
    }
  }

  list(graph=prelim_g, vnames=vertex_names)
}

LP_igraph <- function(g, max_iterations, add_self_loop=FALSE, consensus=FALSE){
  elg <- fast_convert_igraph(g, add_self_loop)
  el <- elg$graph
  elvn <- elg$vnames
  useConsensus <- (is.logical(consensus) && consensus) || is.numeric(consensus)
  if(useConsensus){
    if(is.logical(consensus))
      to_iter <- c(0, 0.2,0.4,0.6,0.8,1,1.33,1.67,2)
    else
      to_iter <- consensus
    #to_iter <- c(seq(0, 0.9, by=0.1), seq(1,2,by=0.25))
    assignments <- matrix(nrow=length(el), ncol=length(to_iter))
    for(i in seq_along(to_iter)){
      l <- sigTransform_network(el, slope=to_iter[i], cutoff=0.2)
      assignments[,i] <- run_clabel_prop(l, max_iterations)
    }

    res <- consensus_clustering(assignments)
  } else {
    res <- run_clabel_prop(el, max_iterations)
  }
  names(res) <- elvn
  if(!any(grepl('[^0-9.]', elvn)))
    res <- res[order(as.numeric(elvn))]
  else
    res <- res[order(elvn)]
  res
}

deparse_communities <- function(commPred){
  comms <- rep(seq_along(commPred), times=lengths(commPred))
  names(comms) <- unlist(commPred)
  return(comms)
}

fastlabel_oom <- function(edgelistfiles, outfile=tempfile(),
                          mode=c("undirected", "directed"), add_self_loops=FALSE,
                          iterations=1L,
                          returnTable=FALSE, verbose=TRUE,
                          sep='\t', linesep='\n',
                          tempfiledir=tempdir(), cleanup_files=TRUE){
  if(!is.numeric(iterations)){
    stop("iterations must be a positive integer or numeric.")
  } else {
    iterations <- as.integer(iterations)
  }
  if(is.na(iterations) || is.null(iterations) || iterations < 1){
    warning("Invalid value of iterations, defaulting to 1.")
    iterations <- 1L
  }
  counter_cluster_binary <- tempfile(tmpdir=tempfiledir)
  csr_table_binary <- tempfile(tmpdir=tempfiledir)
  qfiles <- c(tempfile(tmpdir=tempfiledir),
              tempfile(tmpdir=tempfiledir),
              tempfile(tmpdir=tempfiledir))
  hashdir <- file.path(tempfiledir, "OOMhashes")
  mode <- match.arg(mode)
  is_undirected <- mode == "undirected"
  if(dir.exists(hashdir)){
    for(f in list.files(hashdir, full.names=TRUE))
      file.remove(f)
  } else {
    dir.create(hashdir)
  }

  if(verbose){
    cat("Temporary files stored at ", tempfiledir, "\n")
    cat("\tCSR: ", basename(csr_table_binary), "\n")
    cat("\tClusters: ", basename(counter_cluster_binary), "\n")
    cat("\tQueue 1: ", basename(qfiles[1]), "\n")
    cat("\tQueue 2: ", basename(qfiles[2]), "\n")
    cat("\tQueue counter: ", basename(qfiles[3]), "\n")
    cat("\tHashes: ", basename(hashdir), "\n")
  }

  seps <- paste(sep, linesep, sep='')
  ctr <- 1
  if(!all(file.exists(edgelistfiles))) stop("edgelist file does not exist")
  edgelistfiles <- normalizePath(edgelistfiles, mustWork=TRUE)
  # R_hashedgelist(tsv, csr, clusters, queues, hashdir, seps, 1, iter, verbose)
  .Call("R_hashedgelist", edgelistfiles, length(edgelistfiles), csr_table_binary,
        counter_cluster_binary, qfiles, hashdir, seps, ctr, iterations,
        verbose, is_undirected, add_self_loops)

  # R_write_output_clusters(clusters, hashes, length(hashes), out_tsvpath, seps)
  .Call("R_write_output_clusters", counter_cluster_binary, hashdir,
        outfile, "\t\n", verbose)

  if(cleanup_files){
    for(f in c(csr_table_binary, counter_cluster_binary, qfiles))
      if(file.exists(f)) file.remove(f)
    for(f in list.files(hashdir, full.names=TRUE))
      if(file.exists(f)) file.remove(f)
    file.remove(hashdir)
  }
  if(returnTable){
    tab <- read.table(outfile, sep="\t")
    colnames(tab) <- c("Vertex", "Cluster")
    if(file.exists(outfile)) file.remove(outfile)
    return(tab)
  } else {
    return(outfile)
  }
}
