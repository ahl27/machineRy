init_nn <- function(input_size, layersizes,
                    activfuncs=rep(1L, length(layersizes)),
                    learning_rate=0.05,
                    loss_fxn=1L){
  if(is.numeric(layersizes) && !is.integer(layersizes))
    layersizes <- as.integer(layersizes)
  # infov <- list(layers=c(input_size, layersizes),
  #               activations=c("ReLU", "Sigmoid", "Softplus", "Tanh", "Identity")[c(5L,activfuncs+1L)],
  #               learning_rate=learning_rate,
  #               loss=c("AbsError", "SqError")[loss_fxn+1L])
  .Call("R_initNNptr",
        length(layersizes), layersizes, activfuncs,
        input_size, learning_rate, loss_fxn)
}

train_nn <- function(formula, data, layers, learning_rate=0.05,
                     epoch_size=1L, activation_functions=NULL, verbose=TRUE){
  if(missing(data))
    data <- environment(formula)
  if(is.numeric(epoch_size)){
    if(!is.integer(epoch_size))
      epoch_size <- as.integer(epoch_size)
    if(length(epoch_size) > 1){
      warning("epoch_size has length greater than 1; discarding values after the first...")
      epoch_size <- epoch_size[1]
    }
  } else {
    stop("epoch_size must be an integer!")
  }
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", 'data'), names(mf), 0L)
  mf <- mf[m]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, 'terms')
  y <- model.response(mf, 'numeric')
  x <- model.matrix(mt, mf)
  nrx <- nrow(x)

  input_size <- ncol(x)
  if(!is.matrix(y)) y <- matrix(y, nrow=1L)
  output_size <- nrow(y)

  nn <- init_nn(input_size, c(layers,output_size), learning_rate=learning_rate)
  nepochs <- (nrx %/% epoch_size) + (nrx %% epoch_size != 0)
  distrib <- sample(rep(seq_len(nepochs), length=nrx))
  progbar <- c('|', '/','-','\\')
  for(i in seq_len(nepochs)){
    if(verbose) cat("Epoch ", i, '/', nepochs, ', Loss: \\', sep='')
    ds <- x[distrib==i,,drop=FALSE]
    dsy <- y[,distrib==i,drop=FALSE]
    l <- ld <- matrix(NA_real_, nrow=nrow(ds), ncol=output_size)
    for(j in seq_len(nrow(ds))){
      if(verbose) cat('\b', progbar[i%%4 + 1L], sep='')
      lv <- .Call("R_TrainForInput",
                    ds[j,,drop=FALSE], dsy[,j,drop=FALSE],
                    nn, PACKAGE='machineRy')
      l[j,] <- lv[seq_len(output_size)]
      ld[j,] <- lv[seq_len(output_size)+output_size]
    }
    if(verbose) cat('\b', sprintf('%.02f', mean(l)), ' (', sprintf('%.02f', mean(ld)), ')\n', sep='')
    .Call("R_UpdateWeights", colMeans(ld), nn)
  }

  structure(nn, formula=formula, class='MLP_Model')
}
