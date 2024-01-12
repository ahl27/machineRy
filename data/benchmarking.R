library(machineRy)
library(randomForest)

nentries <- 1000L
train_bench <- function(){
to_iter <- seq(1,5,0.5)
vOrig <- vFort <- vFort2 <- numeric(length(to_iter))
for(i in seq_along(to_iter)){
  nentries <- floor(10**to_iter[i])
  cat("i = ", nentries, '\n')
  df <- data.frame(v1=runif(nentries), v2=runif(nentries),
                   v3=rnorm(nentries), v4=rgamma(nentries, 2L),
                   y=as.factor(sample(1:3, nentries, r=TRUE)))
  cat('\t - Original...\n')
  vOrig[i] <- system.time(randomForest(y~., data=df, nodesize=5, ntree=10L))[3]
  cat('\t - New...\n')
  vFort2[i] <- system.time(RandForest(y~., data=df, nodesize=5, verbose=FALSE))[3]
}
print(rbind(vOrig,vFort2))
}

## Prediction
prediction_bench <- function(){
to_iter <- seq(1,7,0.5)
vPredOrig <- vPredFortran <- numeric(length(to_iter))
nentries <- 1000L
df <- data.frame(v1=runif(nentries), v2=runif(nentries),
                 v3=rnorm(nentries), v4=rgamma(nentries, 2L),
                 y=as.factor(sample(1:3, nentries, r=TRUE)))
origModel <- randomForest(y~., data=df, nodesize=5, ntree=10L)
newModel <- RandForest(y~., data=df, nodesize=5, verbose=FALSE)

for(i in seq_along(to_iter)){
  nentries <- floor(10**to_iter[i])
  cat("i = ", nentries, '\n')
  df <- data.frame(v1=runif(nentries), v2=runif(nentries),
                   v3=rnorm(nentries), v4=rgamma(nentries, 2L),
                   y=as.factor(sample(1:3, nentries, r=TRUE)))

  cat('\t - Original...\n')
  vPredOrig[i] <- system.time(predict(origModel, df))[3]
  cat('\t - New...\n')
  vPredFortran[i] <- system.time(predict(newModel, df))[3]
}
print(rbind(vPredOrig, vPredFortran))
}
# x <- RandForest(y~., data=df[1:25,], nodesize=5, verbose=FALSE)
# x

memory_bench <- function(){
to_iter <- seq(1,7,0.5)
vOrigM <- vFortM <- integer(length(to_iter))
for(i in seq_along(to_iter)){
  nentries <- floor(10**to_iter[i])
  cat("i = ", nentries, '\n')
  df <- data.frame(v1=runif(nentries), v2=runif(nentries),
                   v3=rnorm(nentries), v4=rgamma(nentries, 2L),
                   y=as.factor(sample(1:3, nentries, r=TRUE)))
  cat('\t - Original...\n')
  vOrigM[i] <- object.size(randomForest(y~., data=df, nodesize=5, ntree=1L)$forest)[1]
  cat('\t - New...\n')
  vFortM[i] <- object.size(RandForest(y~., data=df, nodesize=5, ntree=1L, verbose=FALSE))[1]
}
print(rbind(vOrigM,vFortM))
}
