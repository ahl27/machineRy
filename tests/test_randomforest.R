library(machineRy)
#library(SynExtend)

## Regression
set.seed(199L)
n_samp <- 1000L
AA <- rnorm(n_samp, mean=1, sd=5)
BB <- rnorm(n_samp, mean=2, sd=3)
CC <- rgamma(n_samp, shape=1, rate=2)
err <- rnorm(n_samp, sd=0.5)
y <- AA + BB + 2*CC + err

d <- data.frame(AA,BB,CC,y)
train_i <- 1:500
test_i <- 501:1000
train_data <- d[train_i,]
test_data <- d[test_i,]

rf_regr <- RandForest(y~., data=train_data, rf.mode="regression", mtry=3, max_depth=5L)
as.dendrogram(rf_regr[[1]])
plot(rf_regr[[1]])

## Predict and visualize results
preds <- predict(rf_regr, test_data)

test_data2 <- cbind(DD=0, test_data)
preds2 <- predict(rf_regr, test_data2)

## this should always be true!
all(preds == preds2)

test_data3 <- test_data
test_data3$AA <- sample(c(TRUE,FALSE), nrow(test_data), replace=TRUE)
preds3 <- predict(rf_regr, test_data3)


limits <- c(min(c(preds[,1], test_data$y)), max(c(preds[,1], test_data$y)))
plot(x=test_data$y, y=preds[,1], xlim=limits, ylim=limits)
abline(a=0,b=1,col='blue')
sum((preds - test_data$y)**2)


## Classification
set.seed(234L)
v3 <- as.factor(sample(c("A", "B"), n_samp, replace=TRUE))
v1 <- vapply(seq_len(n_samp), \(i) rnorm(1, ifelse(v3[i]=='A', -1, 1)), numeric(1L))
v2 <- vapply(seq_len(n_samp), \(i) rnorm(1, ifelse(v3[i]=='A', 1, -1)), numeric(1L))
y <- ifelse(v3=='A', 1L, 3L)
y[y==1] <- y[y==1] + ifelse(v1[y==1] > -1, 0, 1)
y[y==3] <- y[y==3] + ifelse(v2[y==3] < -1, 0, 1)
y <- as.factor(y)
d <- data.frame(v1,v2,v3,y)
train_i <- 1:500
test_i <- 501:1000
train_data <- d[train_i,]
test_data <- d[test_i,]
rf_classif <- RandForest(y~., data=train_data, rf.mode="classification", mtry=3)


preds <- predict(rf_classif, test_data)

## build a confusion matrix
conf_mat <- matrix(0L, nrow=ncol(preds), ncol=ncol(preds))
for(i in seq_len(nrow(preds))){
	pos1 <- which.max(preds[i,])
	pos2 <- as.integer(test_data$y[i])
	conf_mat[pos1,pos2] <- conf_mat[pos1,pos2] + 1L
}
conf_mat