RuntimeResults <- rbind(
  c(0.21,0.19,0.31),
  c(0.21,0.20,0.31),
  c(0.20,0.20,0.27),
  c(0.20,0.19,0.34),
  c(0.22,0.19,0.34),
  c(0.20,0.20,0.50),
  c(0.23,0.21,0.75),
  c(0.26,0.23,1.34),
  c(0.63,0.49,7.19),
  c(1.14,0.79,17.11),
  c(3.66,1.37,43.00),
  c(4.25,1.67,57.91),
  c(20.92,3.78,165.09),
  c(NA,NA,NA)
)

RuntimeResults <- cbind(c(50,100,250,500,1000,2500,5000,10000,50000,
                          100000, 200000, 250000, 500000, 1000000), RuntimeResults)

colnames(RuntimeResults) <- c("nvertedge", "igraph", 'lp in-mem', 'oom lp')


plot(NULL, xlim=c(1,max(RuntimeResults[,1])), ylim=c(1, max(RuntimeResults[,-1],na.rm=T)),
     log='yx', xlab='num vertices and edges', ylab='runtime (seconds)')
for(i in seq_len(3)){
  lines(x=RuntimeResults[,1], y=RuntimeResults[,1+i],col=i)
  l <- lm(log(RuntimeResults[,1+i]) ~ log(RuntimeResults[,1]))
  cat("Scaling of", colnames(RuntimeResults)[1+i], "is", coefficients(l)[2], '\n')
}
