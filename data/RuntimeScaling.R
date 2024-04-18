RuntimeResults <- rbind(
  c(0.20,0.18,0.23),
  c(0.18,0.19,0.26),
  c(0.19,0.18,0.26),
  c(0.20,0.20,0.27),
  c(0.19,0.19,0.32), # 1000
  c(0.22,0.19,0.37),
  c(0.21,0.20,0.56),
  c(0.24,0.21,0.89), # 10000
  c(0.58,0.38,4.32),
  c(1.36,0.59,9.21), # 100,000
  c(3.48,0.93,18.44),
  c(5.44,1.15,23.50),
  c(16.78,2.33,51.09), # 500,000
  c(55.05,5.22,101.21), # 1,000,000
  c(267.86,17.89,274.16),
  c(768.53,39.84,644.08),
  c(4030,86,1432) #10,000,000
)

RuntimeResults <- cbind(c(50,100,250,500,1000,2500,5000,10000,50000,
                          100000, 200000, 250000, 500000, 1000000,
                          2500000, 5000000, 10000000), RuntimeResults)

colnames(RuntimeResults) <- c("nvertedge", "igraph", 'lp in-mem', 'oom lp')

USE_LAST_ROWS <- 5L
to_use <- seq(nrow(RuntimeResults)-USE_LAST_ROWS, nrow(RuntimeResults))

if(interactive()){
  plot(NULL, xlim=c(1,max(RuntimeResults[,1])), ylim=c(1, max(RuntimeResults[,-1],na.rm=T)),
       log='yx', xlab='num vertices and edges', ylab='runtime (seconds)')
  for(i in seq_len(3)){
    lines(x=RuntimeResults[,1], y=RuntimeResults[,1+i],col=i)
    l <- lm(log(RuntimeResults[to_use,1+i]) ~ log(RuntimeResults[to_use,1]))
    cat("Scaling of", colnames(RuntimeResults)[1+i], "is", coefficients(l)[2], '\n')
  }
}
