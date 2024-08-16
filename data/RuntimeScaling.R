# Local results using `time`
# columns are igraph, in-memLP, OOMLP, MCL I2.0
# note that LP are both doing a single iteration
#   so runtimes slightly better than will be expected
RuntimeResults <- rbind(
  c(0.20,0.18,0.23,0.01),
  c(0.18,0.19,0.26,0.01),
  c(0.19,0.18,0.26,0.01),
  c(0.20,0.20,0.27,0.01),
  c(0.19,0.19,0.32,0.01), # 1000
  c(0.22,0.19,0.37,0.05),
  c(0.21,0.20,0.56,0.13),
  c(0.24,0.21,0.89,0.33), # 10000
  c(0.58,0.38,4.32,2.67),
  c(1.36,0.59,9.21,6.70), # 100,000
  c(5.44,1.15,23.50,26.8),
  c(16.78,2.33,51.09,86.93), # 500,000
  c(55.05,5.22,101.21,303.19), # 1,000,000
  c(267.86,17.89,274.16,1758.34),
  c(768.53,39.84,644.08,6685),
  c(4030,86,1432,26145) #10,000,000
)

RuntimeResults <- cbind(c(50,100,250,500,1000,2500,5000,10000,50000,
                          100000, 250000, 500000, 1000000,
                          2500000, 5000000, 10000000), RuntimeResults)

colnames(RuntimeResults) <- c("nvertedge", "igraph", 'lp in-mem', 'oom lp', "mcl_i20")

USE_LAST_ROWS <- 5L
to_use <- seq(nrow(RuntimeResults)-USE_LAST_ROWS, nrow(RuntimeResults))

if(interactive()){
  plot(NULL, xlim=c(10000,max(RuntimeResults[,1])), ylim=c(1, max(RuntimeResults[,-1],na.rm=T)),
      xlab='num vertices and edges', ylab='runtime (seconds)')
  for(i in seq_len(ncol(RuntimeResults)-1)){
    points(x=RuntimeResults[,1], y=RuntimeResults[,1+i],col=i, pch=19, cex=0.5)
    lines(x=RuntimeResults[,1], y=RuntimeResults[,1+i],col=i, type='l')
    df <- data.frame(nvert = RuntimeResults[to_use,1], time=RuntimeResults[to_use,1+i])
    l <- lm(log(time) ~ log(nvert), data=df)
    cat("Scaling of", colnames(RuntimeResults)[1+i], "is", coefficients(l)[2], '\n')
  }
}
