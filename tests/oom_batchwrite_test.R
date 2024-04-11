library(machineRy)

# Verifying that we're writing out vertices correctly when using batch writes
NVERTS <- 100
VERTNAMELEN <- c(4,7)
VERTNAMES <- character(NVERTS)
while(length(unique(VERTNAMES)) < NVERTS){
  VERTNAMES <- unique(VERTNAMES)
  VERTNAMES <- c(vapply(seq_len(NVERTS-length(VERTNAMES)),
                        \(i) paste(sample(letters, sample(seq(VERTNAMELEN[1], VERTNAMELEN[2]), 1), r=T), collapse=''),
                        character(1L)),
                 VERTNAMES)
}

edgelist <- cbind(VERTNAMES, c(VERTNAMES[-1], VERTNAMES[1]), runif(NVERTS))

cat("Total degree should be:", sum(table(edgelist[,1])) + sum(table(edgelist[,2])), '\n')

edgefile <- tempfile()
ttn <- tempfile()
hashdir <- tempdir()
hashdir <- file.path(hashdir, "OOMhashes")
dir.create(hashdir)
write.table(edgelist, edgefile, sep='\t', eol='\n', quote=FALSE, col.names=FALSE, row.names=FALSE)

cat("Writing names to hash files...\n")
.Call("test_batchnamehash", edgefile, 1L, ttn, hashdir, "\t\n", 1, TRUE, TRUE)
cat("Checking what we wrote...\n")
outputfiles <- list.files(hashdir, full.names = TRUE)
.Call("test_batchnamehashverify", outputfiles, length(outputfiles), as.integer(NVERTS))
