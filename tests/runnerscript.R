library(machineRy)
num_vertices <- 100L
edgelistfile <- '/Users/aidan/Downloads/tmp_clustering_scratch/to_keep/edgelist.tsv'
outtable <- '/Users/aidan/Downloads/tmp_clustering_scratch/to_keep/counts.bin'
finaltable <- '/Users/aidan/Downloads/tmp_clustering_scratch/to_keep/csr.bin'

hashdir <- "/Users/aidan/Downloads/tmp_clustering_scratch/hashes"
seps <- "\t\n"
for(f in list.files(hashdir, full.names=TRUE)){
  file.remove(f)
}

if(file.exists(outtable)) file.remove(outtable)
if(file.exists(finaltable)) file.remove(finaltable)

all_edges <- read.table(edgelistfile)
incidence <- character(nrow(all_edges)*2)
for(i in seq_len(nrow(all_edges))){
  incidence[i*2-1] <- all_edges[[1]][i]
  incidence[i*2] <- all_edges[[2]][i]
}

# this is on purpose to get the ordering the same as in C
ordering <- unique(incidence)
for(i in seq_along(ordering)){
  all_conn <- c(all_edges[[1]][all_edges[[2]]==ordering[i]], all_edges[[2]][all_edges[[1]]==ordering[i]])
  res <- match(all_conn, ordering)-1
  cat(i-1, ":", res, '\n')
}
num_edges <- table(incidence)[unique(incidence)]

#num_edges <- sort(num_edges, decreasing=TRUE)
.Call("R_hashedgelist", edgelistfile, finaltable, outtable, hashdir, seps, 1)
#res <- .Call("test_outputs", finaltable, length(num_edges)+1L)
#res <- sort(res, decreasing=TRUE)

print(num_edges)
#print(res)
#print(res == num_edges)
#print(all(res==num_edges))
