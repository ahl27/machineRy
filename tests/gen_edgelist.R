num_verts <- 10
num_edges <- 200
vert_len <- sample(5:25, num_verts, r=TRUE)
vert_len <- rep(3L, num_verts)
all_verts <- vapply(seq_len(num_verts),
                    \(i) paste(sample(letters, vert_len[i]), collapse=''),
                    character(1L))
#all_verts <- sample(letters, num_verts)
outfile <- '/Users/aidan/Downloads/tmp_clustering_scratch/to_keep/edgelist.tsv'
all_edges <- vapply(seq_len(num_edges),
                    \(i) paste(c(sample(all_verts, 2L), as.character(round(runif(1),3))), collapse='\t'),
                    character(1L))
if(file.exists(outfile)) file.remove(outfile)
writeLines(all_edges, outfile)
