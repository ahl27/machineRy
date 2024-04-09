efile <- commandArgs(TRUE)[1L]
num_verts <- 1000
num_edges <- 2000

vert_len <- sample(5:15, num_verts, r=TRUE)
all_verts <- vapply(seq_len(num_verts),
                    \(i) paste(sample(letters, vert_len[i]), collapse=''),
                    character(1L))

all_edges <- vapply(seq_len(num_edges),
                    \(i) paste(c(sample(all_verts, 2L), as.character(round(runif(1),3))), collapse='\t'),
                    character(1L))
if(file.exists(efile)) file.remove(efile)
writeLines(all_edges, efile)
