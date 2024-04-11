args <- commandArgs(TRUE)
efile <- args[1L]
num_verts <- as.integer(args[2L])
num_edges <- as.integer(args[3L])
cat('writing to', efile, '\n')

#vert_len <- sample(10:20, num_verts, r=TRUE)
#all_verts <- vapply(seq_len(num_verts),
#                    \(i) paste(sample(letters, vert_len[i]), collapse=''),
#                    character(1L))

all_verts <- as.character(seq_len(num_verts))
edge_table <- cbind(sample(all_verts, num_edges, r=T), sample(all_verts, num_edges, r=T),
                    round(runif(num_edges), 3))
#all_edges <- vapply(seq_len(num_edges),
#                    \(i) paste(c(sample(all_verts, 2L), as.character(round(runif(1),3))), collapse='\t'),
#                    character(1L))
if(file.exists(efile)) file.remove(efile)
write.table(edge_table, file=efile, quote=FALSE, sep='\t', eol='\n',
            row.names=FALSE, col.names=FALSE)
#writeLines(all_edges, efile)
