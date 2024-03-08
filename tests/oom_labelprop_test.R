library(machineRy)

basepath <- tempdir()

num_verts <- 10
num_edges <- 20
num_files <- 1L
vert_len <- sample(5:25, num_verts, r=TRUE)
all_verts <- vapply(seq_len(num_verts),
                    \(i) paste(sample(letters, vert_len[i]), collapse=''),
                    character(1L))

all_efiles <- paste0(file.path(basepath, 'edgelist'), seq_len(num_files), '.tsv')
all_efiles <- "~/Downloads/test_edgelist.tsv"
for(i in seq_len(num_files)){
  edgefile <- all_efiles[i]
  all_edges <- vapply(seq_len(num_edges),
                      \(i) paste(c(sample(all_verts, 2L), as.character(round(runif(1),3))), collapse='\t'),
                      character(1L))
  if(file.exists(edgefile)) file.remove(edgefile)
  writeLines(all_edges, edgefile)
}

library(igraph)
set.seed(123L)
res <- fastlabel_oom(all_efiles, outfile=tempfile(), iterations=1000L,
                     returnTable=TRUE, verbose=FALSE,
                     sep='\t', linesep='\n',
                     tempfiledir=tempdir(), cleanup_files=TRUE)

print(res)

graph_edge <- read.delim(all_efiles, header = FALSE)
set.seed(123L)
res2 <- LP_igraph(graph_from_data_frame(graph_edge, directed=FALSE), max_iterations=1000L)



comparison <- matrix(nrow=2, ncol=length(res2))
# get the same indexing that OOM clustering uses
colnames(comparison) <- unique(c(t(graph_edge[,1:2])))
rownames(comparison) <- c("inmem", "outmem")
comparison[1,names(res2)] <- res2
comparison[2,res[,1]] <- res[,2]
comparison[1,] <- match(comparison[1,], unique(comparison[1,]))
comparison[2,] <- match(comparison[2,], unique(comparison[2,]))
comparison
