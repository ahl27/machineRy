library(machineRy)

basepath <- tempdir()

num_verts <- 100
num_edges <- 200
num_files <- 2L
vert_len <- sample(5:25, num_verts, r=TRUE)
all_verts <- vapply(seq_len(num_verts),
                    \(i) paste(sample(letters, vert_len[i]), collapse=''),
                    character(1L))

all_efiles <- paste0(file.path(basepath, 'edgelist'), seq_len(num_files), '.tsv')
for(i in seq_len(num_files)){
  edgefile <- all_efiles[i]
  all_edges <- vapply(seq_len(num_edges),
                      \(i) paste(c(sample(all_verts, 2L), as.character(round(runif(1),3))), collapse='\t'),
                      character(1L))
  if(file.exists(edgefile)) file.remove(edgefile)
  writeLines(all_edges, edgefile)
}

res <- fastlabel_oom(all_efiles, outfile=tempfile(), iterations=1L,
                     returnTable=TRUE, verbose=TRUE,
                     sep='\t', linesep='\n',
                     tempfiledir=tempdir(), cleanup_files=TRUE)
print(res)