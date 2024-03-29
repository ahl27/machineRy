\name{fastlabel_oom}
\alias{fastlabel_oom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Out of Memory Fast Label Propagation
}
\description{
Runs Fast Label Propagation using disk space for constant memory complexity.
}
\usage{
fastlabel_oom(edgelistfiles, outfile=tempfile(),
              mode=c("undirected", "directed"),
              add_self_loops=FALSE, iterations=1L,
              returnTable=FALSE, verbose=TRUE,
              sep='\t', linesep='\n',
              tempfiledir=tempdir(), cleanup_files=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgelistfiles}{
Character vector of files to be processed. Each entry should be a machine-interpretable path to an edgelist file. See Details for expected format.
}

\item{outfile}{
File to write final clusters to. Optional, defaults to a temporary file.
}

\item{mode}{
Character string specifying whether edges should be interpreted as undirected (default) or directed. Can be "undirected", "directed", or an unambiguous abbreviation.
}

\item{add_self_loops}{
Should self loops be added to each node? If \code{TRUE}, adds a self loop with
weight 0.5 to each node.
}

\item{iterations}{
Number of iterations to run fast label propagation algorithm for. Set to a value of 0 or less for infinite iterations.
}

\item{returnTable}{
Should result of clustering be returned as a file, or a \code{data.frame} object? If \code{FALSE}, returns a character vector corresponding to the path of \code{outfile}. If \code{TRUE}, parses \code{outfile} using \code{\link{read.table}} and returns the result. Not recommended for very large graphs.
}

\item{verbose}{
Should status messages (output, progress, etc.) be displayed while running?
}

\item{sep}{
Character that separates entries on a line in each file in \code{edgelistfiles}. Defaults to tab, as would be expected in a \code{.tsv} formatted file. Set to \code{','} for a \code{.csv} file.
}

\item{linesep}{
Character that separates lines in each file in \code{edgelistfiles}. Defaults to a newline, but can be set to any arbitrary character if you have specially formatted files.
}
\item{tempfiledir}{
Character vector corresponding to the location where temporary files used during execution should be stored. Defaults to R's \code{\link{tempdir}}.
}

\item{cleanup_files}{
Should intermediary files be deleted when the process completes? Note that \code{outfile} will only be deleted if \code{returnTable=TRUE} AND \code{cleanup_files=TRUE}.
}
}
\details{
Very large graphs require too much RAM for processing on some machines. In a graph containing billions of nodes and edges, loading the entire structure into RAM is rarely feasible. This implementation uses disk space for storing representations of each graph. While this is slower than computing on RAM, it allows this algorithm to scale to graphs of enormous size while only taking a small amount of RAM. If your graph is small enough to fit into RAM, consider using the \code{LP_igraph()} function instead.

This function expects a set of edge list files. Each entry in the file is expected to be in the following:

\code{VERTEX1<sep>VERTEX2<sep>WEIGHT<linesep>}

This line defines a single edge between vertices \code{VERTEX1} and \code{VERTEX2} with weight \code{WEIGHT}. \code{VERTEX1} and \code{VERTEX2} are strings corresponding to vertex names, \code{WEIGHT} is a numeric value that can be interpreted as a \code{double}. The separators \code{<sep>} and \code{<linesep>} correspond to the arguments \code{sep} and \code{linesep}, respectively. The default arguments work for standard \code{.tsv} formatting, i.e., a file of three columns of tab-separated values.
}
\value{
If \code{returnTable=TRUE}, returns a \code{\link{data.frame}} object with two columns. The first column contains the name of each vertex, and the second column contains the cluster it was assigned to.

If \code{returnTable=FALSE}, returns a character vector of length 1. This vector contains the path to the file where clusters were written to. The file is formatted as a \code{.tsv}, with each line containing two tab separated columns (vertex name, assigned cluster)
}

\note{
Unweighted graphs are currently not directly supported. You can supply an unweighted graph by fixing all weights to a nonzero number. This is planned to be improved in the future.
}

\references{
Traag, V.A., Subelj, L. Large network community detection by fast label propagation. \emph{Sci Rep} \bold{13}, 2701 (2023). https://doi.org/10.1038/s41598-023-29610-z
}
\author{
Aidan Lakshman <AHL27@pitt.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{
While this algorithm can scale to very large graphs, it does have some internal limitations. First, nodes must be comprised of no more than 255 characters. If this limitation is restrictive, please feel free to contact me. Alternatively, you can increase the size yourself by changing the definition of \code{MAX_NODE_NAME_SIZE} in \code{src/outmem_graph.c}. This limitation is provided to decrease memory overhead and improve runtime, but arbitrary values are possible.

Second, nodes are indexed using 64-bit unsigned integers, with 0 reserved for other values. This means that the maximum possible number of nodes available is 2^64-2, which is about 18.5 quintillion.

Third, this algorithm uses disk space to store large objects. As such, please ensure you have sufficient disk space for the graph you intend to process. I've tried to put safeguards in the code itself, but funky stuff can happen when the OS runs out of space. A good conservative approximation of disk consumption is \eqn{10(3v+4e)} bytes, where \eqn{v} is the number of vertices and \eqn{e} the number of edges. If you only know the average node degree \eqn{d}, then it would be \eqn{10v(3+4d)} bytes. That means that a network with one billion nodes and average node degree 4 would take around 190GB of disk space. Note that this is an approximate upper bound--actual consumption is likely less, but may be close to this bound with lots of long vertex names.
}

\examples{
## See tests/oom_labelprop_test.R for a larger example
num_verts <- 20L
num_edges <- 20L
all_verts <- sample(letters, num_verts)
all_edges <- vapply(seq_len(num_edges),
      \(i) paste(c(sample(all_verts, 2L),
                   as.character(round(runif(1),3))),
                 collapse='\t'),
                    character(1L))
edgefile <- tempfile()
if(file.exists(edgefile)) file.remove(edgefile)
writeLines(all_edges, edgefile)
res <- fastlabel_oom(edgefile, returnTable=TRUE)
print(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
