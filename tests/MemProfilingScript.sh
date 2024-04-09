#!/bin/sh

# This script requires either gtime (gnu-time, OSX) or time (Linux)
# Set the below alias accordingly
alias ftime="gtime -f 'TIME:\t%E total\t%Us user\t%Ss sys\nMEMORY:\t%M(B) max'"

# Create the graph
R -f tests/MemProfilingSetup.R --args ${1} >/dev/null

# test igraph's performance
echo "Loading igraph:"
ftime R -e "suppressPackageStartupMessages(library(igraph))" >/dev/null
echo "\nRunning igraph::cluster_label_prop"
ftime R -e "suppressPackageStartupMessages(library(igraph));\
            g <- graph_from_data_frame(read.delim('TestTable.tsv', header = FALSE), directed=FALSE);\
            communities(cluster_label_prop(g))" >/dev/null

# test in-memory implementation
echo "-----\nLoading machineRy and igraph:"
ftime R -e "suppressPackageStartupMessages(library(machineRy));\
            suppressPackageStartupMessages(library(igraph))" >/dev/null
echo "\nRunning Fast LP in-memory:"
ftime R -e "suppressPackageStartupMessages(library(machineRy));\
            suppressPackageStartupMessages(library(igraph));\
            LP_igraph(graph_from_data_frame(read.delim('${1}', header = FALSE), \
                                                      directed=FALSE), \
                      add_self_loop = FALSE, max_iterations=1L)" >/dev/null

# test out of memory implementation
echo "-----\nLoading machineRy:"
ftime R -e "suppressPackageStartupMessages(library(machineRy))" >/dev/null
echo "\nRunning LP out of memory:"
ftime R -e "suppressPackageStartupMessages(library(machineRy));\
            fastlabel_oom('${1}', \
                          add_self_loops=FALSE,\
                          iterations=1L,\
                          returnTable=FALSE,\
                          verbose=FALSE)" >/dev/null
