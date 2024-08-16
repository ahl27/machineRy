#!/bin/sh

# This script requires either gtime (gnu-time, OSX) or time (Linux)
# Set the below alias accordingly
alias ftime="gtime -f '%E %M'"

#for loopctr in 2500000 5000000 10000000
for loopctr in 10000000
do
#echo "$loopctr vertices and edges:"
# Create the graph
R -f tests/MemProfilingSetup.R --args ${1} $loopctr $loopctr >/dev/null

# test igraph's performance
#echo "igraph:"
ftime R -e "suppressPackageStartupMessages(library(igraph))" >/dev/null
# echo "Running igraph::cluster_label_prop"
# ftime R -e "suppressPackageStartupMessages(library(igraph));\
#             g <- graph_from_data_frame(read.delim('${1}', header = FALSE), directed=FALSE);\
#             communities(cluster_label_prop(g))" >/dev/null

#echo "Running igraph::cluster_label_prop"
ftime R -e "suppressPackageStartupMessages(library(igraph));\
            g <- graph_from_data_frame(read.delim('${1}', header = FALSE), directed=FALSE);\
            communities(cluster_louvain(g))" >/dev/null

# test in-memory implementation

# echo "Running Fast LP in-memory:"
# ftime R -e "suppressPackageStartupMessages(library(machineRy));\
#             suppressPackageStartupMessages(library(igraph));\
#             LP_igraph(graph_from_data_frame(read.delim('${1}', header = FALSE), \
#                                                       directed=FALSE), \
#                       add_self_loop = FALSE, max_iterations=1L)" >/dev/null

# # test out of memory implementation
# echo "Running LP out of memory:"
# ftime R -e "suppressPackageStartupMessages(library(machineRy));\
#             fastlabel_oom('${1}', \
#                           add_self_loops=FALSE,\
#                           iterations=1L,\
#                           returnTable=FALSE,\
#                           verbose=FALSE)" >/dev/null
done
