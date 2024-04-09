#!/bin/sh

# This script requires either gtime (gnu-time, OSX) or time (Linux)
# Set the below alias accordingly
alias ftime="gtime -f 'TIME:\t%E total\t%Us user\t%Ss sys\nMEMORY:\t%K(B) avg\t%M(B) max'"

R -f tests/MemProfilingSetup.R --args ${1} >/dev/null
echo "Loading machineRy and igraph:"
ftime R -e "suppressPackageStartupMessages(library(machineRy));\
            suppressPackageStartupMessages(library(igraph))" >/dev/null
echo "\nRunning Fast LP in-memory:"
ftime R -e "suppressPackageStartupMessages(library(machineRy));\
            suppressPackageStartupMessages(library(igraph));\
            LP_igraph(graph_from_data_frame(read.delim('${1}', header = FALSE), \
                                                      directed=FALSE), \
                      add_self_loop = TRUE, max_iterations=1L)" >/dev/null
echo "-----\nLoading machineRy:"
ftime R -e "suppressPackageStartupMessages(library(machineRy))" >/dev/null
echo "\nRunning LP out of memory:"
ftime R -e "suppressPackageStartupMessages(library(machineRy));\
            fastlabel_oom('${1}', \
                          add_self_loops=TRUE,\
                          iterations=1L,\
                          returnTable=FALSE,\
                          verbose=FALSE)" >/dev/null
