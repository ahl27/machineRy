#!/bin/sh

Rscript tests/MemProfilingSetup.R --args ${1}

valgrind --tool=massif `Rscript -e "library(machineRy);LP_igraph(graph_from_data_frame(read.delim(${1}, header = FALSE), directed=FALSE), add_self_loop = TRUE, max_iterations=1000L)"`
valgrind --tool=massif `Rscript -e "library(machineRy);fastlabel_oom(${1}, add_self_loops=TRUE, iterations=1000L,returnTable=TRUE)"`
valgrind --tool=massif `Rscript -e "library(machineRy)"`


valgrind --tool=massif `Rscript -e "library(machineRy);LP_igraph(graph_from_data_frame(read.delim(${1}, header = FALSE), directed=FALSE), add_self_loop = TRUE, max_iterations=1000L)"`
