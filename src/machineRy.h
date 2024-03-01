#ifndef MACHINERY_MAIN
#define MACHINERY_MAIN

#include <stdlib.h>
#include <string.h>

/*
 * These are needed for outmem_graph.c
 */
#include <stdint.h>
#include <stdio.h>
#include <limits.h>

// mmap not supported on windows
#ifdef _WIN32
  #define mmap(addr, len, prot, flags, fildes, off) error("This function is not available on Windows.\n")
  #define munmap(addr, len) error("This function is not available on Windows.\n")
  #define stat(name, st) error("This function is not available on Windows.\n")
  #define struct stat
  #define MAP_SHARED
  #define PROT_READ
  #define PROT_WRITE
#else
  #define HAS_MMAN
  #include <sys/mman.h>
  #include <sys/stat.h>
#endif

/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

/*
 * R_ext/Random.h is needed for internal RNG routines
 */
#include <R_ext/Random.h>

// for math functions
#include <math.h>
#include <time.h>

// R_init_machineRy.c
void R_init_machineRy(DllInfo *info);

/**** NeuralNetwork ****/
SEXP R_initNNptr(SEXP NLAYERS, SEXP LAYERSIZES, SEXP ACTIVFUNCS, SEXP INPUT_SIZE, SEXP LR, SEXP LFXN);
SEXP R_PredictForInput(SEXP INPUTVECTOR, SEXP nnPtr);
SEXP R_TrainForInput(SEXP INPUTVECTOR, SEXP OUTPUTVECTOR, SEXP nnPtr);
SEXP R_UpdateWeights(SEXP LOSSVEC, SEXP nnPtr);

/**** Random Forest ****/
SEXP R_learn_tree_classif(SEXP DATA, SEXP NROWS, SEXP NCOLS, SEXP CLASSES,
                          SEXP NCLASSES, SEXP TO_CHECK, SEXP MAX_DEPTH, SEXP MIN_NODESIZE);
SEXP R_get_treeptr(SEXP VolatilePtr, SEXP INDICES, SEXP THRESHOLDS, SEXP GINIS);
SEXP R_rfpredict(SEXP RF_Obj, SEXP DATA, SEXP L, SEXP NENTRIES);
SEXP test_bfs_q2tree(SEXP INDICES, SEXP THRESHOLDS, SEXP GINIS, SEXP LENGTH);
SEXP test_tabulate(SEXP V, SEXP L);

/**** Label Propagation ****/
SEXP R_fastLP(SEXP NETWORK, SEXP MAX_ITER, SEXP NLEN, SEXP SQUEUE);
SEXP R_convertgraph(SEXP V1, SEXP V2, SEXP WEIGHT, SEXP NROW, SEXP NVERT);
SEXP R_fastcount(SEXP LVEC, SEXP MOD, SEXP OFFSET);

/**** out of memory stuff ****/
SEXP R_hashedgelist(SEXP FILENAME, SEXP NUM_EFILES, SEXP TABNAME, SEXP TEMPTABNAME, SEXP QFILE, SEXP OUTDIR,
                    SEXP SEPS, SEXP CTR, SEXP ITER, SEXP VERBOSE, SEXP IS_UNDIRECTED);
SEXP R_write_output_clusters(SEXP CLUSTERFILE, SEXP HASHEDFILES, SEXP NUM_FILES, SEXP OUTFILE, SEXP SEPS);

#endif
