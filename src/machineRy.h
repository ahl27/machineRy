#ifndef MACHINERY_MAIN
#define MACHINERY_MAIN

#include <stdlib.h>
#include <string.h>

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

#endif
