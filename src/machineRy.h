#ifndef MACHINERY_MAIN
#define MACHINERY_MAIN

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

// for math functions
#include <math.h>

// R_init_machineRy.c
void R_init_machineRy(DllInfo *info);

/**** NeuralNetwork ****/
SEXP R_initNNptr(SEXP NLAYERS, SEXP LAYERSIZES, SEXP ACTIVFUNCS, SEXP INPUT_SIZE, SEXP LR, SEXP LFXN);
SEXP R_PredictForInput(SEXP INPUTVECTOR, SEXP nnPtr);
SEXP R_TrainForInput(SEXP INPUTVECTOR, SEXP OUTPUTVECTOR, SEXP nnPtr);
SEXP R_UpdateWeights(SEXP LOSSVEC, SEXP nnPtr);

static void FreeNetwork(SEXP nnPtr);

#endif
