#include "machineRy.h"
#include "neuralnetwork.h"
#include "randomforest.h"


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define C_DEF(name, n)  {#name, (DL_FUNC) &name, n}

/*
 * -- REGISTRATION OF THE .Call ENTRY POINTS ---
 */
static const R_CallMethodDef callMethods[] = { // method name, num args
  CALLDEF(R_initNNptr, 6),
  CALLDEF(R_PredictForInput, 2),
  CALLDEF(R_TrainForInput, 3),
  CALLDEF(R_UpdateWeights, 2),
  CALLDEF(R_learn_tree, 9),
  CALLDEF(R_get_treeptr, 4),
  CALLDEF(R_rfpredict, 4),
  CALLDEF(R_fastLP, 4),
  CALLDEF(R_convertgraph, 5),
  CALLDEF(R_fastcount, 3),
  CALLDEF(R_hashedgelist, 15),
  CALLDEF(R_write_output_clusters, 5),
  CALLDEF(test_bfs_q2tree, 4),
  CALLDEF(test_tabulate, 2),
  {NULL, NULL, 0}
};

static const R_CMethodDef cMethods[] = {
  C_DEF(testnetwork, 0),
  //C_DEF(FreeNetwork, 1),
  {NULL, NULL, 0}
};

void R_init_machineRy(DllInfo *info){
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
