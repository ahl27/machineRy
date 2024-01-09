#ifndef RANDOMFOREST_INC
#define RANDOMFOREST_INC

#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// these structs will only be used in this file
// they should not be expected to persist

struct DTreeNode {
  struct DTreeNode *left;
  struct DTreeNode *right;
  double threshold;
  double gini_gain;
  int index;
};
typedef struct DTreeNode DTN;

struct DTNqueue{
  struct DTNqueue *next;
  DTN *ptr;
};
typedef struct DTNqueue queue;

DTN *initNode();

// internal functions
void freeDecisionTree(DTN *tree);
void bfs_q2tree(int *indices, double *thresholds, int length);


void learntreeclassif_helper(DTN *node, double *data, int *class_response,
                              int nrows, int ncols, int nclass, int num_to_check,
                              int cur_depth, int max_depth);
void split_decision_node_classif(DTN *node, double *data, int *class_response,
                                  int nrows, int ncols, int nclass, int num_to_check);

/* testing stuff */
struct printQ {
  struct printQ *next;
  DTN *node;
  int index, depth;
  double threshold;
};
void printDecisionTree(DTN *tree);

#endif