#ifndef RANDOMFOREST_INC
#define RANDOMFOREST_INC

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// these structs will only be used in this file
// they should not be expected to persist

struct DTreeNode {
  struct DTreeNode *left;
  struct DTreeNode *right;
  double threshold;
  int index;
};
typedef struct DTreeNode DTN;

struct DTNqueue{
  struct DTNqueue *next;
  DTN *ptr;
};
typedef struct DTNqueue queue;

DTN *initNode();

void freeDecisionTree(DTN *tree);
void bfs_q2tree(int *indices, double *thresholds, int length);

/* testing stuff */
struct printQ {
  struct printQ *next;
  DTN *node;
  int index, depth;
  double threshold;
};
void printDecisionTree(DTN *tree);

#endif