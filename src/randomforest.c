#include "machineRy.h"
#include "randomforest.h"

static DTN *rf_head;

DTN *initNode(){
  DTN *node = malloc(sizeof(DTN));
  node->threshold = 0;
  node->index = -1;
  node->left = NULL;
  node->right = NULL;
  return(node);
}

void bfs_q2tree(int *indices, double *thresholds, int length){
  queue *q = malloc(sizeof(queue));
  queue *end = q;
  queue *tmp_q = q;
  DTN *tmp;

  rf_head = initNode();
  q->ptr = rf_head;
  q->next = NULL;
  int i=0, cur_ind;

  while(q && i<length){
    // load value into queue
    cur_ind = indices[i];
    tmp = q->ptr;
    tmp->threshold = thresholds[i];
    tmp->index = cur_ind;
    if(cur_ind > -1){
      // add both children of the node into the queue
      end->next = malloc(sizeof(queue));
      end = end->next;
      tmp->left = initNode();
      end->ptr = tmp->left;
      end->next = malloc(sizeof(queue));
      end=end->next;
      tmp->right = initNode();
      end->ptr = tmp->right;
      end->next = NULL;
    }

    i++;
    q = q->next;
  }

  // free the entire queue
  while(tmp_q){
    q = tmp_q;
    tmp_q = tmp_q->next;
    free(q);
  }
  return;
}

void freeDecisionTree(DTN *tree){
  if(!tree) return;
  freeDecisionTree(tree->left);
  freeDecisionTree(tree->right);
  free(tree);
}

/*
 * testing functions
 */
SEXP test_bfs_q2tree(SEXP INDICES, SEXP THRESHOLDS, SEXP LENGTH){
  int length = INTEGER(LENGTH)[0];
  int *indices = INTEGER(INDICES);
  double *thresholds = REAL(THRESHOLDS);

  Rprintf("Reading vectors into tree...\n");
  bfs_q2tree(indices, thresholds, length);
  Rprintf("Done!\nPrinting Tree...\n");
  printDecisionTree(rf_head);
  Rprintf("Done!\nFreeing Tree...\n");
  freeDecisionTree(rf_head);
  Rprintf("Done!\n");
  return(R_NilValue);
}

void printDecisionTree(DTN *tree){
  /*
   * print like this:
   *              v,ttt
   *         /              \
   *      v,ttt           v,ttt
   *     /     \         /     \
   *  v,ttt   v,ttt   v,ttt   v,ttt
   *   / \     / \     / \     / \
   * ttt ttt ttt ttt ttt ttt ttt ttt
   *
   */
  struct printQ *pq = malloc(sizeof(struct printQ));
  struct printQ *tmp = pq;
  struct printQ *end = pq;
  int maxDepth = 0, tmpDepth;
  pq->depth = 0;
  pq->node = tree;
  while(pq){
    pq->threshold = pq->node->threshold;
    pq->index = pq->node->index;
    if(pq->index != -1){
      tmpDepth = pq->depth + 1;
      end->next = malloc(sizeof(struct printQ));
      end = end->next;
      end->node = pq->node->left;
      end->depth = tmpDepth;
      end->next = malloc(sizeof(struct printQ));
      end = end->next;
      end->node = pq->node->right;
      end->depth = tmpDepth;
      end->next = NULL;

      if(tmpDepth > maxDepth)
        maxDepth = tmpDepth;
    }
    pq = pq->next;
  }
  Rprintf("\n");

  pq = tmp;
  // now we know how deep the tree is
  // each interior node will take 5 chars: "%d,%1.2f"
  // each leaf will take 3 chars: "%1.2f"
  // we need one space between each leaf, 3 between the next row, 5 between the next, etc.
  // total linewidth is thus the number of leaves*3 + (leaves-1) = 4*leaves - 1
  int curHeight = maxDepth;
  int padding = (((1<<curHeight)*4-1) - 5) / 2;
  int curidx=0;
  int *shouldprint = calloc(1<<(maxDepth-1), sizeof(int));
  shouldprint[0] = 1;
  shouldprint[1] = 1;
  // we can print while we free
  tmpDepth = 0;

  // this will be used to determine if we should skip printing subtrees
  end = pq;

  while(tmp){
    if(tmp->depth != tmpDepth){ // moved to a new level of the tree

      // Figure out which subtrees we should skip printing
      for(int i=(1<<tmpDepth); i>=0; i--){
        curidx = shouldprint[i];
        shouldprint[i*2] = curidx;
        shouldprint[i*2+1] = curidx;
      }
      curidx = 0;
      while(end->depth != tmp->depth){
        if(!shouldprint[curidx*2]){
          curidx++;
          continue;
        }
        shouldprint[curidx*2] *= end->index!=-1;
        shouldprint[curidx*2+1] *= end->index!=-1;
        curidx++;
        end = end->next;
      }
      curidx = 0;

      // new line
      Rprintf("\n");

      // figuring out new padding
      tmpDepth = tmp->depth;
      curHeight = maxDepth - tmpDepth;
      padding = (1<<curHeight)*4-1;
      padding = (padding - 5) / 2;
      padding = curHeight == 0 ? 0 : padding;

      // drawing branches
      Rprintf("%*s", padding + 2 + (curHeight>0), "");
      for(int i=0; i<(1<<tmpDepth); i++){
        if(shouldprint[i])
          Rprintf("%c", "/\\"[i%2]);
        else
          Rprintf(" ");
        if(curHeight == 0){
          Rprintf("%*s", i%2 ? 5 : 1, "");
        } else {
          Rprintf("%*s", 2*padding+3+4*((i)%2), "");
        }
      }
      Rprintf("\n");
    }
    while(!shouldprint[curidx] && curidx < (1<<tmpDepth)){
      Rprintf("%*s", padding, "");
      Rprintf("%*s", tmp->depth == maxDepth ? 3 : 5, "");
      Rprintf("%*s", padding+1, "");
      curidx++;
    }
    Rprintf("%*s", padding, "");
    if(tmp->index==-1){
      if(tmp->depth != maxDepth)
        Rprintf(" %1.1f ", tmp->threshold);
      else
        Rprintf("%1.1f", tmp->threshold);
    } else {
      Rprintf("%d,%1.1f", tmp->index, tmp->threshold);
    }
    Rprintf("%*s", padding, "");
    Rprintf(" ");
    tmp = tmp->next;
    curidx++;
  }
  Rprintf("\n");

  while(pq){
    tmp = pq;
    pq = pq->next;
    free(tmp);
  }
  free(shouldprint);
  return;
}
