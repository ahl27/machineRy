#include "machineRy.h"
#include "fortran_includes.h"
#include "randomforest.h"

static DTN *rf_head;

DTN *initNode(){
  DTN *node = malloc(sizeof(DTN));
  node->threshold = 0;
  node->index = -1;
  node->gini_gain = 0;
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

SEXP R_learn_tree_classif(SEXP DATA, SEXP NROWS, SEXP NCOLS, SEXP CLASSES, SEXP NCLASSES, SEXP TO_CHECK, SEXP MAX_DEPTH){
  // array input
  double *data = REAL(DATA);
  int *class_response = INTEGER(CLASSES);

  // variable inputs
  int nrows = INTEGER(NROWS)[0];
  int ncols = INTEGER(NCOLS)[0];
  int nclasses = INTEGER(NCLASSES)[0];
  int num_to_check = INTEGER(TO_CHECK)[0];
  int max_depth = INTEGER(MAX_DEPTH)[0];

  // internal vars
  DTN *head = initNode();

  // helper function will destroy data and class_response, so duplicate them first
  double *dup_data = malloc(sizeof(double)*nrows*ncols);
  int *dup_class_response = malloc(sizeof(int)*nrows);

  // these do not need to be free'd -- will be free'd in the helper function
  Rprintf("Duplicating memory.\n");
  dup_data = memcpy(dup_data, data, sizeof(double)*nrows*ncols);
  dup_class_response = memcpy(dup_class_response, class_response, sizeof(int)*nrows);

  Rprintf("Training Tree:\n");
  learntreeclassif_helper(head, dup_data, dup_class_response, nrows, ncols, nclasses,
                          num_to_check, 0, max_depth);
  Rprintf("Done training.\n");
  // now we should have our entire tree created, and our duplicated arrays destroyed.

  // for testing, we'll print out the decision tree:
  Rprintf("Printing...\n");
  printDecisionTree(head);

  Rprintf("Freeing...\n");
  // not sure what to do with it yet; for now just free it and return NULL
  freeDecisionTree(head);

  Rprintf("Finished successfully!\n");
  return(R_NilValue);
}

void learntreeclassif_helper(DTN *node, double *data, int *class_response,
                              int nrows, int ncols, int nclasses, int num_to_check,
                              int cur_depth, int max_depth){
  // IMPORTANT: *data is assumed to malloc'd elsewhere
  // *data WILL BE FREE'd IN THIS FUNCTION
  // the same is true of *class_response
  // Calling function should duplicate any memory that cannot be destroyed

  // First we'll check to see if all the entries are the same class
  int curclass = nrows==0 ? -1 : class_response[0];
  int foundDifferent=0;
  for(int i=0; i<nrows; i++){
    if(curclass != class_response[i]){
      foundDifferent=1;
      break;
    }
  }
  if(!foundDifferent) cur_depth = max_depth;

  // if already at max_depth, just assign the most prominent class and return
  if(cur_depth == max_depth){
    // these should be preinitialized but it doesn't hurt to be safe
    node->index=-1;
    node->threshold = 1; // default value is the first class
    node->left = NULL;
    node->right = NULL;
    node->gini_gain = 0;

    if(!foundDifferent){
      // find the most prominent class
      int *classcounts = calloc(nclasses, sizeof(int));
      for(int i=0; i<nrows; i++)
        classcounts[class_response[i]-1]++;

      for(int i=1; i<nclasses; i++)
        if(classcounts[i] > classcounts[(int)(node->threshold-1)])
          node->threshold = (double)(i+1);

      free(classcounts);
    } else {
      // if we already know they're all the same, no need to check all again
      node->threshold = (double)curclass;
    }

    free(class_response);
    free(data);
    return;
  }

  // otherwise we need to split into nodes
  split_decision_node_classif(node, data, class_response,
                              nrows, ncols, nclasses, num_to_check);

  double splitpoint = node->threshold;
  int ind = node->index;
  int nrow_left = 0, nrow_right=0;
  double *v = &data[nrows*ind];

  // How big do we need the new arrays to be?
  for(int i=0; i<nrows; i++){
    if(v[i] <= splitpoint)
      nrow_left++;
    else
      nrow_right++;
  }

  double *left_data = malloc(sizeof(double) * nrow_left*ncols);
  double *right_data = malloc(sizeof(double) * nrow_right*ncols);
  int *left_class = malloc(sizeof(int) * nrow_left);
  int *right_class = malloc(sizeof(int) * nrow_right);
  int ctr_l=0, ctr_r=0;
  for(int i=0; i<nrows*ncols; i++){
    if(v[i%nrows] <= splitpoint){
      left_data[ctr_l] = data[i];
      if(ctr_l < nrow_left)
        left_class[ctr_l] = class_response[i%nrows];
      ctr_l++;
    } else {
      right_data[ctr_r] = data[i];
      if(ctr_r < nrow_right)
        right_class[ctr_r] = class_response[i%nrows];
      ctr_r++;
    }
  }

  // FREEING INPUT DATA/CLASS_RESPONSE
  free(data);
  free(class_response);

  DTN *left_node = initNode();
  DTN *right_node = initNode();

  // left node
  learntreeclassif_helper(left_node, left_data, left_class, nrow_left,
                          ncols, nclasses, num_to_check, cur_depth+1, max_depth);
  // right node
  learntreeclassif_helper(right_node, right_data, right_class, nrow_right,
                          ncols, nclasses, num_to_check, cur_depth+1, max_depth);

  node->left = left_node;
  node->right = right_node;

  return;
}

void split_decision_node_classif(DTN *node, double *data, int *class_response,
                                  int nrows, int ncols, int nclass, int num_to_check){
  // data should always be a numeric
  // response should be an int ranging from 1:n
  // nclass, num_to_check are constant throughout execution of the program

  // data will be a matrix stored by column (first nrows entries are col1, second are col2, etc.)
  // we'll just assume that all the preprocessing is done in R, no need to fiddle with that here
  // processing the SEXPs will be done separately so we can repeatedly call this internally

  // setting up a random sample of ints
  int *cols = malloc(sizeof(int) * ncols);
  for(int i=0; i<ncols; i++) cols[i] = i;
  int choice, tmp;

  // shuffle the columns
  GetRNGstate();
  for(int i=ncols-1; i>0; i--){
    choice = floor(unif_rand()*i);
    tmp = cols[choice];
    cols[choice] = cols[i];
    cols[i] = tmp;
  }
  PutRNGstate();

  Rprintf("\t\t");
  double *results = malloc(sizeof(double) * num_to_check);
  double *gini_gain = malloc(sizeof(double) * num_to_check);
  double curmax = -0.5;
  choice = -1;
  for(int i=0; i<num_to_check; i++){
    F77_CALL(find_gini_split)(&data[nrows*cols[i]], class_response, &nrows, &nclass, &results[i], &gini_gain[i]);
    if(gini_gain[i] > curmax){
      choice = i;
      curmax = gini_gain[i];
    }
  }
  Rprintf("\n");

  node->threshold = results[choice];
  node->index = cols[choice];
  node->gini_gain = curmax;


  free(results);
  free(gini_gain);
  free(cols);

  return;
}

void split_decision_node_regress(double *data, double *dbl_response, double *num_response, int nrows, int ncols){
  // not implemented yet
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

SEXP test_tabulate(SEXP V, SEXP L){
  int l = INTEGER(L)[0];
  double *outval = malloc(sizeof(double) * l);
  int *outcount = malloc(sizeof(int) * l);
  //int *ctr = malloc(sizeof(int));
  int ctr;
  F77_CALL(tabulate_double)(REAL(V), INTEGER(L), outval, outcount, &ctr);

  SEXP outv = PROTECT(allocVector(VECSXP, 2));
  SEXP VALS = PROTECT(allocVector(REALSXP, ctr));
  SEXP COUNTS = PROTECT(allocVector(INTSXP, ctr));

  memcpy(REAL(VALS), outval, sizeof(double) * ctr);
  memcpy(INTEGER(COUNTS), outcount, sizeof(int) * ctr);

  SET_VECTOR_ELT(outv, 0, VALS);
  SET_VECTOR_ELT(outv, 1, COUNTS);

  free(outval);
  free(outcount);

  UNPROTECT(3);
  return(outv);
}

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
