#include "machineRy.h"
#include "fortran_includes.h"
#include "randomforest.h"

/*
 * Tree node initialization
 */
DTN *initNode(){
  DTN *node = malloc(sizeof(DTN));
  node->threshold = 0;
  node->index = -1;
  node->gini_gain = 0;
  node->left = NULL;
  node->right = NULL;
  return node;
}

/*
 * R-exposed functions
 */
SEXP R_learn_tree_classif(SEXP DATA, SEXP NROWS, SEXP NCOLS, SEXP CLASSES, SEXP NCLASSES, SEXP TO_CHECK, SEXP MAX_DEPTH, SEXP MIN_NODESIZE){
  // array input
  double *data = REAL(DATA);
  int *class_response = INTEGER(CLASSES);

  // variable inputs
  int nrows = INTEGER(NROWS)[0];
  int ncols = INTEGER(NCOLS)[0];
  int nclasses = INTEGER(NCLASSES)[0];
  int num_to_check = INTEGER(TO_CHECK)[0];
  int max_depth = INTEGER(MAX_DEPTH)[0];
  int min_nodesize = INTEGER(MIN_NODESIZE)[0];

  // internal vars
  DTN *head = initNode();

  // helper function will destroy data and class_response, so duplicate them first
  double *dup_data = malloc(sizeof(double)*nrows*ncols);
  int *dup_class_response = malloc(sizeof(int)*nrows);

  // these do not need to be free'd -- will be free'd in the helper function
  //Rprintf("Duplicating memory.\n");
  dup_data = memcpy(dup_data, data, sizeof(double)*nrows*ncols);
  dup_class_response = memcpy(dup_class_response, class_response, sizeof(int)*nrows);

  //Rprintf("Training Tree:\n");
  learntreeclassif_helper(head, dup_data, dup_class_response, nrows, ncols, nclasses,
                          num_to_check, 0, max_depth, min_nodesize);
  //Rprintf("Done training.\n");
  // now we should have our entire tree created, and our duplicated arrays destroyed.

  // for testing, we'll print out the decision tree:
  //Rprintf("Printing...\n");
  //printDecisionTree(head);

  // Now let's export it to an R object
  // not going to store the gini gain for now.

  //Rprintf("Exporting...\n");
  // these objects will be allocated in `export_internal_tree`
  int *indices = NULL;
  double *thresholds = NULL, *gini_gain=NULL;
  int l = 0;
  export_internal_tree(head, &indices, &thresholds, &gini_gain, &l);

  // This is one option, I'm instead just going to register the external
  // pointer right away and return it, since I think that's easier.
  // Avoids a double call, and most people will predict right after
  // training anyway.

  // Read values back into R
  SEXP R_retval = PROTECT(allocVector(VECSXP, 4));
  SEXP R_indices = PROTECT(allocVector(INTSXP, l));
  SEXP R_thresholds = PROTECT(allocVector(REALSXP, l));
  SEXP R_gini = PROTECT(allocVector(REALSXP, l));

  memcpy(INTEGER(R_indices), indices, sizeof(int)*l);
  memcpy(REAL(R_thresholds), thresholds, sizeof(double)*l);
  memcpy(REAL(R_gini), gini_gain, sizeof(double)*l);
  free(indices);
  free(thresholds);

  SET_VECTOR_ELT(R_retval, 1, R_indices);
  SET_VECTOR_ELT(R_retval, 2, R_thresholds);
  SET_VECTOR_ELT(R_retval, 3, R_gini);
  UNPROTECT(3);

  // register the external pointer and then return
  SEXP R_ptr = PROTECT(R_MakeExternalPtr(head, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(R_ptr, (R_CFinalizer_t) R_TreeFinalizer, TRUE);
  SET_VECTOR_ELT(R_retval, 0, R_ptr);
  UNPROTECT(2);

  return(R_retval);
}

SEXP R_get_treeptr(SEXP VolatilePtr, SEXP INDICES, SEXP THRESHOLDS, SEXP GINIS){
  // We're going to lazily evaluate these trees
  // if they exist, just return the external pointer
  // note that it seems NULL can be treated as an external pointer address for whatever reason
  if(VolatilePtr != R_NilValue && R_ExternalPtrAddr(VolatilePtr)) return(VolatilePtr);
  int madePtr = 0;

  // otherwise, create the tree and return an external pointer to it
  DTN *tree = bfs_q2tree(INTEGER(INDICES), REAL(THRESHOLDS), REAL(GINIS), LENGTH(INDICES));
  if(VolatilePtr == R_NilValue){
    VolatilePtr = PROTECT(R_MakeExternalPtr(tree, R_NilValue, R_NilValue));
    madePtr = 1;
  } else {
    R_SetExternalPtrAddr(VolatilePtr, tree);
  }
  R_RegisterCFinalizerEx(VolatilePtr, (R_CFinalizer_t) R_TreeFinalizer, TRUE);
  if(madePtr) UNPROTECT(1);
  return VolatilePtr;
}

SEXP R_rfpredict(SEXP RF_Obj, SEXP DATA, SEXP L, SEXP NENTRIES){
  // assume we transpose the input before it's passed in so that each column is one entry
  // pointer needs to be reprotected because it may not be protected after R_get_treeptr
  SEXP R_ptr = PROTECT(R_get_treeptr(VECTOR_ELT(RF_Obj, 0), VECTOR_ELT(RF_Obj, 1), VECTOR_ELT(RF_Obj, 2), VECTOR_ELT(RF_Obj, 3)));
  SET_VECTOR_ELT(RF_Obj, 0, R_ptr);
  DTN *tree = (DTN *) R_ExternalPtrAddr(R_ptr);

  int nentries = INTEGER(NENTRIES)[0];
  int len = INTEGER(L)[0];
  double *data = REAL(DATA);
  SEXP R_return = PROTECT(allocVector(REALSXP, nentries));
  double *retval = REAL(R_return);

  for(int i=0; i<nentries; i++)
    retval[i] = predict_for_input(tree, &data[i*len]);

  UNPROTECT(2);
  return(R_return);
}


/*
 * Import/Export trees between C and R
 */
DTN *bfs_q2tree(int *indices, double *thresholds, double *gini, int length){
  queue *q = malloc(sizeof(queue));
  queue *end = q;
  queue *tmp_q = q;
  DTN *tmp, *head;

  head = initNode();
  q->ptr = head;
  q->next = NULL;
  int i=0, cur_ind;

  while(q && i<length){
    // load value into queue
    cur_ind = indices[i];
    tmp = q->ptr;
    tmp->threshold = thresholds[i];
    tmp->gini_gain = gini[i];
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
  return head;
}

void export_internal_tree(DTN *tree, int **indices, double **thresholds, double **gini_gain, int *outlength){
  // Notable here: `indices` and `thresholds` are NOT ALLOCATED
  // They should be allocated in this function and then returned
  int outlen = 0;
  queue *q = malloc(sizeof(queue));
  queue *end = q;
  queue *tmp_q = q;

  q->ptr = tree;
  q->next = NULL;

  // create a BFS queue of the tree
  while(tmp_q){
    if(tmp_q->ptr->index != -1){
      end->next = malloc(sizeof(queue));
      end = end->next;
      end->ptr = tmp_q->ptr->left;
      end->next = malloc(sizeof(queue));
      end = end->next;
      end->ptr = tmp_q->ptr->right;
      end->next = NULL;
    }
    tmp_q = tmp_q->next;
    outlen++;
  }

  // allocate storage to be returned
  *indices = malloc(sizeof(int)*outlen);
  *thresholds = malloc(sizeof(double)*outlen);
  *gini_gain = malloc(sizeof(double)*outlen);

  // write values into storage (will be freed later)
  *outlength = outlen;
  for(int i=0; i<outlen; i++){
    (*indices)[i] = q->ptr->index;
    (*thresholds)[i] = q->ptr->threshold;
    (*gini_gain)[i] = q->ptr->gini_gain;
    tmp_q = q;
    q = q->next;
    free(tmp_q);
  }
  return;
}

/*
 * Other Helper Functions
 */
double predict_for_input(DTN *tree, double *data){
  DTN *tmp=tree;

  while(tmp->index != -1){
    if(data[tmp->index] <= tmp->threshold)
      tmp = tmp->left;
    else
      tmp = tmp->right;
  }

  return(tmp->threshold);
}

void learntreeclassif_helper(DTN *node, double *data, int *class_response,
                              int nrows, int ncols, int nclasses, int num_to_check,
                              int cur_depth, int max_depth, int min_nodesize){
  R_CheckUserInterrupt();
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

  // if already at max_depth or we have fewer observations than minimum node size,
  // just assign the most prominent class and return
  if(cur_depth == max_depth || nrows <= min_nodesize){
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
  if(ind == -1){
    node->left = NULL;
    node->right = NULL;
    free(data);
    free(class_response);
    return;
  }

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
                          ncols, nclasses, num_to_check, cur_depth+1,
                          max_depth, min_nodesize);
  // right node
  learntreeclassif_helper(right_node, right_data, right_class, nrow_right,
                          ncols, nclasses, num_to_check, cur_depth+1,
                          max_depth, min_nodesize);

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

  double *results = malloc(sizeof(double) * num_to_check);
  double *gini_gain = malloc(sizeof(double) * num_to_check);
  double curmax = 0;
  choice = -1;
  for(int i=0; i<num_to_check; i++){
    F77_CALL(find_gini_split)(&data[nrows*cols[i]], class_response, &nrows, &nclass, &results[i], &gini_gain[i]);
    if(gini_gain[i] >= curmax){
      // using geq in case we find a split that also leaves gini unchanged
      // a split with equal gini is better than just terminating (probably?)
      choice = i;
      curmax = gini_gain[i];
    }
  }
  if(choice == -1){
    node->index = -1;
    node->threshold = 1.0; //fix later
    node->gini_gain = 0.0;
  } else {
    node->threshold = results[choice];
    node->index = cols[choice];
    node->gini_gain = curmax;
  }

  free(results);
  free(gini_gain);
  free(cols);

  return;
}

void split_decision_node_regress(double *data, double *dbl_response, double *num_response, int nrows, int ncols){
  // not implemented yet
  return;
}

/*
 * Cleanup Functions
 */
void R_TreeFinalizer(SEXP TreePointer){
  if (!R_ExternalPtrAddr(TreePointer)) return;
  DTN *head = (DTN *) R_ExternalPtrAddr(TreePointer);
  freeDecisionTree(head);
  R_ClearExternalPtr(TreePointer);
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

SEXP test_bfs_q2tree(SEXP INDICES, SEXP THRESHOLDS, SEXP GINIS, SEXP LENGTH){
  int length = INTEGER(LENGTH)[0];
  int *indices = INTEGER(INDICES);
  double *thresholds = REAL(THRESHOLDS);
  double *ginis = REAL(GINIS);

  Rprintf("Reading vectors into tree...\n");
  DTN *head = bfs_q2tree(indices, thresholds, ginis, length);
  Rprintf("Done!\nPrinting Tree...\n");
  printDecisionTree(head);
  Rprintf("Done!\nFreeing Tree...\n");
  freeDecisionTree(head);
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
