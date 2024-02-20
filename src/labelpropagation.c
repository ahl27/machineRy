#include "machineRy.h"

struct llnode {
	struct llnode *next;
	int value;
};

struct queue {
	struct llnode *head;
	struct llnode *tail;
	int size;
};

typedef struct tabsort {
	int index;
	double weight;
} t_intdouble;

typedef struct graphnode {
	int *neighbors;
	double *weights;
} t_gnode;

typedef struct llnode ll;
typedef struct queue queue;

double roll_RRNG(){
	GetRNGstate();
	double rv = unif_rand();
	PutRNGstate();
	return rv;
}

void enqueue(queue *q, int ipt){
	ll *tmp = malloc(sizeof(ll));
	tmp->next = NULL;
	tmp->value = ipt-1; // subtract -1 since R indices are 1-indexed
	if(q->size){
		// if queue is nonempty, add element to tail position
		q->tail->next = tmp;
		q->tail = tmp;
		q->size++;
	} else {
		// if queue empty, initialize it
		q->head = tmp;
		q->tail = tmp;
		q->size = 1;
	}
	return;
}

int dequeue(queue *q){
	int v = -1;
	if(q->size){
		ll *tmp = q->head;
		v = tmp->value;
		q->head = tmp->next;
		q->size--;
		if(!(q->size))
			q->tail = NULL;
		free(tmp);
	}
	return v;
}

queue* alloc_queue(int *input, int size){
	// allocate the queue
	queue *q = malloc(sizeof(queue));
	q->head = NULL;
	q->tail = NULL;
	q->size = 0;

	for(int i=0; i<size; i++)
		enqueue(q, input[i]);

	return q;
}

void free_q(queue *q){
	while(q->size) dequeue(q);
	free(q);
}

int tid_cmpfunc(const void *a, const void *b){
	return (((t_intdouble *)a)->index - ((t_intdouble *) b)->index);
}

int tabulate(int *indices, double *weights, int *clusters, int l){
	t_intdouble sortingarray[l];

	for(int i=0; i<l; i++){
		sortingarray[i].index = clusters[indices[i]-1];
		sortingarray[i].weight = weights[i];
	}

  qsort(sortingarray, l, sizeof(t_intdouble), tid_cmpfunc);
  int max_group=-1, cur_group=sortingarray[0].index;
  double max_weight=-1, cur_weight=0;
  for(int i=0; i<l; i++){
  	cur_weight += sortingarray[i].weight;

  	if((i+1)==l || sortingarray[i+1].index != cur_group){
  		// end of group, check if we need to update
  		if(cur_weight > max_weight || (cur_weight == max_weight && roll_RRNG() < 0.5)){
  			max_weight = cur_weight;
  			max_group = cur_group;
  		}
  		cur_group = (i+1)==l ? -1 : sortingarray[i+1].index;
  		cur_weight = 0;
  	}
  }

  return max_group;
}

SEXP R_fastLP(SEXP NETWORK, SEXP MAX_ITER, SEXP NLEN, SEXP SQUEUE){
	int max_iterations = INTEGER(MAX_ITER)[0];
	int n = INTEGER(NLEN)[0];
	int *input_perm = INTEGER(SQUEUE);

	// dynamic allocation because worried about heap limits
	char *isNotInQueue = calloc(n, sizeof(char));
	int *clusters = malloc(n*sizeof(int));
	for(int i=0; i<n; i++) clusters[i] = i+1;

  // queue allocation
	queue *q = alloc_queue(input_perm, n);
	int iter = 0, subiter=0;
	int nextnode=-1;

	// getting node attributes
	SEXP R_node;
	int *connectedNodes, num_connected, new_group;
	double *weights;

	// main loop
	while(q->size && (iter != max_iterations)){
		// using neq means we can use max_iterations = -1 for infinite
		nextnode = dequeue(q);
		isNotInQueue[nextnode] = TRUE;

		// node elements are in R_node[[1]], weights are in R_node[[2]]
		R_node = VECTOR_ELT(NETWORK, nextnode);
		weights = REAL(VECTOR_ELT(R_node, 1));
		R_node = VECTOR_ELT(R_node, 0);
		num_connected = LENGTH(R_node);
		if(!num_connected){
			if(++subiter == n){
				subiter = 0;
				iter++;
			}
			continue;
		}

		connectedNodes = INTEGER(R_node);
		new_group = tabulate(connectedNodes, weights, clusters, num_connected);

		if(new_group != clusters[nextnode]){
			// add all neighbors in different clusters to the queue (if not already in)
			clusters[nextnode] = new_group;
			for(int i=0; i < num_connected; i++){
				nextnode = connectedNodes[i]-1;
				if(clusters[nextnode] != new_group && isNotInQueue[nextnode]){
					enqueue(q, nextnode+1);
					isNotInQueue[nextnode] = FALSE;
				}
			}
		}

		if(++subiter == n){
			subiter = 0;
			iter++;
			R_CheckUserInterrupt();
		}
	}

	// now the total clustering is stored in `clusters`
	SEXP RETVAL = PROTECT(allocVector(INTSXP, n));
	int *rvptr = INTEGER(RETVAL);
	for(int i=0; i<n; i++) rvptr[i] = clusters[i];

	free(isNotInQueue);
	free(clusters);
	free_q(q);

	UNPROTECT(1);
	return(RETVAL);
}

SEXP R_convertgraph(SEXP V1, SEXP V2, SEXP WEIGHT, SEXP NROW, SEXP NVERT){
	// we could easily change this to read/write to files instead of work in-memory
	// ensure that V1, V2 are 0-indexed
	int *v_from = INTEGER(V1);
	int *v_to = INTEGER(V2);
	double *weights = REAL(WEIGHT);
	int n_edges = INTEGER(NROW)[0];
	int n_vert = INTEGER(NVERT)[0];

	int from_v, to_v;
	double w;

	// count the number of neighbors for each node
	int *vert_edges_count = calloc(n_vert, sizeof(int));
	for(int i=0; i<n_edges; i++){
		vert_edges_count[v_from[i]]++;
		vert_edges_count[v_to[i]]++;
	}

	// initialize all the node objects
	t_gnode *graph = malloc(n_vert*sizeof(t_gnode));
	for(int i=0; i<n_vert; i++){
		graph[i].neighbors = malloc(vert_edges_count[i]*sizeof(int));
		graph[i].weights = malloc(vert_edges_count[i]*sizeof(double));
	}

	int *tmpcount = malloc(n_vert*sizeof(int));
	memcpy(tmpcount, vert_edges_count, n_vert*sizeof(int));

	for(int i=0; i<n_edges; i++){
		from_v = v_from[i];
		to_v = v_to[i];
		w = weights[i];

		graph[from_v].neighbors[--tmpcount[from_v]] = to_v;
		graph[to_v].neighbors[--tmpcount[to_v]] = from_v;

		graph[from_v].weights[tmpcount[from_v]] = w;
		graph[to_v].weights[tmpcount[to_v]] = w;
	}
	free(tmpcount);

	SEXP OUT_GRAPH = PROTECT(allocVector(VECSXP, n_vert));
	SEXP NODELIST, EDGELIST, ELEM;
	int nct;
	for(int i=0; i<n_vert; i++){
		nct = vert_edges_count[i];
		ELEM = PROTECT(allocVector(VECSXP, 2));
		NODELIST = PROTECT(allocVector(INTSXP, nct));
		EDGELIST = PROTECT(allocVector(REALSXP, nct));
		memcpy(INTEGER(NODELIST), graph[i].neighbors, nct*sizeof(int));
		memcpy(REAL(EDGELIST), graph[i].weights, nct*sizeof(double));
		SET_VECTOR_ELT(ELEM, 0, NODELIST);
		SET_VECTOR_ELT(ELEM, 1, EDGELIST);
		SET_VECTOR_ELT(OUT_GRAPH, i, ELEM);
		UNPROTECT(3);
		free(graph[i].neighbors);
		free(graph[i].weights);
	}
	free(graph);
	free(vert_edges_count);
	UNPROTECT(1);
	return(OUT_GRAPH);
}

SEXP R_fastcount(SEXP LVEC, SEXP NVERT, SEXP OFFSET){
	int *lv = INTEGER(LVEC);
	int mod = INTEGER(NVERT)[0];
	int l = LENGTH(LVEC);
	int offset = INTEGER(OFFSET)[0]-1;
	int to_check = -1;

	SEXP OUTVEC = PROTECT(allocVector(INTSXP, mod));
	int *outvec = INTEGER(OUTVEC);
	memset(outvec, 0, mod*sizeof(int));
	for(int i=0; i<l; i++){
		if(!(i%mod)) to_check = lv[i+offset];
		outvec[i%mod] += lv[i]==to_check;
	}

	UNPROTECT(1);
	return(OUTVEC);
}
