/*
 * Storing graphs in CSR format
 * Weighted graph will have associated files:
 * 	1. Index: n'th entry contains the start location of information for vertex n
 *						Length n+1, last entry is the end of the file (fixed width int)
 *	2. Connections: Contains vertices the n'th vertex is connected to.
 *									Indexed by `Index` file. (fixed width int)
 *	3. Weights: Same as `Connections`, but weights (double)
 *
 * Example:
 *		graph: {0,1,2,3}
 *    edges: {(0,1,0.1), (0,2,0.2), (1,3,0.3), (2,3,0.4)}
 *
 * 				  Index:  0 2 4 6 8
 *		Connections:  1   2   0   3   0   3   1   2
 *        Weights: 0.1 0.2 0.1 0.3 0.2 0.4 0.3 0.4
 *
 * Algorithm:
 *	1. count number of edges for each node
 *		- Probably can't be done in-memory
 *		- Initialize two files (conn, conn_tmp) to all zeros, r/w to first increment over time
 *	2. loop over edges again, writing neighbors and weights
 *		- for each edge, write info to file[conn[n]+conn_tmp]
 *		- increment conn_tmp at both vertices
 *
 * Plan for data:
 * 	- make each function / execution pipeline work on a single file
 *  - out of memory means that we don't have to process everything at once
 *	- for edgelist files in a directory, call C for each from R
 *	- do it again for other files
 */

#include "machineRy.h"

#define uint uint32_t
#define l_uint uint64_t

/*
 * common limits are defined in limits.h
 * 	PAGE_SIZE: size in bytes of a page
 *   PATH_MAX: max size in bytes of a path name
 */

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif

#define MAX_NODE_NAME_SIZE 256

const char DELIM = 23; // 23 = end of transmission block, not really used for anything nowadays
const int L_SIZE = sizeof(l_uint);
const l_uint MAX_EDGES_EXACT = 5000; // this is a soft cap -- if above this, we sample edges probabalistically

typedef struct ll {
	l_uint id;
	double w;
	struct ll* next;
} ll;

ll* insert_ll(ll* head, l_uint id, double w){
	ll *tmp = head;
	if(!tmp){
		tmp = malloc(sizeof(ll));
		tmp->id = id;
		tmp->w = w;
		tmp->next=NULL;
		return(tmp);
	}
	while(tmp->id != id && tmp->next) tmp = tmp->next;
	if(tmp->id!=id){
		tmp->next = malloc(sizeof(ll));
		tmp = tmp->next;
		tmp->id = id;
		tmp->w = w;
		tmp->next = NULL;
	} else {
		tmp->w += w;
	}

	return head;
}

l_uint indexed_insert(ll *head, l_uint id){
	// head is always just going to be dummy start pos
	ll* tmp = head;
	l_uint ctr = 0;
	while(tmp->next && tmp->next->id != id){
		tmp = tmp->next;
		ctr++;
	}

	// two scenarios: 1) tmp->next is NULL, didn't find; 2) tmp->next is the id
	if(!tmp->next){
		tmp->next = malloc(sizeof(ll));
		tmp = tmp->next;
		tmp->id = id;
		tmp->next = NULL;
	}

	return ctr+1;
}


void errorclose_file(FILE *f1, FILE *f2, const char* message){
	fclose(f1);
	if(f2) fclose(f2);
	error(message);
}

uint hash_string_fnv(const char *str){
	/*
	 * this is a Fowler-Noll-Vo hash function, it's fast and simple -- see wikipedia for constants
	 * https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
	 *
	 * hashes to a 32-bit uint, so we can have up to 65,536 files in the folder
	*/
	const uint fnv_prime = 0x01000193;
	uint hash = 0x811c9dc5;
	const char *s = str;
	char tmp;

	while(*s){
		tmp = *s++;
		if(tmp != DELIM){
			// overflow here is fine
			hash *= fnv_prime;
			hash ^= tmp;
		}
	}

	return hash;
}

l_uint rw_vertname(const char *vname, const char *dir, l_uint ctr){
	/*
	 * ctr is used so we don't have to duplicate too much code
	 * if ctr==0 we search for the vertex, else we try to write it
	 * there is a concern with endianness here that we will figure out later
	 */
	char fname[PATH_MAX];
	const char *vtmp = vname;
	uint hash = hash_string_fnv(vname);
	char found_name, c;

	// build filename using the hash
	sprintf(fname, "%s/%08x", dir, hash);
	FILE *f;
	if(ctr){
		// create file if it doesn't exist
		f = fopen(fname, "ab+");
		fclose(f);
		f = fopen(fname, "rb+");
	}
	else
		f = fopen(fname, "rb");


	// file doesn't exist or failed to open
	if(!f){
		if(!ctr) return 0;
		error("Error opening file for writing\n");
	}

	c = getc(f);
	while(c != EOF){
		// read in the vertex name, comparing as we go
		vtmp = vname;
		found_name = 1;
		while(c != DELIM){
			if(found_name)
				found_name = c == *vtmp++;
			c = getc(f);
		}

		// edge case: written name is shorter than search name (e.g. "test", "testtest")
		found_name = found_name && !(*vtmp);

		// if we found the name, just return
		if(found_name){
			fread(&ctr, sizeof(l_uint), 1, f);
			fclose(f);
			return ctr;
		}

		// if we didn't find it, c==DELIM and we have to skip the ctr (64 bit number)
		fseek(f, 8, SEEK_CUR);
		c = getc(f);
	}

	// if we made it here, we haven't found the vertex and we're at the end of the file
	if(!ctr){
		// didn't find in read-only mode, return 0
		fclose(f);
		return ctr;
	}

	c = DELIM;
	fwrite(vname, sizeof(char), strlen(vname), f);
	fwrite(&c, sizeof(char), 1, f);
	fwrite(&ctr, L_SIZE, 1, f);
	fclose(f);

	return ctr+1;
}

l_uint hash_file_vnames(const char* fname, const char* dname, const char *ftable, const char sep, const char line_sep, l_uint ctr){
	/*
	 * fname: .tsv list of edges
	 * dname: directory of hash codes
	 * ftable: output counts file, should be L_SIZE*num_v bytes
	 */
	FILE *f = fopen(fname, "rb");
	FILE *tab = fopen(ftable, "wb+");
	char vname[MAX_NODE_NAME_SIZE];
	int cur_pos = 0;
	char c = getc(f);
	l_uint found_vert = 0, cur_ctr = ctr, num_edges = 0;

	while(c != EOF){
		// going to assume we're at the beginning of a line
		// lines should be of the form `start end weight`
		// separation is by char `sep`
		for(int iter=0; iter<2; iter++){
			cur_pos = 0;
			while(c != sep){
				vname[cur_pos++] = c;
				c = getc(f);
				if(cur_pos >= MAX_NODE_NAME_SIZE)
					errorclose_file(f, tab, "Node name is larger than max allowed name size.\n");

				if(c == EOF)
					errorclose_file(f, tab, "Unexpected end of file.\n");
			}

			vname[cur_pos] = '\0';
			c = getc(f);
			found_vert = rw_vertname(vname, dname, cur_ctr);
			if(found_vert > cur_ctr){
				fseek(tab, 0, SEEK_END);
				num_edges = 1;
				cur_ctr = found_vert;
			} else {
				fseek(tab, (found_vert-1)*L_SIZE, SEEK_SET);
				fread(&num_edges, L_SIZE, 1, tab);
				num_edges++;
				fseek(tab, -1*L_SIZE, SEEK_CUR);
			}
			fwrite(&num_edges, L_SIZE, 1, tab);
		}

		while(c != line_sep && c != EOF) c = getc(f);
		if(c == line_sep) c=getc(f);
	}

	fclose(f);
	fclose(tab);
	return cur_ctr-1;
}

int read_edge_to_table(FILE *edgefile, FILE *mastertab, FILE *countstab, const char* hashdir,
												const char sep, const char linesep, uint entrysize, l_uint num_v){
	/*
	 * function reads in a single edge
	 * should assume that the edgefile pointer is already at the correct location
	 * other filepointers have undefined locations
	 */
	char tmp[MAX_NODE_NAME_SIZE], c;
	uint ctr=0;
	l_uint indices[2], offset[2], locs[2];
	double weight;

	// index both the names
	for(int i=0; i<2; i++){
		ctr = 0;
		c = getc(edgefile);

		// return if we get to the end of the file -- simplest way to do this
		if(c==EOF) return 0;

		while(c != sep){
			tmp[ctr++] = c;
			c = getc(edgefile);
		}
		tmp[ctr] = '\0';
		indices[i] = rw_vertname(tmp, hashdir, 0)-1;

		// get offset for location we'll write to in the counts file
		fseek(countstab, (indices[i])*L_SIZE, SEEK_SET);
		fread(&offset[i], L_SIZE, 1, countstab);
		offset[i]--;
		fseek(countstab, -1*L_SIZE, SEEK_CUR);
		fwrite(&offset[i], L_SIZE, 1, countstab);

		// get start location in table file we'll write to
		fseek(mastertab, (indices[i])*L_SIZE, SEEK_SET);
		fread(&locs[i], L_SIZE, 1, mastertab);
	}

	// read in the weights
	ctr = 0;
	c = getc(edgefile);
	while(c != linesep){
		tmp[ctr++] = c;
		c = getc(edgefile);
	}
	tmp[ctr] = '\0';
	weight = atof(tmp);

	// write to file
	for(int i=0; i<2; i++){
		fseek(mastertab, (num_v+1)*L_SIZE, SEEK_SET);
		fseek(mastertab, (locs[i]+offset[i])*entrysize, SEEK_CUR);
		fwrite(&indices[(i+1)%2], L_SIZE, 1, mastertab);
		fwrite(&weight, sizeof(double), 1, mastertab);
	}

	return 1;
}

void reformat_counts(const char* curcounts, const char* mastertable, l_uint n_vert){
	/*
	 * Creates a new table with cumulative counts
	 * leaves the old table unchanged, this will act as a temporary counts file later
	 */
	const uint l_size = L_SIZE;
	l_uint cumul_total = 0, curcount;
	FILE *tmptab = fopen(curcounts, "rb");
	FILE *mtab = fopen(mastertable, "wb+");

	for(l_uint i=0; i<n_vert; i++){
		fwrite(&cumul_total, l_size, 1, mtab);
		fread(&curcount, l_size, 1, tmptab);
		cumul_total += curcount;
	}

	// ending position of file
	fwrite(&cumul_total, l_size, 1, mtab);

	fclose(tmptab);
	fclose(mtab);
	return;
}

void csr_compress_edgelist(const char* edgefile, const char* dname, const char* curcountfile, const char* ftable,
														const char sep, const char linesep, l_uint num_v){
	/*
	 * This should be called after we've already read in all our files
	 * critically, ensure we're rewritten our ftable file such that it is cumulative counts and not vertex counts
	 */
	const uint entry_size = L_SIZE + sizeof(double);
	FILE *mastertable, *tmptable, *edgelist;
	mastertable = fopen(ftable, "rb+");
	if(!mastertable){
		mastertable = fopen(ftable, "ab+");
		fclose(mastertable);
		mastertable = fopen(ftable, "rb+");
		if(!mastertable) error("error opening temporary counts file.\n");
	}

	tmptable = fopen(curcountfile, "rb+");
	if(!tmptable) errorclose_file(mastertable, NULL, "error opening master table file.\n");

	edgelist = fopen(edgefile, "rb");
	if(!edgelist) errorclose_file(tmptable, mastertable, "error opening edgelist file.\n");

	int status = 1;
	while(status){
		status = read_edge_to_table(edgelist, mastertable, tmptable, dname, sep, linesep, entry_size, num_v);
	}

	fclose(mastertable);
	fclose(tmptable);
	fclose(edgelist);

	return;
}

l_uint update_node_cluster(l_uint ind, l_uint offset, FILE *mastertab, FILE *clusterings){
	/*
	 * Determine number of edges using the table file (next - cur)
	 * If number of edges too large, use some sort of hash to bin edges, then rerun with less
	 * Inputs:
	 * 	- 	      ind: 0-indexed vertex id
	 * 	-      offset: location where edges begin in mastertab
	 *	-   mastertab: file pointer to CSR-compressed graph
	 *	- clusterings: file of current cluster numbers (0=unassigned)
	 */

	l_uint start, end, num_edges, tmp_cl, tmp_id, zeromaxid=ind+1;
	double tmp_w, acceptance_prob, zeromax=0;

	// move to information for the vertex and read in number of edges
	fseek(mastertab, L_SIZE*ind, SEEK_SET);
	fread(&start, L_SIZE, 1, mastertab);
	fread(&end, L_SIZE, 1, mastertab);


	// if we're above the max edges, subsample the edges we read
	num_edges = end - start;
	// if it has no edges we can't do anything
	if(num_edges == 0) return ind+1;

	acceptance_prob = fmin((double)MAX_EDGES_EXACT / num_edges, 1.0);
	ll *searchpath = NULL;
	GetRNGstate();
	fseek(mastertab, L_SIZE*offset, SEEK_SET);
	fseek(mastertab, (L_SIZE+sizeof(double))*start, SEEK_CUR);
	for(int i=0; i<num_edges; i++){
		// skip with probability
		if(unif_rand() > acceptance_prob){
			fseek(mastertab, L_SIZE+sizeof(double), SEEK_CUR);
			continue;
		}

		// read in neighbor and weight
		fread(&tmp_id, L_SIZE, 1, mastertab);
		fread(&tmp_w, sizeof(double), 1, mastertab);

		// get which cluster it belongs to
		fseek(clusterings, L_SIZE*tmp_id, SEEK_SET);
		fread(&tmp_cl, L_SIZE, 1, clusterings);

		/*
		 * this solves a special edge case where uninitialized nodes with a self-loop
		 * can be counted incorrectly if their neighbors have already been assigned to
		 * their node. Thus, if it's a self loop and the cluster is uninitialized,
		 * set it to its own cluster (what it would be initialized to)
		 */
		if(tmp_id == ind && tmp_cl == 0)
			tmp_cl = ind+1;

		// store value
		if(tmp_cl){
			searchpath = insert_ll(searchpath, tmp_cl, tmp_w);
		} else if(tmp_w > zeromax) {
			zeromax = tmp_w;
			zeromaxid = tmp_id+1; // +1 because vertices 0-indexed and clusters 1-indexed
		}
	}

	// now we've read in all the edges, find the new community
	// (and free along the way)
	// max weight and cluster will be stored in tmp_w and tmp_cl (resp.)
	ll *tmp_ll = searchpath;
	tmp_w = -1;
	while(tmp_ll){
		searchpath = tmp_ll;
		tmp_ll = tmp_ll->next;

		if(searchpath->w > tmp_w || (searchpath->w == tmp_w && unif_rand() < 0.5) ){
			// ties are broken randomly
			tmp_w = searchpath->w;
			tmp_cl = searchpath->id;
		}
		free(searchpath);
	}
	PutRNGstate();
	if(zeromax > tmp_w) tmp_cl = zeromaxid;

	/*
	 * now we have the cluster to choose stored in tmp_cl
	 * Edge cases:
	 * - crazy chance results in every edge being skipped
	 *		-> since zeromax=0 and zeromaxid=id+1 at start, we just cluster it with itself and continue
	 * - all clusters are initialized, searchpath ends up NULL
	 * 		-> we won't free a NULL, and since zeromax > -1 we'll set it to the strongest uninitialized connection
	 */

	// last step is to write it to the cluster file
	fseek(clusterings, L_SIZE*ind, SEEK_SET);
	fwrite(&tmp_cl, L_SIZE, 1, clusterings);

	return tmp_cl;
}

void reformat_clusters(FILE *clusterfile, l_uint num_v){
	l_uint tmp_cind;
	ll *head = malloc(sizeof(ll));
	head->next=NULL; // other values can just be garbage

	rewind(clusterfile);
	for(l_uint i=0; i<num_v; i++){
		fread(&tmp_cind, L_SIZE, 1, clusterfile);
		tmp_cind = indexed_insert(head, tmp_cind);
		fseek(clusterfile, -1*L_SIZE, SEEK_CUR);
		fwrite(&tmp_cind, L_SIZE, 1, clusterfile);
	}

	ll *tmp = head;
	while(tmp){
		head = tmp;
		tmp = tmp->next;
		free(head);
	}
	return;
}

void cluster_file(const char* mastertab_fname, const char* clust_fname, l_uint num_v){
	// temporary implementation for now, will adjust later
	// main runner function to cluster nodes
	FILE *masterfile = fopen(mastertab_fname, "rb");
	FILE *clusterfile = fopen(clust_fname, "rb+");

	for(l_uint i=0; i<num_v; i++){
		update_node_cluster(i, num_v+1, masterfile, clusterfile);
	}

	reformat_clusters(clusterfile, num_v);

	fclose(masterfile);
	fclose(clusterfile);
	return;
}

void verify_filecontents(const char *filename, l_uint num_v){
	FILE *f = fopen(filename, "rb");
	l_uint tmp;
	for(l_uint i=0; i<num_v; i++){
		fread(&tmp, L_SIZE, 1, f);
		Rprintf("%lu ", tmp);
		if(i % 20 == 19) Rprintf("\n");
	}
	Rprintf("\n");
	fclose(f);
	return;
}

void verify_edgelist(const char *filename, l_uint num_v){
	FILE *f = fopen(filename, "rb");
	l_uint tmp, all_indices[num_v+1];
	int curpos = 0, curvert=0;
	double dtmp;
	fread(all_indices, L_SIZE, (num_v+1), f);
	for(l_uint i=0; i<all_indices[num_v]; i++){
		if(curpos == all_indices[curvert]){
			Rprintf("\n%lu: ", curvert++);
		}
		curpos++;
		fread(&tmp, L_SIZE, 1, f);
		fread(&dtmp, sizeof(double), 1, f);
		Rprintf("%lu (%0.2f); ", tmp, dtmp);
	}
	Rprintf("\n");
	fclose(f);
	return;
}

void verify_clusters(const char *filename, l_uint num_v){
	// this is pretty slow, but the point is more just to demonstrate that it works
	// report by cluster once we implement re-indexing
	FILE *f = fopen(filename, "rb");
	l_uint tmp;
	for(l_uint i=0; i<num_v; i++){
		fread(&tmp, L_SIZE, 1, f);
		Rprintf("%lu: %lu\n", i, tmp);
	}
	fclose(f);
	return;
}



SEXP R_hashedgelist(SEXP FILENAME, SEXP TABNAME, SEXP TEMPTABNAME, SEXP OUTDIR, SEXP SEPS, SEXP CTR){
	/*
	 * I always forget how to handle R strings so I'm going to record it here
	 * R character vectors are STRSXPs, which is the same as a list (VECSXP)
	 * Each entry in a STRSXP is a CHARSXP, accessed with STRING_ELT(str, ind)
	 * You can get a `const char*` from a CHARSXP with CHAR()
	 * You can also re-encode with `const char* Rf_translateCharUTF8()`
	 */

	const char* dir = CHAR(STRING_ELT(OUTDIR, 0));
	const char* edgefile = CHAR(STRING_ELT(FILENAME, 0));
	const char* tabfile = CHAR(STRING_ELT(TABNAME, 0));
	const char* temptabfile = CHAR(STRING_ELT(TEMPTABNAME, 0));
	const char* seps = CHAR(STRING_ELT(SEPS, 0));
	l_uint ctr = (l_uint)(REAL(CTR)[0]), num_v;

	// first, index all vertex names and record how many edges each has
	Rprintf("Hashing file...\n");
 	num_v = hash_file_vnames(edgefile, dir, temptabfile, seps[0], seps[1], ctr);
 	Rprintf("Hashed! Total nodes is %lu\nEdge Counts:\n", num_v);

 	verify_filecontents(temptabfile, num_v);

 	// next, create an indexed table file of where data for each vertex will be located
 	// note that we can modify this later to hash multiple edge files before reformatting
 	reformat_counts(temptabfile, tabfile, num_v);


 	Rprintf("Counts reformatted!\n");
 	// then, we'll create the CSR compression of all our edges
 	csr_compress_edgelist(edgefile, dir, temptabfile, tabfile, seps[0], seps[1], num_v);
 	Rprintf("Edges read in!\n");

 	// temptabfile now becomes our clustering file
 	//verify_edgelist(tabfile, num_v);

 	Rprintf("Clustering nodes...\n");
 	cluster_file(tabfile, temptabfile, num_v);
 	verify_clusters(temptabfile, num_v);

	SEXP RETVAL = PROTECT(allocVector(REALSXP, 1));
	REAL(RETVAL)[0] = (double) num_v;
	UNPROTECT(1);
	return RETVAL;
}

SEXP test_writing(){
	l_uint ctr = pow(2,32) + 1;
	char dirpath[] = "/Users/aidan/Downloads/tmp_clustering_scratch";
	char *vnames[] = {"test", "othertest", "morevalues", "testesttest", "AAAAAA", "testC"};

	for(int i=0; i<5; i++){
		Rprintf("Writing '%s'\n", vnames[i]);
		ctr = rw_vertname(vnames[i], dirpath, ctr);
	}
	Rprintf("done!\n");

	Rprintf("reading values...\n");
	for(int i=0; i<6; i++){
		Rprintf("Reading '%s'...", vnames[i]);
		ctr = rw_vertname(vnames[i], dirpath, 0);
		if(ctr)
			Rprintf("found index %lu!\n", ctr);
		else
			Rprintf("not found!\n");
	}

	return R_NilValue;
}

SEXP test_outputs(SEXP TABFILE, SEXP NVERT){
	const char* tabfile = CHAR(STRING_ELT(TABFILE, 0));
	int nvert = INTEGER(NVERT)[0];
	SEXP RETVAL = PROTECT(allocVector(REALSXP, nvert));
	double *d = REAL(RETVAL);
	l_uint tmp;
	FILE *f = fopen(tabfile, "rb");
	for(int i=0; i<nvert; i++){
		fread(&tmp, L_SIZE, 1, f);
		d[i] = (double)tmp;
	}

	fclose(f);
	UNPROTECT(1);
	return RETVAL;
}