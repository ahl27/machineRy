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

SEXP test_writing(){
	l_uint ctr = pow(2,32) + 1;
	char dirpath[] = "./tmp_clustering_scratch";
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

SEXP test_batchnamehash(SEXP FILENAME, SEXP NUM_EFILES, SEXP TEMPTABNAME, SEXP OUTDIR,
                        SEXP SEPS, SEXP CTR, SEXP VERBOSE, SEXP IS_UNDIRECTED){
  /*
   * reads in values and checks if we correctly hashed them
   * .Call("test_batchnamehash", edgefiles, length(edgefiles), (tf <- tempfile()), (td <- tempdir()), "\t\n", 0, TRUE, TRUE)
   * .Call("test_batchnameverify", list.files(td, full.names=TRUE), length(list.files(td)))
   */

  const char* dir = CHAR(STRING_ELT(OUTDIR, 0));
  const char* edgefile;
  const char* temptabfile = CHAR(STRING_ELT(TEMPTABNAME, 0));
  const char* seps = CHAR(STRING_ELT(SEPS, 0));
  const int num_edgefiles = INTEGER(NUM_EFILES)[0];
  const int verbose = LOGICAL(VERBOSE)[0];
  const int is_undirected = LOGICAL(IS_UNDIRECTED)[0];
  l_uint num_v = (l_uint)(REAL(CTR)[0]);

  // first, index all vertex names
  for(int i=0; i<num_edgefiles; i++){
    edgefile = CHAR(STRING_ELT(FILENAME, i));
    num_v += hash_file_vnames_batch(edgefile, dir, temptabfile, seps[0], seps[1], num_v, verbose, is_undirected);
  }

  // set up the table file with 0s for all vertices
  // have to while loop this, not guaranteed to write this many values in one call
  l_uint v = 0;
  l_uint num_wrote = num_v;
  FILE *f = fopen(temptabfile, "ab");
  while(num_wrote > 0){
    num_wrote -= fwrite(&v, L_SIZE, num_wrote, f);
  }
  fclose(f);

  // next, record how many edges each has
  for(int i=0; i<num_edgefiles; i++){
    edgefile = CHAR(STRING_ELT(FILENAME, i));
    get_edge_counts_batch(edgefile, dir, temptabfile, seps[0], seps[1], verbose, is_undirected);
  }

  return(R_NilValue);
}

SEXP test_batchnamehashverify(SEXP HASHEDFILES, SEXP NUM_FILES, SEXP NUM_NODES){
  const int nfiles = INTEGER(NUM_FILES)[0];
  const char* fname;
  uint16_t name_len;
  int LEN_SIZE = sizeof(uint16_t);
  int nnodes = INTEGER(NUM_NODES)[0];
  char buf[MAX_NODE_NAME_SIZE];
  char *indic = malloc(nnodes);
  l_uint index;
  FILE *f;

  memset(indic, 0, nnodes);
  for(int i=0; i<nfiles; i++){
    fname = CHAR(STRING_ELT(HASHEDFILES,i));
    f = fopen(fname, "rb");
    while(fread(&name_len, LEN_SIZE, 1, f)){
      // read in data from node names
      memset(buf, '\0', MAX_NODE_NAME_SIZE);
      safe_fread(buf, 1, name_len, f);
      safe_fread(&index, L_SIZE, 1, f);
      index--;
      //Rprintf("Node %s received counter %llu\n", buf, index);
      indic[index]++;
    }
    fclose(f);
  }

  int failed = 0;
  for(int i=0; i<nnodes; i++){
    if(indic[i] == 0){
      Rprintf("No node assigned to %d\n", i);
      failed = 1;
    }
    if(indic[i] > 1){
      Rprintf("Multiple nodes assigned to %d\n", i);
      failed = 1;
    }
  }

  if(!failed) Rprintf("All nodes assigned to unique labels!\n");
  return R_NilValue;
}

/*****************/
/* Old functions */
/*****************/

l_uint rw_vertname(const char *vname, const char *dir, l_uint ctr){
	/*
	 * ctr is used so we don't have to duplicate too much code
	 * if ctr==0 we search for the vertex, else we try to write it
	 * there is a concern with endianness here that we will figure out later
	 *
	 * Update: Going to write each value as LENGTH STRING NUM
	 *  This way we can just compare the lengths to check for equality
	 *	if the lengths are different, just fseek forward (length + L_SIZE)
	 *	else compare.
	 *  Since the string is assumed to be at most 256 char, it will fit into an uint16_t
	 */
	//const char *vtmp = vname;
	const uint16_t vname_len = strlen(vname);
	const int LEN_SIZE = sizeof(uint16_t);

	uint16_t read_len;
	uint hash = hash_string_fnv(vname);
	l_uint tmp_ctr;

	char read_name[MAX_NODE_NAME_SIZE];
	char fname[PATH_MAX];

	// build filename using the hash
	snprintf(fname, strlen(dir) + 6, "%s/%04x", dir, hash);
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
		error("%s", "Error opening file for writing\n");
	}

	// returns a size_t of the number of elements read
	while(fread(&read_len, LEN_SIZE, 1, f)){
		if(read_len != vname_len){
			fseek(f, read_len+L_SIZE, SEEK_CUR);
		} else {
			memset(read_name, '\0', MAX_NODE_NAME_SIZE);
			tmp_ctr = 0;
			safe_fread(read_name, 1, read_len, f);
			safe_fread(&tmp_ctr, L_SIZE, 1, f);
			if(!strcmp(vname, read_name)){
				fclose(f);
				return tmp_ctr;
			}
		}
	}

	// if we made it here, we haven't found the vertex and we're at the end of the file
	if(!ctr){
		// didn't find in read-only mode, return 0
		fclose(f);
		return ctr;
	}

	// else write the string
	fwrite(&vname_len, LEN_SIZE, 1, f);
	fwrite(vname, 1, vname_len, f);
	fwrite(&ctr, L_SIZE, 1, f);
	fclose(f);

	return ctr+1;
}

l_uint hash_file_vnames(const char* fname, const char* dname, const char *ftable,
	const char sep, const char line_sep, l_uint ctr, int v, int is_undirected){
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
	l_uint print_counter = 0;

	if(v) Rprintf("Reading file %s...\n", fname);

	while(!feof(f)){
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

				if(feof(f)) errorclose_file(f, tab, "Unexpected end of file.\n");
			}

			vname[cur_pos] = '\0';
			c = getc(f);


			found_vert = rw_vertname(vname, dname, cur_ctr);
			// if directed, only increment the outgoing edge
			if(found_vert > cur_ctr){
				fseek(tab, 0, SEEK_END);
				num_edges = is_undirected ? 1 : 1-iter;
				cur_ctr = found_vert;
			} else {
				if(!is_undirected && iter == 1) continue;
				fseek(tab, (found_vert-1)*L_SIZE, SEEK_SET);
				safe_fread(&num_edges, L_SIZE, 1, tab);
				num_edges++;
				fseek(tab, -1*L_SIZE, SEEK_CUR);
			}

			fwrite(&num_edges, L_SIZE, 1, tab);
		}

		while(c != line_sep && !feof(f)) c = getc(f);
		if(c == line_sep) c=getc(f);
		print_counter++;
		if(!(print_counter % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%lu lines read\r", print_counter);
			else R_CheckUserInterrupt();
		}
	}
	if(v) Rprintf("\t%lu lines read\n\t%lu total nodes.\n", print_counter, cur_ctr-1);
	fclose(f);
	fclose(tab);
	return cur_ctr-1;
}

int read_edge_to_table(FILE *edgefile, FILE *mastertab, FILE *countstab, const char* hashdir,
												const char sep, const char linesep, size_t entrysize, l_uint num_v, int is_undirected, const int self_loop_inc){
	/*
	 * function reads in a single edge
	 * should assume that the edgefile pointer is already at the correct location
	 * other filepointers have undefined locations
	 */
	char tmp[MAX_NODE_NAME_SIZE], c;
	uint ctr=0;
	l_uint indices[2], offset[2], locs[2];
	double weight;
	int itermax = is_undirected ? 2 : 1;

	// index both the names
	for(int i=0; i<2; i++){
		memset(tmp, '\0', MAX_NODE_NAME_SIZE);
		ctr = 0;
		c = getc(edgefile);

		// return if we get to the end of the file -- simplest way to do this
		if(feof(edgefile)) return 0;

		while(c != sep){
			tmp[ctr++] = c;
			c = getc(edgefile);
		}
		indices[i] = rw_vertname(tmp, hashdir, 0)-1;

		// get offset for location we'll write to in the counts file
		fseek(countstab, (indices[i])*L_SIZE, SEEK_SET);
		safe_fread(&offset[i], L_SIZE, 1, countstab);

		// get start location in table file we'll write to
		fseek(mastertab, (indices[i])*L_SIZE, SEEK_SET);
		safe_fread(&locs[i], L_SIZE, 1, mastertab);
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
	for(int i=0; i<itermax; i++){
		// note if !is_undirected then we only write from->to direction
		offset[i]--;
		fseek(mastertab, (num_v+1)*L_SIZE, SEEK_SET);
		fseek(mastertab, (locs[i]+offset[i]+self_loop_inc)*entrysize, SEEK_CUR);
		fwrite(&indices[(i+1)%2], L_SIZE, 1, mastertab);
		fwrite(&weight, sizeof(double), 1, mastertab);

		// decrement the counts file
		fseek(countstab, indices[i]*L_SIZE, SEEK_SET);
		fwrite(&offset[i], L_SIZE, 1, countstab);
	}

	return 1;
}

void csr_compress_edgelist(const char* edgefile, const char* dname, const char* curcountfile, const char* ftable,
														const char sep, const char linesep, l_uint num_v, int v, int is_undirected, int has_self_loops){
	/*
	 * This should be called after we've already read in all our files
	 * critically, ensure we're rewritten our ftable file such that it is cumulative counts and not vertex counts
	 */
	const size_t entry_size = L_SIZE + sizeof(double);
	const int self_loop_inc = has_self_loops ? 1 : 0;

	FILE *mastertable, *tmptable, *edgelist;
	mastertable = fopen(ftable, "rb+");
	if(!mastertable){
		mastertable = fopen(ftable, "ab+");
		fclose(mastertable);
		mastertable = fopen(ftable, "rb+");
		if(!mastertable) error("%s", "error opening temporary counts file.\n");
	}

	tmptable = fopen(curcountfile, "rb+");
	if(!tmptable) errorclose_file(mastertable, NULL, "error opening master table file.\n");

	edgelist = fopen(edgefile, "rb");
	if(!edgelist) errorclose_file(tmptable, mastertable, "error opening edgelist file.\n");

	int status = 1;
	l_uint print_counter = 0;
	if(v) Rprintf("Reading edges from file %s...\n", edgefile);

	// super slow, need to read edges in batches
	while(status){
		// call unique_strings_with_sideeffects and then find_node_indices_batch
		// create two other caches to hold indices and weights
		// convert each edge to its index in cache1
		// loop over and add all edges (both directions if undirected)
		status = read_edge_to_table(edgelist, mastertable, tmptable, dname, sep, linesep, entry_size, num_v, is_undirected, self_loop_inc);
		print_counter++;
		if(!(print_counter % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%lu edges read\r", print_counter);
			else R_CheckUserInterrupt();
		}
	}
	print_counter--;
	if(v) Rprintf("\t%lu edges read\n", print_counter);
	fclose(mastertable);
	fclose(tmptable);
	fclose(edgelist);
	return;
}


void simple_luint_fcopy(const char* f1, const char* f2, int mode, l_uint mod){
	// if add_index, add the index into the file to sort
	// otherwise reindex the sorted values and flip the values (for later)
	double_lu dlu = {0,0};
	l_uint prev_ind=0, ctr=0;
	int nread;
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");

	switch(mode){
	case 0: // write and add index
		while(1){
			dlu.ctr1++;
			nread = fread(&dlu.ctr2, L_SIZE, 1, orig);
			if(!nread) break;
			fwrite(&dlu, sizeof(double_lu), 1, copy);
		}
		break;
	case 1: // reindex clusters and write
		while(1){
			nread = fread(&dlu, sizeof(double_lu), 1, orig);
			if(prev_ind != dlu.ctr2){
				prev_ind = dlu.ctr2;
				dlu.ctr2 = ++ctr;
			} else {
				dlu.ctr2 = ctr;
			}
			if(!nread) break;
			fwrite(&dlu, sizeof(double_lu), 1, copy);
		}
		break;
	case 2: // read structs, write flipped (struct becomes index, cluster)
		while(1){
			nread = fread(&dlu, sizeof(double_lu), 1, orig);
			if(!nread) break;
			prev_ind = dlu.ctr1;
			dlu.ctr1 = dlu.ctr2;
			dlu.ctr2 = prev_ind;
			fwrite(&dlu, sizeof(double_lu), 1, copy);
		}
		break;
	case 3: // read structs, write only second value (cluster)
		ctr = 1;
		while(1){
			nread = fread(&dlu, sizeof(double_lu), 1, orig);
			if(!nread) break;
			if(ctr){
				prev_ind = mod - dlu.ctr1;
				ctr = 0;
			}
			dlu.ctr1 = ((dlu.ctr1 + prev_ind) % mod) + 1;
			fwrite(&dlu.ctr1, L_SIZE, 1, copy);
		}
		break;
	}

	fclose(orig);
	fclose(copy);
	return;
}

void mergesort_clust_file(const char* f, const char* dir, int firstpass){
	/*
	 * reindex the clustering file with mergesort
	 * first pass: affix an index and sort as {clust, ind}
	 * then swap the values so that we have {ind, clust}
	 * second pass: re-sort the values and drop the indices
	 */

	// two read pointers, one write pointer
	FILE *f1_r1, *f1_r2, *f2_w;
	char file1[PATH_MAX], file2[PATH_MAX];
	char *finalfile;
	size_t dlu_size = sizeof(double_lu);

	// create the junk files we'll use
	safe_filepath_cat(dir, "msfile1", file1, PATH_MAX);
	safe_filepath_cat(dir, "msfile2", file2, PATH_MAX);

	// first, we'll use the cache to read in preprocessed sorted blocks of size BYTES_PER_ENTRY
	l_uint block_size = FILE_READ_CACHE_SIZE;
	l_uint total_lines = 0;
	l_uint max_clust = 0;
	double_lu read_cache[FILE_READ_CACHE_SIZE];

	// copy the original file into file 1
	// on first pass append counters
	// on second pass, flip the counters
	simple_luint_fcopy(f, file1, firstpass ? 0 : 2, max_clust);

	// open file, read in chunks, sort locally, write to file
	l_uint cachectr = 0;
	f1_r1 = fopen(file1, "rb");
	f2_w = fopen(file2, "wb");
	while(fread(&read_cache[cachectr++], dlu_size, 1, f1_r1)){
		if(!firstpass){
			if(read_cache[cachectr-1].ctr1 > max_clust)
				max_clust = read_cache[cachectr-1].ctr1;
		}
		total_lines++;
		if(cachectr == block_size){
			qsort(read_cache, cachectr, dlu_size, l_uint_compar);
			for(int i=0; i<cachectr; i++)
				fwrite(&read_cache[i], dlu_size, 1, f2_w);
			cachectr=0;
		}
	}
	if(cachectr){
		qsort(read_cache, cachectr, dlu_size, l_uint_compar);
		for(int i=0; i<cachectr; i++)
			fwrite(&read_cache[i], dlu_size, 1, f2_w);
	}

	fclose(f1_r1);
	fclose(f2_w);
	finalfile = file2;


	l_uint cur_lines = 0;
	int iter1, iter2, previt1, previt2;
	int flip = 0, cmp;
	double_lu tmp1, tmp2;
	char *f1, *f2;
	while(block_size < total_lines){
		// f1 is always the reading file, f2 the writing file
		if(flip){
			f1 = file1;
			f2 = file2;
		} else {
			f1 = file2;
			f2 = file1;
		}

		flip = !flip;

		f1_r1 = fopen(f1, "rb");
		f1_r2 = fopen(f1, "rb");
		f2_w = fopen(f2, "wb");
		// move second pointer forward to second block
		fseek(f1_r2, dlu_size*block_size, SEEK_CUR);

		// sort file 1 into file 2
		while(cur_lines < total_lines){
			// sort one block from file 1 into file 2 -- iter stores lines remaining
			iter1=total_lines - cur_lines;
			if(iter1 > block_size) iter1 = block_size;
			cur_lines += iter1;

			iter2=total_lines - cur_lines;
			if(iter2 > block_size) iter2 = block_size;
			cur_lines += iter2;

			previt1=iter1+1;
			previt2=iter2+1;
			while(iter1 || iter2){
				// there is a bug here
				if(iter1 && iter1 != previt1){
					safe_fread(&tmp1, dlu_size, 1, f1_r1);
					previt1 = iter1;
				}
				if(iter2 && iter2 != previt2){
					safe_fread(&tmp2, dlu_size, 1, f1_r2);
					previt2 = iter2;
				}

				cmp = l_uint_compar(&tmp1, &tmp2);
				if(iter1 && (!iter2 || cmp <= 0 )){
					fwrite(&tmp1, dlu_size, 1, f2_w);
					iter1--;
				} else {
					fwrite(&tmp2, dlu_size, 1, f2_w);
					iter2--;
				}
			}
			// advance pointers one block past where we just read:
			// if we move too far it doesn't really matter, we'll catch it on the next part
			fseek(f1_r1, dlu_size*block_size, SEEK_CUR);
			fseek(f1_r2, dlu_size*block_size, SEEK_CUR);
		}

		fclose(f1_r1);
		fclose(f1_r2);
		fclose(f2_w);
		cur_lines = 0;
		block_size *= 2;
		finalfile = f2;
	}
	simple_luint_fcopy(finalfile, f, firstpass ? 1 : 3, max_clust);

	return;
}

// old reformat_clusters

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

void reformat_clusters(FILE *clusterfile, l_uint num_v){
	l_uint tmp_cind;
	ll *head = malloc(sizeof(ll));
	head->next=NULL; // other values can just be garbage

	rewind(clusterfile);
	for(l_uint i=0; i<num_v; i++){
		safe_fread(&tmp_cind, L_SIZE, 1, clusterfile);
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