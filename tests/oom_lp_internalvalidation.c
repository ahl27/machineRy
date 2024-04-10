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
