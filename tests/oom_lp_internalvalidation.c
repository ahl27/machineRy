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