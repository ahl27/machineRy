/*
 * Out of memory clustering with fast label propagation
 * Author: Aidan Lakshman
 *
 * This set of functions creates 5 files and a directory (using dummy names):
 *	-   OOMhashes: (directory) mapping between hash value and vertex name, along with index assignment
 * 	-  counts.bin: binary file containing edge counts for each vertex, and later cluster assignments
 * 	-  queue1.bin: out of memory queue for traversing nodes
 *  -  queue2.bin: second queue to improve r/w over jumping around queue1
 *  -     csr.bin: csr-compressed graph structure. Contains n+1 uint64_t values corresponding to vertex
 *		 						 indices, where the k'th value denotes the start position in the file for vertex k.
 *		 						 After the indices follow edge information of the form `d1 w1 d2 w2 d3 w3...`, where
 *		 						 each `d` corresponds to the destination of that edge, and each `w` the weight of
 *		 						 that edge. Indices index into this, e.g., if the first two values are 0 100 then the
 *		 						 outgoing edges from the first vertex are the first 100 edge entries (0-99).
 *  - outfile.tsv: .tsv file returned to R, contains two tab-separated columns (vertex name, cluster)
 *
 * Additional notes and TODOs:
 *	- At some point it's probably worth refactoring all file accesses into some kind of struct w/ accessors
 *  	-> something like a virtual array object that's actually r/w to disk, could be useful in future
 *		-> this implementation should use mmap (and Windows equivalent when necessary) to improve random r/w
 *  - Error checking needs to be improved
 *		-> no checks to ensure the edgelists are formatted the way the user claims
 *  - No calls to R_CheckUserInterrupt(), likely blocking for a while on large graphs
 *  - Unweighted graphs aren't supported. Unweighted graphs could save a ton of space (2x less in csr file).
 *		-> this would be a fairly big rewrite
 *  - Change the node name lookup to use Trie structures so that we get logarithmic lookup time
 *  - sizeof(char) is guaranteed to be 1 (see https://en.wikipedia.org/wiki/Sizeof). 1 is used instead to
 *		eliminate the extra function call and simplify code somewhat.
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
#define NODE_NAME_CACHE_SIZE 4096
#define FILE_READ_CACHE_SIZE 2048

const char DELIM = 23; // 23 = end of transmission block, not really used for anything nowadays
const int L_SIZE = sizeof(l_uint);
const int MAX_READ_RETRIES = 10;

// set this to 1 if we should sample edges rather than use all of them
const int use_limited_nodes = 0;
const l_uint MAX_EDGES_EXACT = 20000; // this is a soft cap -- if above this, we sample edges probabalistically
const int PRINT_COUNTER_MOD = 10;

typedef struct {
	uint len;
	char *name;
	l_uint ctr;
} full_read_line;

static void safe_fread(void *buffer, size_t size, size_t count, FILE *stream){
	size_t found_values = fread(buffer, size, count, stream);
	if(found_values != count){
		// two scenarios:

		// 1. read past the end of the file (throw error and return)
		if(feof(stream))
			error("%s", "Internal error: fread tried to read past end of file.\n");

		// 2. some undefined reading error (retry a few times and then return)
		for(int i=0; i<MAX_READ_RETRIES; i++){
			// if we read a partial value, reset the counter back
			if(found_values) fseek(stream, -1*((int)found_values), SEEK_CUR);

			// try to read again
			found_values = fread(buffer, size, count, stream);
			if(found_values == count) return;
		}

		// otherwise throw an error
		error("Internal error: fread read %zu values (expected %zu).\n", found_values, count);
	}
	return;
};

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
	error("%s", message);
}

/*
 * Functions for reading in node names
 */

uint hash_string_fnv(const char *str){
	/*
	 * this is a Fowler-Noll-Vo hash function, it's fast and simple -- see wikipedia for constants
	 * https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
	 *
	 * hashes to a 32-bit uint that we truncate to a 16-bit, so we can have up to 65,536 files in the folder
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

	// XOR fold to 16 bits
	hash = (hash & 0x0000FFFF) ^ ((hash & 0xFFFF0000) >> 16);

	/*
	 * take the lowest 4 bits -- the fewer files, the better due to batch processing
	 * this means we instead only have 16 possible files. 8 bit may end up being better long-term.
	 *
	 * The important consideration is that reading is fast, but writing and open/closing is slow
	 * Thus having a few big files is a lot faster than lots of small ones.
	 * However, having a few big files affects lookup time of rw_vertname -- likely want to batch/cache this as well
	 */
	hash &= 0xF;
	return hash;
}

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
			if(strcmp(vname, read_name) == 0){
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

int node_name_cmpfunc(const void *a, const void *b){
	// sort based on hash, then length, then strcmp
	const char *aa = *(const char **)a;
	const char *bb = *(const char **)b;

	// sort first by hash
	int v1 = hash_string_fnv(aa);
	int v2 = hash_string_fnv(bb);
	if(v1 != v2) return v1 - v2;

	// sort second by string length
	v1 = strlen(aa);
	v2 = strlen(bb);
	if(v1 != v2) return v1 - v2;

	// sort second by string comparison (this includes string length)
	return strcmp(aa, bb);
}

int nohash_name_cmpfunc(const void *a, const void *b){
	// same as above, but skip the hash comparison
	const char *aa = *(const char **)a;
	const char *bb = *(const char **)b;

	// sort first by string length
	int v1 = strlen(aa), v2 = strlen(bb);
	if (v1 != v2) return v1 - v2;

	return strcmp(aa, bb);
}

int struct_frl_cmpfunc(const void *a, const void *b){
	full_read_line aa = *(full_read_line *)a;
	full_read_line bb = *(full_read_line *)b;

	// cast to prevent unsigned arithmetic
	if(aa.len != bb.len) return ((int)aa.len) - ((int)bb.len);

	return strcmp(aa.name, bb.name);
}

void unique_strings_with_sideeffects(char **names, int num_to_sort, uint *filecounts,
																		uint *hashes, uint *NumUniqueHashes, int *InsertPoint,
																		uint *counts, int useCounts){
	/*
	 * This code is duplicated a lot, so just putting it here for consistency
	 * This uniques the set of strings **names and stores additional information
	 *	- NumUniqueHashes: number of unique hash values
	 *	-     	   hashes: unique hash values in set
	 * 	-      filecounts: number of strings per unique hash value
	 *  -     InsertPoint: number of unique strings
	 *  -          counts: number of each unique string (ignored if !useCounts)
	 */
	int insert_point = 0;
	int num_unique_hashes = 0;
	int hashctr = 1;
	uint cur_hash, tmp_hash;

	// first sort the array
	qsort(names, num_to_sort, sizeof(char*), node_name_cmpfunc);

	// next, unique the values
	hashes[0] = hash_string_fnv(names[0]);
	if(useCounts) counts[0] = 1;
	hashctr = 1;
	cur_hash = strlen(names[0]);
	for(int i=1; i<num_to_sort; i++){
		if(cur_hash != strlen(names[i]) || (strcmp(names[i], names[insert_point]) != 0)){
			// if the string is different, save it
			insert_point++;
			if(useCounts) counts[insert_point] = 1;
			if(insert_point != i) memcpy(names[insert_point], names[i], MAX_NODE_NAME_SIZE);
			cur_hash = strlen(names[insert_point]);

			// if the hash is different, increment the hash counter
			tmp_hash = hash_string_fnv(names[i]);
			if(hashes[num_unique_hashes] != tmp_hash){
				filecounts[num_unique_hashes] = hashctr;
				hashctr = 0;
				hashes[++num_unique_hashes] = tmp_hash;
			}
			hashctr++;
		} else if(useCounts) {
			// else it's the same, so increment the corresponding count
			counts[insert_point]++;
		}
	}
	filecounts[num_unique_hashes] = hashctr;
	num_unique_hashes++;
	insert_point++;

	// side effect writes
	*InsertPoint = insert_point;
	*NumUniqueHashes = num_unique_hashes;
	return;
}

int check_inputs_against_hashes(char **input_strings, char **file_strings, char *bitarray, uint ninput, uint nfile){
	// input_strings is sorted (length, strcmp), file_strings is not -- first we sort file_strings quickly
	qsort(file_strings, nfile, sizeof(char*), nohash_name_cmpfunc);

	// then we can find all matches in linear time
	// retval will equal the sum of the bitarray minus any matches, meaning 0 if everything found and positive else
	uint i=0, j=0, len1=strlen(input_strings[i]), len2=strlen(file_strings[j]);
	int cmp;

	while(i < ninput && j < nfile){
		// skip past anything we've already found
		if(!bitarray[i] || len1 < len2){
			i++;
			if(i < ninput)
				len1=strlen(input_strings[i]);
		} else if(len1 > len2){
			j++;
			if(j < nfile) len2=strlen(file_strings[j]);
		} else {
			cmp = strcmp(input_strings[i], file_strings[j]);
			if(cmp < 0){
				i++;
			} else if (cmp > 0){
				j++;
			} else {
				bitarray[i] = 0;
				i++;
				j++;
			}
		}
	}

	// return 0 if everything has already been seen, we can short circuit to skip all writes
	for(int i=0; i<ninput; i++)
	 if(bitarray[i]) return 1;
	return 0;
}

l_uint batch_write_nodes(char **names, int num_to_sort, const char *dir, l_uint ctr){
	// assuming all the nodes are in a const char* array
	const int LEN_SIZE = sizeof(uint16_t);

	// holds lines from the file
	char *cached[FILE_READ_CACHE_SIZE];
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) cached[i] = malloc(MAX_NODE_NAME_SIZE);

	// filename, has_seen bit array
	char fname[PATH_MAX], bitarray[num_to_sort];

	// unique hash values, number of entries per hash counter
	uint hashes[num_to_sort], filecounts[num_to_sort];

	// strlens
	uint16_t tmp_len, cur_len;

	// these are sometimes used for temp storage
	uint hashctr, num_unique_hashes;
	l_uint foundctr;
	int status, insert_point;
	FILE *f;

	unique_strings_with_sideeffects(names, num_to_sort, filecounts, hashes, &num_unique_hashes, &insert_point, NULL, 0);

	char **tmp_charptr = names;

	// now the new length of the array is insert_point
	// note that these are first sorted by hash, so once the hash value changes it will never appear again
	for(int i=0; i<num_unique_hashes; i++){
		// advance string array to lines for this file
		if(i) tmp_charptr = &(tmp_charptr[filecounts[i-1]]);

		// build filename using the hash
		snprintf(fname, strlen(dir) + 6, "%s/%04x", dir, hashes[i]);

		// create file if it doesn't exist
		f = fopen(fname, "ab+");
		fclose(f);
		f = fopen(fname, "rb+");

		// file failed to open
		if(!f) error("%s", "Error opening file for writing\n");

		// File is open, we now have to do the following:

		// 1. reset the bit array to 1 ("should write")
		memset(bitarray, 1, filecounts[i]);
		status = 1;
		hashctr = 0;

		// 2. set up a loop to read in all the lines of the file
		while(fread(&tmp_len, LEN_SIZE, 1, f)){
			// 2a. read in the node name
			memset(cached[hashctr], 0, MAX_NODE_NAME_SIZE);
			safe_fread(cached[hashctr], 1, tmp_len, f);
			safe_fread(&foundctr, L_SIZE, 1, f);
			hashctr++;

			// 2b. once we've read in enough lines...
			if(hashctr == FILE_READ_CACHE_SIZE){
				// ...hash all the names, then check if any of the to_read names match
				// if so, mark it in the bitarray array
				status = check_inputs_against_hashes(tmp_charptr, cached, bitarray, filecounts[i], hashctr);
				hashctr = 0;
			}
			if(!status) break;
		}
		// once we're at the end of the file, we have to do 2b again
		if(hashctr) status = check_inputs_against_hashes(tmp_charptr, cached, bitarray, filecounts[i], hashctr);
		if(status){
			// 3. Now we skip to the end of the file and write any line that still has bitarray == 1
			fseek(f, 0, SEEK_END);
			for(int j=0; j<filecounts[i]; j++){
				if(!bitarray[j]) continue;
				cur_len = strlen(tmp_charptr[j]);
				fwrite(&cur_len, LEN_SIZE, 1, f);
				fwrite(tmp_charptr[j], 1, cur_len, f);
				fwrite(&ctr, L_SIZE, 1, f);
				ctr++;
			}
		}
		// 4. close the file
		fclose(f);
	}

	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) free(cached[i]);

	return ctr;
}

int find_node_indices_batch(char **input_strings, full_read_line *file_lines, l_uint *indices, uint ninput, uint nfile){
	// input_strings is sorted (length, strcmp), file_strings is not -- first we sort file_strings quickly
	qsort(file_lines, nfile, sizeof(full_read_line), struct_frl_cmpfunc);

	// then we can find all matches in linear time
	uint i=0, j=0, previ=0, prevj=0;
	uint len1=strlen(input_strings[i]), len2=file_lines[j].len;
	int cmp;

	while(i < ninput && j < nfile){
		if(i != previ){
			len1 = strlen(input_strings[i]);
			previ = i;
		}
		if(j != prevj){
			len2 = file_lines[j].len;
			prevj = j;
		}

		// skip past anything we've already found
		if(indices[i] || len1 < len2){
			i++;
		} else if(len1 > len2){
			j++;
		} else {
			cmp = strcmp(input_strings[i], file_lines[j].name);
			if(cmp < 0){
				i++;
			} else if (cmp > 0){
				j++;
			} else {
				indices[i] = file_lines[j].ctr;
				i++;
				j++;
			}
		}
	}

	// everything already seen in this case means all indices non-zero
	for(int i=0; i<ninput; i++)
		if(!indices[i]) return 1;
	return 0;
}

void lookup_indices_batch(char** names, uint num_to_lookup, const uint num_unique_hashes, const char* dir,
													uint* hashes, uint* filecounts, l_uint* lookup_indices){
	const int LEN_SIZE = sizeof(uint16_t);

	// holds lines from the file
	FILE *f;
	int status;
	uint hashctr;
	uint16_t len;
	full_read_line cached[FILE_READ_CACHE_SIZE];
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++)
		cached[i].name = malloc(MAX_NODE_NAME_SIZE);
	char fname[PATH_MAX];
	char **tmp_charptr = names;
	l_uint *tmp_indices = lookup_indices;

	// now the new length of the array is insert_point
	// note that these are first sorted by hash, so once the hash value changes it will never appear again
	for(int i=0; i<num_unique_hashes; i++){
		// advance string and counts array to lines for this file
		if(i){
			tmp_charptr = &(tmp_charptr[filecounts[i-1]]);
			tmp_indices = &(tmp_indices[filecounts[i-1]]);
		}

		// build filename using the hash and open (should always exist)
		snprintf(fname, strlen(dir) + 6, "%s/%04x", dir, hashes[i]);
		f = fopen(fname, "rb");

		// file failed to open
		if(!f) error("%s", "Error opening file for writing\n");

		// File is open, we now have to do the following:

		// 1. reset variables
		status = 1;
		hashctr = 0;

		// set up a loop to read in all the lines of the file
		while(fread(&len, LEN_SIZE, 1, f)){
			// 2a. read in the node name
			cached[hashctr].len = len;
			memset(cached[hashctr].name, 0, MAX_NODE_NAME_SIZE);
			safe_fread(cached[hashctr].name, 1, len, f);
			safe_fread(&(cached[hashctr].ctr), L_SIZE, 1, f);
			hashctr++;

			// 2b. once we've read in enough lines...
			if(hashctr == FILE_READ_CACHE_SIZE){
				// ...hash all the names, then check if any of the to_read names match
				// if so, mark it in the bitarray array
				status = find_node_indices_batch(tmp_charptr, cached, tmp_indices, filecounts[i], hashctr);
				hashctr = 0;
			}
			if(!status) break;
		}

		// once we're at the end of the file, we have to do 2b again
		if(hashctr) status = find_node_indices_batch(tmp_charptr, cached, tmp_indices, filecounts[i], hashctr);
		fclose(f);
	}
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) free(cached[i].name);

	return;
}

l_uint batch_write_edgecounts(char **names, int num_to_sort, const char *dir, l_uint ctr, FILE *tab){

	// unique hash values, number of entries per hash counter, counts for each entry
	uint hashes[num_to_sort], filecounts[num_to_sort];
	uint counts[num_to_sort];
	l_uint *lookup_indices = NULL;

	// these are sometimes used for temp storage
	uint num_unique_hashes;
	l_uint foundctr;
	int insert_point;

	unique_strings_with_sideeffects(names, num_to_sort, filecounts, hashes, &num_unique_hashes, &insert_point, counts, 1);
	lookup_indices = calloc(insert_point, L_SIZE);
	lookup_indices_batch(names, insert_point, num_unique_hashes, dir, hashes, filecounts, lookup_indices);

	// 3. Now we write all our counts using lookup_indices and counts
	long offset = 0;
	rewind(tab);
	for(int i=0; i<insert_point; i++){
		offset = i ? (lookup_indices[i] - lookup_indices[i-1] - 1) : (lookup_indices[i]-1);
		fseek(tab, offset*L_SIZE, SEEK_CUR);
		safe_fread(&foundctr, L_SIZE, 1, tab);
		fseek(tab, -1*L_SIZE, SEEK_CUR);
		foundctr += counts[i];
		fwrite(&foundctr, L_SIZE, 1, tab);
	}
	free(lookup_indices);
	return ctr;
}

l_uint hash_file_vnames_batch(const char* fname, const char* dname, const char *ftable,
	const char sep, const char line_sep, l_uint ctr, int v, int is_undirected, int count_edges){
	/*
	 * fname: .tsv list of edges
	 * dname: directory of hash codes
	 * ftable: output counts file, should be L_SIZE*num_v bytes
	 */
	FILE *f = fopen(fname, "rb");
	FILE *tab = NULL;
	if(count_edges) tab = fopen(ftable, "rb+");
	char *namecache[NODE_NAME_CACHE_SIZE];
	int itermax = (count_edges && !is_undirected) ? 1 : 2;

	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) namecache[i] = malloc(MAX_NODE_NAME_SIZE);

	char vname[MAX_NODE_NAME_SIZE];
	int cur_pos = 0, cachectr=0;
	char c = getc(f);
	l_uint cur_ctr = ctr;
	l_uint print_counter = 0;

	if(v) Rprintf("Reading file %s...\n", fname);

	while(!feof(f)){
		// going to assume we're at the beginning of a line
		// lines should be of the form `start end weight`
		// separation is by char `sep`
		for(int iter=0; iter<itermax; iter++){
			cur_pos = 0;
			while(c != sep){
				vname[cur_pos++] = c;
				c = getc(f);
				if(cur_pos >= MAX_NODE_NAME_SIZE)
					errorclose_file(f, tab, "Node name is larger than max allowed name size.\n");

				if(feof(f)) errorclose_file(f, tab, "Unexpected end of file.\n");
			}

			vname[cur_pos] = '\0';
			memset(namecache[cachectr], 0, MAX_NODE_NAME_SIZE);
			memcpy(namecache[cachectr], vname, strlen(vname));
			cachectr++;
			if(cachectr == NODE_NAME_CACHE_SIZE){
				if(count_edges)
					cur_ctr = batch_write_edgecounts(namecache, cachectr, dname, cur_ctr, tab);
				else
					cur_ctr = batch_write_nodes(namecache, cachectr, dname, cur_ctr);
				cachectr = 0;
			}

			c = getc(f);
		}

		while(c != line_sep && !feof(f)) c = getc(f);
		if(c == line_sep) c=getc(f);
		print_counter++;
		if(print_counter % PRINT_COUNTER_MOD == 0){
			if(v) Rprintf("\t%lu lines read\r", print_counter);
			else R_CheckUserInterrupt();
			R_CheckUserInterrupt();
		}
	}
	if(cachectr){
		if(count_edges)
			cur_ctr = batch_write_edgecounts(namecache, cachectr, dname, cur_ctr, tab);
		else
			cur_ctr = batch_write_nodes(namecache, cachectr, dname, cur_ctr);
	}

	if(v){
		Rprintf("\t%lu lines read\n", print_counter);
		if(!count_edges)
			Rprintf("\t%lu total nodes.\n", cur_ctr-1);
	}
	fclose(f);
	if(count_edges) fclose(tab);

	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) free(namecache[i]);
	return cur_ctr-1;
}

l_uint write_counts_batch(FILE *ftabptr, l_uint* indices, uint* counts){
	l_uint cur_count, total_count=0;
	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++){
		// once we reach a 0, we're done
		if(indices[i] == 0) return total_count;
		fseek(ftabptr, (indices[i]-1)*L_SIZE, SEEK_SET);
		safe_fread(&cur_count, L_SIZE, 1, ftabptr);
		fseek(ftabptr, -1*L_SIZE, SEEK_CUR);
		cur_count += counts[i];
		total_count += counts[i];
		fwrite(&cur_count, L_SIZE, 1, ftabptr);
	}
	return total_count;
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
		if(print_counter % PRINT_COUNTER_MOD == 0){
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

void reformat_counts(const char* curcounts, const char* mastertable, l_uint n_vert, int add_self_loops){
	/*
	 * Creates a new table with cumulative counts
	 * leaves the old table unchanged, this will act as a temporary counts file later
	 */
	const uint l_size = L_SIZE;
	l_uint cumul_total = 0, curcount;
	FILE *tmptab = fopen(curcounts, "rb");
	FILE *mtab = fopen(mastertable, "wb+");
	int self_loop = add_self_loops ? 1 : 0;

	for(l_uint i=0; i<n_vert; i++){
		fwrite(&cumul_total, l_size, 1, mtab);
		safe_fread(&curcount, l_size, 1, tmptab);
		cumul_total += curcount + self_loop; // add an extra count for each node if we add self loops
	}

	// ending position of file
	fwrite(&cumul_total, l_size, 1, mtab);
	fclose(tmptab);
	fclose(mtab);
	return;
}

void add_self_loops_to_csrfile(const char *ftable, l_uint num_v){
	// If self loops are included, the first entry for each node remains empty
	// here we'll fill it in with the node itself
	// thus, we can just write to whatever the first index of the value is and set the value to 0
	const uint entry_size = L_SIZE + sizeof(double);
	const double self_weight = 0.5;

	l_uint tmp_pos;

	FILE *mastertab = fopen(ftable, "rb+");
	if(!mastertab) error("%s", "error opening CSR file.\n");

	for(l_uint i=0; i<num_v; i++){
		fseek(mastertab, i*L_SIZE, SEEK_SET);
		safe_fread(&tmp_pos, L_SIZE, 1, mastertab);
		// now we're at position i+1, need to go to num_v+2
		// so we just have to move forward (num_v-i+1) positions to get to the end
		fseek(mastertab, (num_v-i)*L_SIZE, SEEK_CUR);

		// then move forward another [entry] amounts and write the current index to get a self loop
		fseek(mastertab, tmp_pos*entry_size, SEEK_CUR);
		fwrite(&i, L_SIZE, 1, mastertab);
		fwrite(&self_weight, sizeof(double), 1, mastertab);
	}

	fclose(mastertab);
}

void normalize_csr_edgecounts(const char* ftable, l_uint num_v){
	const int entry_size = L_SIZE + sizeof(double);
	double tmp_val, normalizer;
	l_uint start, end;
	FILE *mastertab = fopen(ftable, "rb+");
	if(!mastertab) error("%s", "error opening CSR file.\n");

	start = 0;
	for(l_uint i=1; i<num_v; i++){
		normalizer = 0;
		fseek(mastertab, i*L_SIZE, SEEK_SET);
		safe_fread(&end, L_SIZE, 1, mastertab);

		// moving one extra L_SIZE so we're offset to the weights only
		fseek(mastertab, (num_v-i+2)*L_SIZE, SEEK_CUR);
		fseek(mastertab, entry_size*start, SEEK_CUR);
		for(l_uint j=0; j<(end-start); j++){
			safe_fread(&tmp_val, sizeof(double), 1, mastertab);
			normalizer += tmp_val;
			fseek(mastertab, L_SIZE, SEEK_CUR);
		}

		// now we're at the double entry at (end + 1), need to move back to start
		// that means moving back (end+1-start) spaces
		fseek(mastertab, -1*entry_size*(end-start+1), SEEK_CUR);
		normalizer = normalizer == 0 ? 1 : normalizer;

		// finally we overwrite each of the values
		for(l_uint j=0; j<(end-start); j++){
			safe_fread(&tmp_val, sizeof(double), 1, mastertab);
			tmp_val /= normalizer;
			fseek(mastertab, -1*sizeof(double), SEEK_CUR);
			fwrite(&tmp_val, sizeof(double), 1, mastertab);
			fseek(mastertab, L_SIZE, SEEK_CUR);
		}
		start = end;
	}

	fclose(mastertab);
	return;
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
		if(print_counter % PRINT_COUNTER_MOD == 0){
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


void batch_write_edgelocs(char** names_cache, double* weights_cache, int cache_size,
													const char* dname, FILE *countstab, FILE *mastertab,
													size_t entrysize, l_uint num_v,
													const int is_undirected, const int self_loop_inc){
	l_uint *indices, *lookup_indices;
	// temporary copy of the names cache, since the sort algorithm is destructive
	char *names_cache_copy[NODE_NAME_CACHE_SIZE];
	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++){
		names_cache_copy[i] = malloc(MAX_NODE_NAME_SIZE);
		memcpy(names_cache_copy[i], names_cache[i], MAX_NODE_NAME_SIZE);
	}

	// sort and unique the names we found
	uint num_unique_hashes;
	int num_unique_names;
	uint hashes[cache_size], filecounts[cache_size];
	unique_strings_with_sideeffects(names_cache_copy, cache_size, filecounts, hashes,
																	&num_unique_hashes, &num_unique_names, NULL, 0);

	// batch lookup all indices, store in lookup_indices
	lookup_indices = calloc(num_unique_names, L_SIZE);
	lookup_indices_batch(names_cache_copy, num_unique_names, num_unique_hashes, dname, hashes, filecounts, lookup_indices);

	indices = calloc(cache_size, L_SIZE);
	for(int i=0; i<cache_size; i++){
		for(int j=0; j<num_unique_names; j++){
			if(strcmp(names_cache[i], names_cache_copy[j]) == 0){
				indices[i] = lookup_indices[j];
				break;
			}
		}
	}

	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) free(names_cache_copy[i]);
	free(lookup_indices);

	int itermax = is_undirected ? 2 : 1;
	l_uint inds[2], offset, loc;
	double w;
	for(int i=0; i<cache_size/2; i++){
		w = weights_cache[i];
		inds[0] = indices[i*2]-1;
		inds[1] = indices[i*2+1]-1;
		// run twice if undirected to get both directions
		for(int j=0; j<itermax; j++){
			//if(!j) Rprintf("%lu - %lu (%0.2f)\n", inds[0], inds[1], w);
			// get offset for location we'll write to in the counts file
			fseek(countstab, (inds[j])*L_SIZE, SEEK_SET);
			safe_fread(&offset, L_SIZE, 1, countstab);

			// get start location in table file we'll write to
			fseek(mastertab, (inds[j])*L_SIZE, SEEK_SET);
			safe_fread(&loc, L_SIZE, 1, mastertab);

			// write the edge
			fseek(mastertab, (num_v+1)*L_SIZE, SEEK_SET);
			fseek(mastertab, (loc+offset+self_loop_inc)*entrysize, SEEK_CUR);
			fwrite(&inds[(j+1)%2], L_SIZE, 1, mastertab);
			fwrite(&w, sizeof(double), 1, mastertab);

			// decrement the counts file
			fseek(countstab, inds[j]*L_SIZE, SEEK_SET);
			offset--;
			fwrite(&offset, L_SIZE, 1, countstab);
		}
	}
	free(indices);
	return;
}

void csr_compress_edgelist_batch(const char* edgefile, const char* dname, const char* curcountfile, const char* ftable,
														const char sep, const char linesep, l_uint num_v, int v, const int is_undirected, int has_self_loops){
	/*
	 * This should be called after we've already read in all our files
	 * critically, ensure we're rewritten our ftable file such that it is cumulative counts and not vertex counts
	 *
	 * Error checking can be reduced because we would have caught it earlier
	 * cache is going to store v1 v2 at pos i, i+1; weight stored at i/2
	 */
	const size_t entry_size = L_SIZE + sizeof(double);
	const int self_loop_inc = has_self_loops ? 1 : 0;
	char *names_cache[NODE_NAME_CACHE_SIZE];
	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) names_cache[i] = malloc(MAX_NODE_NAME_SIZE);
	double weights_cache[NODE_NAME_CACHE_SIZE/2];

	char vname[MAX_NODE_NAME_SIZE];
	int cachectr=0, stringctr=0;
	l_uint print_counter = 0;

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

	if(v) Rprintf("Reading edges from file %s...\n", edgefile);

	char c = getc(edgelist);
	while(!feof(edgelist)){
		// read in the two vertex names
		for(int i=0; i<2; i++){
			stringctr = 0;
			memset(vname, 0, MAX_NODE_NAME_SIZE);
			while(c != sep){
				vname[stringctr++] = c;
				c = getc(edgelist);
			}
			memset(names_cache[cachectr], 0, MAX_NODE_NAME_SIZE);
			memcpy(names_cache[cachectr], vname, strlen(vname));
			cachectr++;

			// advance one past the separator
			c = getc(edgelist);
		}
		// read in the weight
		stringctr = 0;
		memset(vname, 0, MAX_NODE_NAME_SIZE);
		c = getc(edgelist);
		while(c != linesep){
			vname[stringctr++] = c;
			c = getc(edgelist);
		}
		// advance one past the separator
		c = getc(edgelist);

		weights_cache[(cachectr/2)-1] = atof(vname);
		print_counter++;
		if(print_counter % PRINT_COUNTER_MOD == 0){
			if(v) Rprintf("\t%lu edges read\r", print_counter);
			else R_CheckUserInterrupt();
		}

		if((cachectr + 1) >= FILE_READ_CACHE_SIZE){
			// if size is odd, we need to stop early (only using every other index)
			batch_write_edgelocs(names_cache, weights_cache, cachectr, dname, tmptable, mastertable,
													entry_size, num_v, is_undirected, self_loop_inc);
			cachectr = 0;
		}

	}
	if(cachectr){
		batch_write_edgelocs(names_cache, weights_cache, cachectr, dname, tmptable, mastertable,
												entry_size, num_v, is_undirected, self_loop_inc);
	}

	if(v) Rprintf("\t%lu edges read\n", print_counter);
	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) free(names_cache[i]);
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
	safe_fread(&start, L_SIZE, 1, mastertab);
	safe_fread(&end, L_SIZE, 1, mastertab);

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
		if(use_limited_nodes && unif_rand() > acceptance_prob){
			fseek(mastertab, L_SIZE+sizeof(double), SEEK_CUR);
			continue;
		}

		// read in neighbor and weight
		safe_fread(&tmp_id, L_SIZE, 1, mastertab);
		safe_fread(&tmp_w, sizeof(double), 1, mastertab);

		// get which cluster it belongs to
		fseek(clusterings, L_SIZE*tmp_id, SEEK_SET);
		safe_fread(&tmp_cl, L_SIZE, 1, clusterings);

		R_CheckUserInterrupt();
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

void add_to_queue(l_uint clust, l_uint ind, l_uint n_node, FILE *clust_f, FILE *master_f, FILE *q_f, FILE *ctrq_f){
	l_uint start, end, tmp_ind, tmp_cl, nedge;
	l_uint buf[MAX_EDGES_EXACT];
	double dummy;
	int ctr = 0, found;

	fseek(master_f, ind*L_SIZE, SEEK_SET);
	safe_fread(&start, L_SIZE, 1, master_f);
	safe_fread(&end, L_SIZE, 1, master_f);
	fseek(master_f, (n_node+1)*L_SIZE, SEEK_SET);
	fseek(master_f, start*(L_SIZE+sizeof(double)), SEEK_CUR);

	nedge = end - start;
	for(l_uint i=0; i<nedge; i++){
		safe_fread(&tmp_ind, L_SIZE, 1, master_f);
		safe_fread(&dummy, sizeof(double), 1, master_f);
		fseek(clust_f, L_SIZE*tmp_ind, SEEK_SET);
		safe_fread(&tmp_cl, L_SIZE, 1, clust_f);

		if(tmp_cl && tmp_cl == clust) continue;
		tmp_ind++;
		found = 0;
		for(int j=0; j<ctr; j++){
			if(buf[j] == tmp_ind){
				found = 1;
				break;
			}
		}
		if(!found){
			buf[ctr++] = tmp_ind;
			if(ctr == MAX_EDGES_EXACT) break;
		}
	}

	// iterate over queue file, adding numbers if not already there
	rewind(q_f);
	for(int j=0; j<ctr; j++){
		fseek(ctrq_f, buf[j], SEEK_SET);
		found = getc(ctrq_f);
		if(!found){
			fseek(ctrq_f, -1, SEEK_CUR);
			putc(1, ctrq_f);
		} else {
			buf[j] = 0;
		}
	}

	// this is just in case, it adds a little runtime but it's safer to guard fread errors
	fseek(q_f, 0, SEEK_END);
	for(int j=0; j<ctr; j++){
		if(buf[j]){
			buf[j]--;
			fwrite(&buf[j], L_SIZE, 1, q_f);
		}
	}

	return;
}

l_uint get_qsize(FILE *q){
	l_uint scratch, ctr=0;
	while(fread(&scratch, L_SIZE, 1, q)) ctr++;
	rewind(q);
	return ctr;
}

void initialize_queue(FILE *q, l_uint maxv, FILE *ctr_file){
	GetRNGstate();
	l_uint j, tmp;
	for(l_uint i=0; i<maxv; i++){
		putc(1, ctr_file);
		j = (l_uint) trunc((i+1) * (unif_rand()));
		if(j < i){
			// guarding edge case where unif_rand() returns 1.0

			// tmp = arr[j]
			fseek(q, L_SIZE*j, SEEK_SET);
			safe_fread(&tmp, L_SIZE, 1, q);

			// arr[j] = i
			fseek(q, -1*L_SIZE, SEEK_CUR);
			fwrite(&i, L_SIZE, 1, q);
		} else {
			tmp = i;
		}

		// arr[i] = tmp
		fseek(q, L_SIZE*i, SEEK_SET);
		fwrite(&tmp, L_SIZE, 1, q);
	}
	PutRNGstate();

	return;
}

void cluster_file(const char* mastertab_fname, const char* clust_fname,
									const char *qfile_f1, const char *qfile_f2, const char *qfile_log,
									l_uint num_v, int max_iterations, int v){
	// main runner function to cluster nodes
	FILE *masterfile = fopen(mastertab_fname, "rb");
	FILE *clusterfile = fopen(clust_fname, "rb+");
	FILE *cur_q, *next_q, *ctr_q;
	const char* queues[] = {qfile_f1, qfile_f2};
	const char* progress = "\\|/-\\|/-";

	ctr_q = fopen(qfile_log, "wb+");
	l_uint cluster_res, qsize, tmp_ind;

	// randomly initialize queue and ctr file
	if(v) Rprintf("Clustering: |\r");
	cur_q = fopen(queues[0], "wb+");
	initialize_queue(cur_q, num_v, ctr_q);
	fclose(cur_q);

	for(int i=0; i<max_iterations; i++){
		cur_q = fopen(queues[i%2], "rb+");
		next_q = fopen(queues[(i+1)%2], "wb+");
		qsize = get_qsize(cur_q);
		if(qsize == 0) break;

		while(fread(&tmp_ind, L_SIZE, 1, cur_q)){
			fseek(ctr_q, tmp_ind, SEEK_SET);
			putc(0, ctr_q);
			cluster_res = update_node_cluster(tmp_ind, num_v+1, masterfile, clusterfile);
			add_to_queue(cluster_res, tmp_ind, num_v, clusterfile, masterfile, next_q, ctr_q);
		}

		fclose(cur_q);
		fclose(next_q);
		if(v) Rprintf("Clustering: %.0f%% complete %c\r", ((double)(i+1) / max_iterations)*100, progress[i%8]);
	}
	if(v) Rprintf("Clustering: 100%% complete.   \n");
	fclose(ctr_q);
	fclose(masterfile);

	reformat_clusters(clusterfile, num_v);
	fclose(clusterfile);
	remove(queues[0]);
	remove(queues[1]);
	remove(qfile_log);
	return;
}

SEXP R_hashedgelist(SEXP FILENAME, SEXP NUM_EFILES, SEXP TABNAME, SEXP TEMPTABNAME, SEXP QFILES, SEXP OUTDIR,
										SEXP SEPS, SEXP CTR, SEXP ITER, SEXP VERBOSE, SEXP IS_UNDIRECTED, SEXP ADD_SELF_LOOPS){
	/*
	 * I always forget how to handle R strings so I'm going to record it here
	 * R character vectors are STRSXPs, which is the same as a list (VECSXP)
	 * Each entry in a STRSXP is a CHARSXP, accessed with STRING_ELT(str, ind)
	 * You can get a `const char*` from a CHARSXP with CHAR()
	 * You can also re-encode with `const char* Rf_translateCharUTF8()`
	 */

	/*
	 * Input explanation:
	 *		 FILENAME: file of edges in format `v1 v2 w`
	 *      TABNAME: stores CSR compression of graph structure
	 *	TEMPTABNAME: first used to count edges, then used to store clusters
	 * 			 QFILES: two files, both used for queues
	 * 		   OUTDIR: directory to store hashed strings
	 *
	 * R_hashedgelist(tsv, csr, clusters, queues, hashdir, seps, 1, iter, verbose)
	 */

	const char* dir = CHAR(STRING_ELT(OUTDIR, 0));
	const char* edgefile;
	const char* tabfile = CHAR(STRING_ELT(TABNAME, 0));
	const char* temptabfile = CHAR(STRING_ELT(TEMPTABNAME, 0));
	const char* seps = CHAR(STRING_ELT(SEPS, 0));
	const char* qfile1 = CHAR(STRING_ELT(QFILES, 0));
	const char* qfile2 = CHAR(STRING_ELT(QFILES, 1));
	const char* qfile3 = CHAR(STRING_ELT(QFILES, 2));
	const int num_edgefiles = INTEGER(NUM_EFILES)[0];
	const int num_iter = INTEGER(ITER)[0];
	const int verbose = LOGICAL(VERBOSE)[0];
	const int is_undirected = LOGICAL(IS_UNDIRECTED)[0];
	l_uint num_v = (l_uint)(REAL(CTR)[0]);
	const int add_self_loops = LOGICAL(ADD_SELF_LOOPS)[0];

	/* Old method
	// first, index all vertex names and record how many edges each has
	for(int i=0; i<num_edgefiles; i++){
		edgefile = CHAR(STRING_ELT(FILENAME, i));
		num_v += hash_file_vnames(edgefile, dir, temptabfile, seps[0], seps[1], num_v, verbose, is_undirected);
	}
 	// num_v will always have the *next* location that we would insert a node at
 	// thus, once we're all done we need to decrement it by 1.
 	num_v--;
 	*/

 	/* New method: uses caching to speed up vertex name reading */
 	// first, index all vertex names
	for(int i=0; i<num_edgefiles; i++){
		edgefile = CHAR(STRING_ELT(FILENAME, i));
		num_v += hash_file_vnames_batch(edgefile, dir, NULL, seps[0], seps[1], num_v, verbose, is_undirected, 0);
	}
	// num_v will always have the *next* location that we would insert a node at
 	// thus, once we're all done we need to decrement it by 1.
	num_v--;

	// set up the table file with 0s for all vertices
	// have to while loop this, not guaranteed to write this many values in one call
	char buf[4096];
	memset(buf, 0, 4096);
	l_uint to_write = num_v * L_SIZE;
	FILE *f = fopen(temptabfile, "wb");
	while(to_write > 0){
		to_write -= fwrite(buf, 1, to_write > 4096 ? 4096 : to_write, f);
	}
	fclose(f);

	// next, record how many edges each has
	for(int i=0; i<num_edgefiles; i++){
		edgefile = CHAR(STRING_ELT(FILENAME, i));
		hash_file_vnames_batch(edgefile, dir, temptabfile, seps[0], seps[1], 1, verbose, is_undirected, 1);
	}

 	// next, create an indexed table file of where data for each vertex will be located
 	if(verbose) Rprintf("Reformatting counts file...\n");
 	reformat_counts(temptabfile, tabfile, num_v, add_self_loops);

 	// then, we'll create the CSR compression of all our edges
 	for(int i=0; i<num_edgefiles; i++){
 		edgefile = CHAR(STRING_ELT(FILENAME, i));
 		csr_compress_edgelist_batch(edgefile, dir, temptabfile, tabfile, seps[0], seps[1], num_v, verbose, is_undirected, add_self_loops);
 	}

 	if(add_self_loops && verbose) Rprintf("Adding self loops...\n");
 	if(add_self_loops) add_self_loops_to_csrfile(tabfile, num_v);

 	//if(verbose) Rprintf("Normalizing node edge weights...\n");
 	//normalize_csr_edgecounts(tabfile, num_v);

 	// temptabfile now becomes our clustering file
 	cluster_file(tabfile, temptabfile, qfile1, qfile2, qfile3, num_v, num_iter, verbose);

	SEXP RETVAL = PROTECT(allocVector(REALSXP, 1));
	REAL(RETVAL)[0] = (double) num_v;
	UNPROTECT(1);
	return RETVAL;
}

SEXP R_write_output_clusters(SEXP CLUSTERFILE, SEXP HASHEDFILES, SEXP NUM_FILES, SEXP OUTFILE, SEXP SEPS){
	// we're going to leverage list.files() on the R side for this so I don't have to worry
	// about platform-specific file handling for directory parsing

	/*
	 * Inputs:
	 * 	CLUSTERFILE: file of clusters (cluster_counts)
	 *  HASHEDFILES: character vector of all hash files
	 *    NUM_FILES: number of hash files
	 *      OUTFILE: file to write to
	 *         SEPS: %c%c, entry_sep, line_sep
	 *
	 * R_write_output_clusters(cluster_counts, list.files(hashdir), length(...), out_tsvpath, seps)
	 */

	const int nfiles = INTEGER(NUM_FILES)[0];
	const char* cfile = CHAR(STRING_ELT(CLUSTERFILE, 0));
	const char* outfile = CHAR(STRING_ELT(OUTFILE, 0));
	const char* seps = CHAR(STRING_ELT(SEPS, 0));
	const char* fname;
	uint16_t name_len;
	int LEN_SIZE = sizeof(uint16_t);
	char buf[MAX_NODE_NAME_SIZE];
	char write_buf[PATH_MAX];
	l_uint index, clust;
	FILE *outf = fopen(outfile, "w");
	FILE *f;

	// remember that indices will be off by one
	FILE *fclusters = fopen(cfile, "rb");
	for(int i=0; i<nfiles; i++){
		fname = CHAR(STRING_ELT(HASHEDFILES,i));
		f = fopen(fname, "rb");
		while(fread(&name_len, LEN_SIZE, 1, f)){
			// read in data from node names
			memset(buf, '\0', MAX_NODE_NAME_SIZE);
			memset(write_buf, '\0', PATH_MAX);
			safe_fread(buf, 1, name_len, f);
			safe_fread(&index, L_SIZE, 1, f);

			// get the cluster for the node name
			fseek(fclusters, L_SIZE*(index-1), SEEK_SET);
			safe_fread(&clust, L_SIZE, 1, fclusters);

			// prepare data for output
			snprintf(write_buf, (name_len+3)+L_SIZE, "%s%c%llu%c", buf, seps[0], clust, seps[1]);
			fwrite(write_buf, 1, strlen(write_buf), outf);
		}
		fclose(f);
	}

	fclose(fclusters);
	fclose(outf);

	return R_NilValue;
}