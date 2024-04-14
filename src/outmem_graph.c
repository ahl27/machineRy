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
#define FILE_READ_CACHE_SIZE 4096
#define NUM_BITS_HASH 4 // number of bits for hash file usage, max 32

const char DELIM = 23; // 23 = end of transmission block, not really used for anything nowadays
const int L_SIZE = sizeof(l_uint);
const int MAX_READ_RETRIES = 10;

// set this to 1 if we should sample edges rather than use all of them
const int use_limited_nodes = 0;
const l_uint MAX_EDGES_EXACT = 20000; // this is a soft cap -- if above this, we sample edges probabalistically
const int PRINT_COUNTER_MOD = 170;

typedef struct {
	uint len;
	char *name;
	l_uint ctr;
} full_read_line;

typedef struct {
	l_uint ctr1;
	l_uint ctr2;
} double_lu;

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

static void safe_filepath_cat(const char* dir, const char* f, char *fname, size_t fnamesize){
	// fname should be preallocated
	char directory_separator;
#ifdef WIN32
	directory_separator = '\\';
#else
	directory_separator = '/';
#endif
	memset(fname, 0, fnamesize);
	// length is len(dir) + len(f) + 1 (separator) + 1 (terminator)
	snprintf(fname, strlen(dir)+strlen(f)+2, "%s%c%s", dir, directory_separator, f);
	return;
}

int l_uint_compar(const void* a, const void* b){
	double_lu aa = **(double_lu **)(a);
	double_lu bb = **(double_lu **)(b);
	if(aa.ctr2 - bb.ctr2)
		return aa.ctr2 - bb.ctr2;
	return aa.ctr1 - bb.ctr1;
}

void precopy_dlu1(const char* f1, const char* f2){
	// write and add index
	double_lu dlu = {1,0};
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");
	while(fread(&dlu.ctr2, L_SIZE, 1, orig)){
		fwrite(&dlu, sizeof(double_lu), 1, copy);
		dlu.ctr1++;
	}
	fclose(orig);
	fclose(copy);
	return;
}

void precopy_dlu2(const char* f1, const char* f2){
	// write flipped version (index, clust => clust, index)
	double_lu dlu = {0,0};
	l_uint prev_ind=0;
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");
	while(fread(&dlu, sizeof(double_lu), 1, orig)){
			prev_ind = dlu.ctr1;
			dlu.ctr1 = dlu.ctr2;
			dlu.ctr2 = prev_ind;
			fwrite(&dlu, sizeof(double_lu), 1, copy);
	}
	fclose(orig);
	fclose(copy);
	return;
}

void postcopy_dlu1(const char* f1, const char* f2){
	double_lu dlu = {0,0};
	l_uint prev_ind=0, ctr=0;
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");
	while(fread(&dlu, sizeof(double_lu), 1, orig)){
			if(prev_ind != dlu.ctr2){
				prev_ind = dlu.ctr2;
				dlu.ctr2 = ++ctr;
			} else {
				dlu.ctr2 = ctr;
			}
			fwrite(&dlu, sizeof(double_lu), 1, copy);
	}
	fclose(orig);
	fclose(copy);
	return;
}

void postcopy_dlu2(const char* f1, const char* f2){
	// write only the cluster into file
	// also can reindex such that the first cluster listed is cluster 1
	// (this can be removed)
	double_lu dlu = {0,0};
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");

	// uncomment these lines to make the first vertex have cluster 1
	/*
	l_uint max_found = 0, offset;
	while(fread(&dlu, sizeof(double_lu), 1, orig))
			if(dlu.ctr1 > max_found) max_found = dlu.ctr1;
	rewind(orig);
	fread(&dlu, sizeof(double_lu), 1, orig);
	offset = max_found - dlu.ctr1;
	rewind(orig);
	*/

	while(fread(&dlu, sizeof(double_lu), 1, orig)){
			// dlu.ctr1 = ((dlu.ctr1 + offset) % max_found) + 1;
			fwrite(&dlu.ctr1, L_SIZE, 1, copy);
	}
	fclose(orig);
	fclose(copy);
	return;
}



void mergesort_clust_file(const char* f, const char* dir, size_t element_size,
															int (*compar)(const void *, const void *),
															void (*precopy)(const char*, const char*),
															void (*postcopy)(const char*, const char*)){
	/*
	 * general file mergesort function
	 * arguments:
	 *	-            f: file to sort
	 *	-          dir: directory to store junk files
	 *  - element_size: size of each element to read/write
	 *  -      *compar: function pointer used in qsort / mergesort. ensure proper casting.
	 *  -     *precopy: function to copy f into the first junk file
	 *  -    *postcopy: function to write final values back into f
	 *  notes:
	 *  - *precopy should open the file (assume it does not exist)
	 *  - *compar will provide void** values, make sure to double dereference
	 */

	// two read pointers, one write pointer
	FILE *f1_r1, *f1_r2, *f2_w;
	char file1[PATH_MAX], file2[PATH_MAX];
	char *finalfile;
	//size_t dlu_size = sizeof(double_lu);

	// create the junk files we'll use
	safe_filepath_cat(dir, "tmp_ms1", file1, PATH_MAX);
	safe_filepath_cat(dir, "tmp_ms2", file2, PATH_MAX);

	// first, we'll use the cache to read in preprocessed sorted blocks of size `element_size`
	l_uint block_size = FILE_READ_CACHE_SIZE;
	l_uint total_lines = 0;

	// allocate space for data in a void*
	void *read_cache[FILE_READ_CACHE_SIZE];
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) read_cache[i] = malloc(element_size);

	// copy the original file into file 1
	precopy(f, file1);

	// open file, read in chunks, sort locally, write to file
	uint cachectr = 0;
	f1_r1 = fopen(file1, "rb");
	if(!f1_r1) error("%s", "Error opening file obtained from mergesort precopy");
	f2_w = fopen(file2, "wb");
	if(!f2_w) error("%s", "Error opening temporary mergesort file for writing");
	while(fread(read_cache[cachectr++], element_size, 1, f1_r1)){
		total_lines++;
		if(cachectr == block_size){
			qsort(read_cache, cachectr, sizeof(void*), compar);
			for(int i=0; i<cachectr; i++)
				fwrite(read_cache[i], element_size, 1, f2_w);
			cachectr=0;
		}
	}
	if(cachectr){
		qsort(read_cache, cachectr, sizeof(void*), compar);
		for(int i=0; i<cachectr; i++)
			fwrite(read_cache[i], element_size, 1, f2_w);
	}

	fclose(f1_r1);
	fclose(f2_w);
	finalfile = file2;

	l_uint cur_lines = 0;
	int iter1, iter2, previt1, previt2;
	int flip = 0, cmp;
	void *tmp1 = malloc(element_size);
	void *tmp2 = malloc(element_size);
	char *f1, *f2;
	while(block_size < total_lines){
		// f1 is always the reading file, f2 the writing file
		f1 = flip ? file1 : file2;
		f2 = flip ? file2 : file1;
		flip = !flip;

		f1_r1 = fopen(f1, "rb");
		f1_r2 = fopen(f1, "rb");
		f2_w = fopen(f2, "wb");
		if(!(f1_r1 && f1_r2 && f2_w))
			error("%s", "Error opening temporary files in mergesort");
		// move second pointer forward to second block
		fseek(f1_r2, element_size*block_size, SEEK_CUR);

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
				if(iter1 && iter1 != previt1){
					safe_fread(tmp1, element_size, 1, f1_r1);
					previt1 = iter1;
				}
				if(iter2 && iter2 != previt2){
					safe_fread(tmp2, element_size, 1, f1_r2);
					previt2 = iter2;
				}

				cmp = compar(&tmp1, &tmp2);
				if(iter1 && (!iter2 || cmp <= 0 )){
					fwrite(tmp1, element_size, 1, f2_w);
					iter1--;
				} else {
					fwrite(tmp2, element_size, 1, f2_w);
					iter2--;
				}
			}
			// advance pointers one block past where we just read:
			// if we move too far it doesn't really matter, we'll catch it on the next part
			fseek(f1_r1, element_size*block_size, SEEK_CUR);
			fseek(f1_r2, element_size*block_size, SEEK_CUR);
		}

		fclose(f1_r1);
		fclose(f1_r2);
		fclose(f2_w);
		cur_lines = 0;
		block_size *= 2;
		finalfile = f2;
	}
	// free memory allocations
	free(tmp1);
	free(tmp2);
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) free(read_cache[i]);

	// copy result back into f
	postcopy(finalfile, f);

	return;
}

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
	 * hashes to a 32-bit uint that we truncate depending on the value of NUM_BITS_HASH
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
	//hash = (hash & 0x0000FFFF) ^ ((hash & 0xFFFF0000) >> 16);

	/*
	 * take the lowest NUM_BITS_HASH bits -- the fewer files, the better due to batch processing
	 *
	 * The important consideration is that reading is fast, but writing and open/closing is slow
	 * Thus having a few big files is a lot faster than lots of small ones.
	 * However, having a few big files affects lookup time of rw_vertname -- likely want to batch/cache this as well
	 */
	hash &= ((1ULL << NUM_BITS_HASH) - 1);
	return hash;
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
	uint cur_len, tmp_hash;

	// first sort the array
	qsort(names, num_to_sort, sizeof(char*), node_name_cmpfunc);

	// next, unique the values
	hashes[0] = hash_string_fnv(names[0]);
	if(useCounts) counts[0] = 1;
	cur_len = strlen(names[0]);
	for(int i=1; i<num_to_sort; i++){
		if(cur_len != strlen(names[i]) || (strcmp(names[i], names[insert_point]) != 0)){
			// if the string is different, save it
			insert_point++;
			if(useCounts) counts[insert_point] = 1;
			if(insert_point != i) memcpy(names[insert_point], names[i], MAX_NODE_NAME_SIZE);
			cur_len = strlen(names[insert_point]);

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
	uint i=0, j=0, previ=0, prevj=j;
	uint len1=strlen(input_strings[i]), len2=strlen(file_strings[j]);
	int cmp;
	while(i < ninput && j < nfile){
		if(i != previ){
			previ=i;
			len1 = strlen(input_strings[i]);
		}
		if(j != prevj){
			prevj=j;
			len2 = strlen(file_strings[j]);
		}
		// skip past anything we've already found
		if(!bitarray[i] || len1 < len2){
			i++;
		} else if(len1 > len2){
			j++;
		} else {
			cmp = strcmp(input_strings[i], file_strings[j]);
			bitarray[i] = !!cmp;
			i += cmp <= 0;
			j += cmp >= 0;
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
	char fname[PATH_MAX], bitarray[num_to_sort], hashcontainer[5];

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

	uint cum_offsets[num_unique_hashes+1];
	cum_offsets[0] = 0;
	for(int i=0; i<num_unique_hashes; i++) cum_offsets[i+1] = cum_offsets[i] + filecounts[i];

	char **tmp_charptr;

	// now the new length of the array is insert_point
	// note that these are first sorted by hash, so once the hash value changes it will never appear again
	for(int i=0; i<num_unique_hashes; i++){
		// advance string array to lines for this file
		tmp_charptr = &(names[cum_offsets[i]]);

		// build filename using the hash
		snprintf(hashcontainer, 5, "%04x", hashes[i]);
		safe_filepath_cat(dir, hashcontainer, fname, PATH_MAX);
		//snprintf(fname, strlen(dir) + 6, "%s/%04x", dir, hashes[i]);

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
				if(bitarray[j]){
					cur_len = strlen(tmp_charptr[j]);
					fwrite(&cur_len, LEN_SIZE, 1, f);
					fwrite(tmp_charptr[j], 1, cur_len, f);
					fwrite(&ctr, L_SIZE, 1, f);
					ctr++;
				}
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
			if(!cmp) indices[i] = file_lines[j].ctr;
			i += cmp<=0;
			j += cmp>=0;
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
	char fname[PATH_MAX], hashcontainer[5];
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
		snprintf(hashcontainer, 5, "%04x", hashes[i]);
		safe_filepath_cat(dir, hashcontainer, fname, PATH_MAX);
		//snprintf(fname, strlen(dir) + 6, "%s/%04x", dir, hashes[i]);
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

	if(v) Rprintf("\tReading file %s...\n", fname);

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
		if(!(print_counter % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%lu lines read\r", print_counter);
			else R_CheckUserInterrupt();
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
		if(!(indices[i])) return total_count;
		fseek(ftabptr, (indices[i]-1)*L_SIZE, SEEK_SET);
		safe_fread(&cur_count, L_SIZE, 1, ftabptr);
		fseek(ftabptr, -1*L_SIZE, SEEK_CUR);
		cur_count += counts[i];
		total_count += counts[i];
		fwrite(&cur_count, L_SIZE, 1, ftabptr);
	}
	return total_count;
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
		normalizer = !normalizer;

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
			if(!strcmp(names_cache[i], names_cache_copy[j])){
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
		if(!(print_counter % PRINT_COUNTER_MOD)){
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
	R_CheckUserInterrupt();
	l_uint start, end, num_edges, tmp_cl, tmp_id, zeromaxid=ind+1;
	double tmp_w, acceptance_prob, zeromax=0;

	// move to information for the vertex and read in number of edges
	fseek(mastertab, L_SIZE*ind, SEEK_SET);
	safe_fread(&start, L_SIZE, 1, mastertab);
	safe_fread(&end, L_SIZE, 1, mastertab);

	// if we're above the max edges, subsample the edges we read
	num_edges = end - start;
	// if it has no edges we can't do anything
	if(!num_edges) return ind+1;

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

		/*
		 * this solves a special edge case where uninitialized nodes with a self-loop
		 * can be counted incorrectly if their neighbors have already been assigned to
		 * their node. Thus, if it's a self loop and the cluster is uninitialized,
		 * set it to its own cluster (what it would be initialized to)
		 */
		if(tmp_id == ind && !tmp_cl)
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
		if(!qsize) break;

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

	//reformat_clusters(clusterfile, num_v);
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

 	// first, index all vertex names
 	if(verbose) Rprintf("Building hash table for vertex names...\n");
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
	if(verbose) Rprintf("Calculating degree for each node...\n");
	for(int i=0; i<num_edgefiles; i++){
		edgefile = CHAR(STRING_ELT(FILENAME, i));
		hash_file_vnames_batch(edgefile, dir, temptabfile, seps[0], seps[1], 1, verbose, is_undirected, 1);
	}

 	// next, create an indexed table file of where data for each vertex will be located
 	if(verbose) Rprintf("Reformatting node degree file...\n");
 	reformat_counts(temptabfile, tabfile, num_v, add_self_loops);

 	// then, we'll create the CSR compression of all our edges
 	if(verbose) Rprintf("Reading in edges...\n");
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

 	// reindex the clusters from 1 to n
 	mergesort_clust_file(temptabfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu1, postcopy_dlu1);
 	mergesort_clust_file(temptabfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu2, postcopy_dlu2);

	SEXP RETVAL = PROTECT(allocVector(REALSXP, 1));
	REAL(RETVAL)[0] = (double) num_v;
	UNPROTECT(1);
	return RETVAL;
}

SEXP R_write_output_clusters(SEXP CLUSTERFILE, SEXP HASHEDDIR, SEXP OUTFILE, SEXP SEPS, SEXP VERBOSE){
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

	const int v = LOGICAL(VERBOSE)[0];
	const char* cfile = CHAR(STRING_ELT(CLUSTERFILE, 0));
	const char* hashdir = CHAR(STRING_ELT(HASHEDDIR, 0));
	const char* outfile = CHAR(STRING_ELT(OUTFILE, 0));
	const char* seps = CHAR(STRING_ELT(SEPS, 0));
	const uint fnamelen = strlen(hashdir) + 6;
	uint16_t name_len;
	int LEN_SIZE = sizeof(uint16_t);
	char buf[MAX_NODE_NAME_SIZE], fnamebuf[5];
	char write_buf[PATH_MAX], fname[PATH_MAX];
	l_uint index, clust, num_written=0;
	FILE *outf = fopen(outfile, "w");
	FILE *f;
	if(v) Rprintf("Writing node clusters to output file...\n");

	// remember that indices will be off by one
	FILE *fclusters = fopen(cfile, "rb");
	for(uint32_t i=0; i<(1ULL << NUM_BITS_HASH); i++){
		// build all possible file paths, if they don't exist then skip
		snprintf(fnamebuf, 5, "%04x", i);
		safe_filepath_cat(hashdir, fnamebuf, fname, fnamelen);
		f = fopen(fname, "rb");
		if(!f) continue;
		while(fread(&name_len, LEN_SIZE, 1, f)){
			// read in data from node names
			memset(buf, 0, MAX_NODE_NAME_SIZE);
			memset(write_buf, 0, PATH_MAX);
			safe_fread(buf, 1, name_len, f);
			safe_fread(&index, L_SIZE, 1, f);

			// get the cluster for the node name
			fseek(fclusters, L_SIZE*(index-1), SEEK_SET);
			safe_fread(&clust, L_SIZE, 1, fclusters);

			// prepare data for output
			snprintf(write_buf, (name_len+3)+L_SIZE, "%s%c%llu%c", buf, seps[0], clust, seps[1]);
			fwrite(write_buf, 1, strlen(write_buf), outf);
			num_written++;

			if(!(num_written % PRINT_COUNTER_MOD)){
				if(v) Rprintf("%lu nodes written.\r", num_written);
				else R_CheckUserInterrupt();
			}
		}
		fclose(f);
	}
	if(v) Rprintf("%lu nodes written.\n", num_written);
	fclose(fclusters);
	fclose(outf);

	return R_NilValue;
}