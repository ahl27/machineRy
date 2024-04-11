# Algorithm: FLP out of memory

## Big Picture and Definitions

Input: tsv file of n edges with v vertices. Each row is `v1 v2 weight`.

Cache is fixed to k entries, k>1.

The algorithm has five overall parts:

1. Read in all the unique vertex names and reindex them
2. Determine the degree of each vertex
3. Record each edge in CSR format
4. Run clustering
5. Lookup each vertex name to report results

Steps 1-3 each require reading through the entire edgelist once, so best possible runtime is `O(3n)`.


## Runtime of reindexing vertices

k is the cache size, which is fixed to a constant value at the beginning of the program.
n is the number of lines in the edgefile.

```
function MainFunction:
	cache vCache
	for line in edgefile:
		while vCache is not full:
			add vertex name to vCache

	if vCache is full:
		IndexCache(vCache)
		clear(vCache)
end function

function IndexCache(cache vCache):
	// WLOG we can assume all entries hash to the same hash file
	// having additional files means we reduce the runtime by a constant factor
	// recording hashes is done during the call to unique(vCache)
	hash = hash_name(vCache[0])
	vCache = unique(vCache) // O(k)
	vCache = sort(vCache) // O(k log k)
	bool shouldWrite[length(vCache)]

	cache fileCache
	open(hash_file)
	for line in fileCache:
		while fileCache is not full:
			add vertex name to fileCache

		if fileCache is full:
			fileCache = sort(fileCache) // O(k log k) (already unique'd)
			find_matching_entries(vCache, fileCache, shouldWrite) // see below for discussion
			clear(fileCache)

		if not any(bitarray):
			return

	for i in seq_along(vCache) where shouldWrite[i]:
		append vCache[i] to file

	close(file)
end function

```

Now to break down the runtimes.

`IndexCache` makes calls to a `find_matching_entries` function, which takes in the two caches and a bitarray. Since both `vCache` and `fileCache` are sorted, we can traverse both arrays simultaneously, marking `shouldWrite[i] = FALSE` iff `vCache[i] == fileCache[j]` for some `i,j`. The runtime for this operation is `O(l1 + l2)` since we only traverse each array once (`l1,l2` the lengths of `vCache, fileCache`). However, both caches are of fixed maximum size, and thus the worst case runtime is always `O(2k)` for this operation. Since we also sort both caches, the total runtime to compare two caches together is `O(2k log k)`.

However, we don't just compare two caches. Instead, we compare `vCache` to as many caches as it takes to iterate through the entire indexed file. This file is at most as long as the number of unique vertices `v`. This means that we need to compare `vCache` to `ceil(v/k)` caches, meaning the overall runtime is `O((1+ceil(v/k))(k log k))`. This function is called for as many times as it takes to iterate over the entire edge list, meaning another `ceil(n/k)` calls. Note that `ceil(n/k) < n/k + 1`. Thus, the worst case runtime is bounded by `O((n/k + 1 + v/k + 1)(k log k)) = O((n/k + v/k)(k log k))`.

Since the cache size `k` is a constant, `k log k` is also a constant at runtime. Thus `O((n/k + v/k)(k log k)) = O(n/k + v/k) = O(n + v)`. Distributing nodes amongst multiple files does not affect runtime scaling as long as each file is sufficiently large to fill a cache.

It's also important to note that file open/close operations are orders of magnitude slower than writes, which are themselves significantly slower than read and seek operations. Running `IndexCache` on a single file takes roughly the same amount of time to execute whether the cache has 1 element or 4000. This disparity is even worse because all writes are sequential appends in this format, which are faster than random writes. It's unclear at what number of writes this ceases to be true, but future scaling tests can investigate it once I have a completely working implementation.

## Runtime of steps 2-3

Edge counts and CSR compression follow the same rough idea as step 1. The same cache system is used; minor differences change the coefficient of scaling, but the complexity of the operations should scale identically. Step 2 could theoretically be bundled in with step 1, but the implementation to do that is relatively challenging and didn't seem to be worth the additional effort at the moment. Maybe I can come back to it later.

## Other notes

We had previously discussed using a hash table to compare the file names against the vertex names. That implementation likely makes sense, but it was a lot more difficult than this implementation, so I went with the current one for the time being. Much like bundling steps 1 and 2, this can be done in future work.