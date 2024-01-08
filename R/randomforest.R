# Random forests:
# 1. Train decision tree in C
# 2. Save decision tree to R
# 3. Predict w decision tree in R/C
# 4. on the fly pointer checking
# 5. Create multiple trees
# 6. Parallelism? Fork processes to get vectors, then regrow the trees? Depends on regrow runtime
#
# Do we need C? Use nx4 matrix
#
# For each tree:
#   - Struct has 4 elements: left, right, index (int), threshold (double)
# - If value at index is less than threshold, go left, else right
# - If index is -1, return threshold
# - Save with two vectors: int and numeric
# - Assume a BFS traversal, populate starting at root, if index >= 0, add to queue 2x
# - RF object will be a list of decision tree structures
# - Decision tree structure is a list with two vectors (above) and external pointer
# - If external pointer does not exist, reinitialize the object using vectors
