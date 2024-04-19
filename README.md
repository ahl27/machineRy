# machineRy: Machine Learning Models for R

This is somewhat of a sandbox for me to write machine learning models
that will be incorporated into other packages. You're welcome to
clone/download this repo and use it in R, but bear in mind that this is
**not** a stable build. Use at your own risk.

Final versions will be added to the SynExtend package when they're at a
point where I'm confident they work as advertised and won't crash R.

Current implementations:

-   Feed-forward multilayer perceptron: sort of works, pretty buggy
-   Random Forest: works for classification with numerical inputs,
    working on regression and categorical features
-   Fast Label Propagation: Works for `igraph` graphs. Consensus
    clustering across weight differentials also working.
-   Out of memory clustering: Works for arbitrary sets of edgelists
    stored in `v1 v2 w` or `v1 v2` format. `tsv` format is preferred, but support
    for arbitrary encodings exists. Performance roughly matches FLP in
    accuracy.
    
Things I'm working on next:
-   OOM clustering speedups: slowest operations are reading in edges and sorting
    vertex names. I'm not yet sure how to optimize this step further.
-   Consensus clustering for OOMLP.

# Current Stats:

## Random Forest vs. `randomForest`:

Testing done by training/predicting on a matrix of four numeric
variables and one categorical response variable. So far, `machineRy`
trains slower but predicts faster than `randomForest`. Working on
further optimization.

Memory usage is about the same, but the final implementation will likely
be better than `randomForest`.

### Training (runtime, seconds):

```         
# of entries (10^x)   1.0   1.5   2.0   2.5   3.0   3.5   4.0   4.5   5.0
randomForest        0.003 0.001 0.001 0.003 0.007 0.025 0.092 0.367 1.415
machineRy           0.001 0.002 0.004 0.010 0.034 0.123 0.429 1.462 4.672
```

### Prediction (runtime, seconds):

```         
# of entries (10^x)   1.0   1.5   2.0   2.5   3.0   3.5   4.0   4.5   5.0
randomForest        0.001 0.001 0.001 0.001 0.002 0.006 0.016 0.052 0.164
machineRy           0.000 0.000 0.001 0.001 0.002 0.007 0.016 0.050 0.155

# of entries (10^x)   5.5   6.0   6.5    7.0    7.5     8.0
randomForest        0.522 1.871 6.177 20.545 68.069 301.170
machineRy           0.512 1.651 5.052 15.901 53.976 220.985
```

## Label Propagation vs. `igraph::cluster_label_prop`

Performance is roughly identical on weighted LFR graph benchmarks
(measured using AMI). Runtime of this implementation is roughly 2x
faster than `igraph`. Most of the runtime is dedicated to converting
`igraph`-style graphs into something I can work with; final
implementation for `SynExtend` should be significantly faster since I
won't be working with `igraph` graphs. **Note**: This is comparing my
implementation of the Fast Label Propagation algorithm to the standard
Label Propagation implemented in `igraph`. Fast Label Propagation exists
on an `igraph` fork, but not in the distributed library.

Consensus clustering runs slower than either `igraph` or my label
propagation, mainly because multiple runs are required. This
implementation beats `igraph` and my LP algorithm in accuracy. Runtime
is about 10x slower because it does about 10 LP runs.

## In-memory LP vs. Out-of-memory LP

Testing done on MacBook Pro with M1 Pro CPU and 32GB RAM. Out of memory runtime is dominated by data preprocessing rather than the clustering
itself. For example, clustering a network with a million nodes and edges takes 107s running for a single iteration, and 129s running for infinite. 
Consensus clustering is not yet implemented, but should require significantly less additional runtime than the in-memory implementation since subsequent OOM runs will only need to read in the graph a single time.

Results below are shown for a single iteration of FLP algorithms. Running for more will increase runtime but should not meaningfully impact the memory or runtime scaling.

Computational scaling is shown below. Values are exponent `n` for `O(x^n)`, where `x` is number of nodes/edges. In all trials used for this computation, the number of nodes equals the number of edges--future investigation can look into the scaling with respect to number of nodes/edges individually rather than together.

```
                         Runtime Scaling           Memory Scaling
       MCL, I=2.0              1.87                     0.79         
           igraph              1.75                     0.73
machineRy,  inmem              1.19                     0.87
machineRy, outmem              1.11                     0.00
```

1,000 node graph with 8,000 edges:

```         
                      Memory Usage (Max, MB)   Total Elapsed Time (sec)
           igraph              6.7                      0.2
machineRy,  inmem              7.7                      0.2
machineRy, outmem              4.8                      0.6
```

10,000 node graph with 50,000 edges:
```         
                      Memory Usage (Max, MB)   Total Elapsed Time (sec)
           igraph             18.9                      0.3
machineRy,  inmem             34.5                      0.3
machineRy, outmem             23.3                      3.1
```

100,000 node graph with 100,000 edges:
```         
                      Memory Usage (Max, MB)   Total Elapsed Time (sec)
       MCL, I=2.0            105.0                      6.7
           igraph             62.7                      1.3
machineRy,  inmem            106.1                      0.8
machineRy, outmem             65.9                      9.0
```

250,000 node graph with 250,000 edges:
```         
                      Memory Usage (Max, MB)   Total Elapsed Time (sec)
       MCL, I=2.0            262.0                     26.8
           igraph            123.0                      5.3
machineRy,  inmem            197.1                      1.6
machineRy, outmem             84.2                     23.9
```


1,000,000 node graph with 1,000,000 edges:
```         
                      Memory Usage (Max, MB)   Total Elapsed Time (min:sec)
       MCL, I=2.0            778.9                     5:03.2
           igraph            416.6                     1:15.6
machineRy,  inmem            756.9                     0:06.3
machineRy, outmem            103.5                     1:37.0
```

10,000,000 node graph with 10,000,000 edges:
```         
                      Memory Usage (Max, GB)   Total Elapsed Time (hr:min:sec)
       MCL, I=2.0             4.20                     7:15:45
           igraph             1.77                     1:07:10
machineRy,  inmem             5.57                     0:01:26
machineRy, outmem             0.08                     0:23:52
```
