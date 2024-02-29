# machineRy: Machine Learning Models for R

This is somewhat of a sandbox for me to write machine learning models that will be incorporated into other packages. You're welcome to clone/download this repo and use it in R, but bear in mind that this is **not** a stable build. Use at your own risk.

Final versions will be added to the SynExtend package when they're at a point where I'm confident they work as advertised and won't crash R. 

Current implementations:

- Feed-forward multilayer perceptron: sort of works, pretty buggy
- Random Forest: works for classification with numerical inputs, working on regression and categorical features
- Fast Label Propagation: Works for `igraph` graphs. Consensus clustering across weight differentials 
also working.
- Out of memory clustering: Works for arbitrary sets of edgelists stored in `v1 v2 w` format. `tsv` format is preferred, but support for arbitrary encodings exists. Help file will be written soon.

# Current Stats:

## Random Forest vs. `randomForest`:

Testing done by training/predicting on a matrix of four numeric variables and one categorical response variable.
So far, `machineRy` trains slower but predicts faster than `randomForest`. Working on further optimization.

Memory usage is about the same, but the final implementation will likely be better than `randomForest`.

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

Performance is roughly identical on weighted LFR graph benchmarks (measured using AMI). Runtime of this implementation is roughly 2x faster than `igraph`.
Most of the runtime is dedicated to converting `igraph`-style graphs into something I can work with; final implementation for `SynExtend` should be significantly
faster since I won't be working with `igraph` graphs.

Consensus clustering runs slower than either `igraph` or my label propagation, mainly because multiple 
runs are required. This implementation beats `igraph` and my LP algorithm in accuracy. Runtime is about 
10x slower because it does about 10 LP runs. Scaling on all algorithms is approximately linear (0.76 for 
my LP, 0.94 for `igraph`, 1.33 for consensus clustering).
