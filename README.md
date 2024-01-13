# machineRy: Machine Learning Models for R

This is somewhat of a sandbox for me to write machine learning models that will be incorporated into other packages. You're welcome to clone/download this repo and use it in R, but bear in mind that this is **not** a stable build. Use at your own risk.

Final versions will be added to the SynExtend package when they're at a point where I'm confident they work as advertised and won't crash R. 

Current implementations:

- Feed-forward multilayer perceptron: sort of works, pretty buggy
- Random Forest: works for classification with numerical inputs, working on regression and categorical features

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

