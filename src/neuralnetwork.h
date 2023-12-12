#ifndef NN_HEADER
#define NN_HEADER 0
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

typedef struct ActivationFunction{
  double (*func)(double v);
  double (*deriv)(double v);
  char* label;
} AFunc;

typedef struct LossFunction{
  // These functions operate in-place on the first vector
  void (*func)(double *obs, double *exp, int len);
  void (*deriv)(double *obs, double *exp, int len);
} Loss;

typedef struct layer{
  int in_size;
  int out_size;

  // necessarily out_size columns, in_size rows
  double *weights; // 2D implied array
  double *state;
  double *bias;

  // function pointer for activation function
  AFunc* activation;
} NNLayer;

typedef struct network{
  int num_layers;
  int input_size;
  int* sizes; // output sizes of each layer
  double learning_rate;
  double* loss;
  double* lossderiv;
  Loss* loss_function;
  NNLayer** layers;
} NN;

/**** Activation Functions ****/

/** Functions and Derivatives **/
inline static double fReLU(double v){
  return v < 0 ? 0 : v;
}

inline static double dReLU(double v){
  return (v > 0 ? 1 : 0);
}

inline static double fSigmoid(double v){
  return (1.0 / (1+exp(-1*v)));
}

inline static double dSigmoid(double v){
  double sv = fSigmoid(v);
  return (sv * (1-sv));
}

inline static double fSoftplus(double v){
  return log(1+exp(v));
}

inline static double dSoftplus(double v){
  return fSigmoid(v);
}

inline static double fTanh(double v){
  double pexp = exp(v);
  double nexp = exp(-1*v);
  return((pexp-nexp) / (pexp+nexp));
}

inline static double dTanh(double v){
  double ft = fTanh(v);
  return(1 - (ft*ft));
}

inline static double fID(double v){
  return v;
}

inline static double dID(double v){
  return 0;
}

/** Activation Structs **/
/* Key *
 * R = ReLU
 * S = Sigmoid
 * P = Softplus
 * T = Tanh
 * I = Identity
 */
static AFunc ReLU = {fReLU, dReLU, "ReLU"};
static AFunc Sigmoid = {fSigmoid, dSigmoid, "Sigmoid"};
static AFunc Softplus = {fSoftplus, dSoftplus, "Softplus"};
static AFunc Tanh = {fTanh, dTanh, "Tanh"};
static AFunc ID = {fID, dID, "Identity"};
/********************/

/** Loss Functions **/

// SHOULD BE TRUE - OUTPUT
static void fSqError(double *obs, double *exp, int len){
  for (int i=0; i<len; i++)
    obs[i] = 0.5*pow((exp[i] - obs[i]), 2);
  return;
}

static void dSqError(double *obs, double *exp, int len){
  for (int i=0; i<len; i++)
    obs[i] = exp[i] - obs[i];

  return;
}

static void fAbsError(double *obs, double *exp, int len){
  for (int i=0; i<len; i++)
    obs[i] = fabs(exp[i] - obs[i]);

  return;
}

static void dAbsError(double *obs, double *exp, int len){
  for (int i=0; i<len; i++)
    obs[i] = exp[i] > obs[i] ? 1.0 : -1.0;

  return;
}

/** Loss Structs **/
static Loss SqError = {fSqError, dSqError};
static Loss AbsError = {fAbsError, dAbsError};

/******************************/




/**** Allocation Functions ****/
NNLayer* allocLayer(int in, int out, char fxnChoice);

NN* allocNetwork(int num_layers, int* layer_sizes, char* a_functions, int input_size, double lr, char lfxn);

void deallocLayer(NNLayer* layer);

void deallocNetwork(NN* network);
/******************************/

/********** Training **********/
void RunLayerForward(double *ipt, NNLayer* layer);

void BackpropLayer(NNLayer* l1, NNLayer* l2, double learning_rate);

void BackpropNetwork(double* expOutput, NN* nn);

void BackpropNetworkBatch(double* avgLoss, NN* nn);

void PredictSample(double* input, NN* nn);

void TrainSample(double* input, double* output, NN* nn);

void UpdateLoss(double* exp, NN* nn);

/******************************/

/***** Testing *****/
void PrintLayer(NNLayer* layer);
void PrintNetwork(NN* nn);
void testnetwork(void);
/*******************/

#endif