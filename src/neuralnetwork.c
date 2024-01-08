#include "machineRy.h"
#include "neuralnetwork.h"

/************** TODOs **************/
/* Add way to save weights from R  */
/* Maybe remove loss function      */
/***********************************/

/**** R Interface Functions ****/
SEXP R_initNNptr(SEXP NLAYERS, SEXP LAYERSIZES, SEXP ACTIVFUNCS, SEXP INPUT_SIZE, SEXP LR, SEXP LFXN){
  int nlayers = INTEGER(NLAYERS)[0];
  int *afuncR = INTEGER(ACTIVFUNCS);
  char afuncs[nlayers];
  char tmp;
  for(int i=0; i<nlayers; i++){
    switch(afuncR[i]){
    case 0:
      tmp = 'R';
      break;
    case 1:
      tmp = 'S';
      break;
    case 2:
      tmp = 'P';
      break;
    case 3:
      tmp = 'T';
      break;
    case 4:
    default:
      tmp = 'I';
    }
    afuncs[i] = tmp;
  }

  switch(INTEGER(LFXN)[0]){
  case 0:
    tmp = 'A';
    break;
  case 1:
  default:
    tmp = 'S';
  }
  NN* network = allocNetwork(nlayers, INTEGER(LAYERSIZES), afuncs,
                            INTEGER(INPUT_SIZE)[0], REAL(LR)[0], tmp);

  SEXP retval = PROTECT(R_MakeExternalPtr(network, mkString("machineRy_MLP_Model"), R_NilValue));
  R_RegisterCFinalizerEx(retval, (R_CFinalizer_t) FreeNetwork, TRUE);

  UNPROTECT(1);
  return(retval);
}

SEXP R_TrainForInput(SEXP INPUTVECTOR, SEXP OUTPUTVECTOR, SEXP nnPtr){
  if (!R_ExternalPtrAddr(nnPtr)) return R_NilValue;
  NN *network = (NN *) R_ExternalPtrAddr(nnPtr);
  double *inputvector = REAL(INPUTVECTOR);
  double *outputvector = REAL(OUTPUTVECTOR);

  int lastLayer = network->num_layers-1;
  int out_size = network->sizes[lastLayer];

  PredictSample(inputvector, network);
  UpdateLoss(outputvector, network);

  SEXP retval = PROTECT(allocVector(REALSXP, out_size*2));
  double *vptr = REAL(retval);
  for(int i=0; i<out_size; i++){
    vptr[i] = network->loss[i];
    vptr[i+out_size] = network->lossderiv[i];
  }

  UNPROTECT(1);
  return retval;
}

SEXP R_PredictForInput(SEXP INPUTVECTOR, SEXP nnPtr){
  if (!R_ExternalPtrAddr(nnPtr)) return R_NilValue;
  NN *network = (NN *) R_ExternalPtrAddr(nnPtr);
  double *inputvector = REAL(INPUTVECTOR);

  int lastLayer = network->num_layers-1;
  int outval = network->sizes[lastLayer];

  SEXP retval = PROTECT(allocVector(REALSXP, outval));
  double *vptr = REAL(retval);

  PredictSample(inputvector, network);
  memcpy(vptr, network->layers[lastLayer]->state, sizeof(double)*outval);
  UNPROTECT(1);
  return retval;
}

SEXP R_UpdateWeights(SEXP LOSSVEC, SEXP nnPtr){
  if (!R_ExternalPtrAddr(nnPtr)) return R_NilValue;
  NN *network = (NN *) R_ExternalPtrAddr(nnPtr);
  double *lossvector = REAL(LOSSVEC);

  BackpropNetworkBatch(lossvector, network);

  return R_NilValue;
}

void FreeNetwork(SEXP nnPtr){
  if (!R_ExternalPtrAddr(nnPtr)) return;
  NN *network = (NN *) R_ExternalPtrAddr(nnPtr);
  deallocNetwork(network);
  R_ClearExternalPtr(nnPtr);
  return;
}


/**** Allocation Functions ****/
NNLayer* allocLayer(int in, int out, char fxnChoice){
  NNLayer* newLayer = calloc(1, sizeof(NNLayer));
  newLayer->in_size = in;
  newLayer->out_size = out;
  AFunc *f;
  switch(fxnChoice){
  case 'R':
    f = &ReLU;
    break;
  case 'S':
    f = &Sigmoid;
    break;
  case 'P':
    f = &Softplus;
    break;
  case 'T':
    f = &Tanh;
    break;
  case 'I':
    f = &ID;
    break;
  }
  newLayer->activation = f;

  newLayer->weights = calloc(in*out, sizeof(double));
  newLayer->state = calloc(out, sizeof(double));
  newLayer->bias = calloc(out, sizeof(double));

  // randomly initialize weights in [-1,1]
  for(int i=0; i<in*out; i++)
    newLayer->weights[i] = (rand() % 2 ? -1 : 1) * (float)rand() / (float)(RAND_MAX);

  for(int i=0; i<out; i++)
    newLayer->bias[i] = (rand() % 2 ? -1 : 1) * (float)rand() / (float)(RAND_MAX);

  return newLayer;
}

NN* allocNetwork(int num_layers, int* layer_sizes, char* a_functions, int input_size, double lr, char lfxn){
  NN *network = calloc(1, sizeof(NN));

  num_layers++; // add an extra input layer at the beginning
  network->num_layers = num_layers;
  network->input_size = input_size;
  network->learning_rate = lr;

  switch(lfxn){
  case 'A':
    network->loss_function = &AbsError;
    break;
  case 'S':
  default:
    network->loss_function = &SqError;
  }

  network->sizes = calloc(num_layers, sizeof(int));
  network->layers = calloc(num_layers, sizeof(NNLayer*));

  int tmpin, tmpout;
  tmpin = input_size;
  for(int i=0; i<num_layers; i++){
    if(i == 0){
      // input layer
      network->sizes[i] = tmpin;
      network->layers[i] = allocLayer(tmpin, tmpin, 'I');
    } else {
      // all other layers
      tmpout = layer_sizes[i-1];
      network->sizes[i] = tmpout;
      network->layers[i] = allocLayer(tmpin, tmpout, a_functions[i-1]);
      tmpin = tmpout;
    }
  }

  network->loss = calloc(tmpout, sizeof(double));
  network->lossderiv = calloc(tmpout, sizeof(double));

  return network;
}

void deallocLayer(NNLayer* layer){
  free(layer->weights);
  free(layer->state);
  free(layer->bias);
  free(layer);
  return;
}

void deallocNetwork(NN* network){
  //Rprintf("Deallocating...\n");
  for(int i=0; i<network->num_layers; i++)
    deallocLayer(network->layers[i]);
  //Rprintf("Deallocating layers array...\n");
  free(network->layers);
  //Rprintf("Deallocating sizes array...\n");
  free(network->sizes);
  //Rprintf("Deallocating network...\n");
  free(network);
  return;
}
/******************************/

/********** Training **********/
void RunLayerForward(double *ipt, NNLayer* layer){
  int cols = layer->out_size;
  int rows = layer->in_size;
  double *w = layer->weights;
  double v;

  // Can be parallelized if necessary
  for(int i=0; i<cols; i++){
    v = 0;
    for(int j=0; j<rows; j++){
      v += ipt[j] * w[rows*i + j];
    }
    //Rprintf("%f %f\n", v, layer->activation->func(v));
    v = layer->activation->func(v+layer->bias[i]);
    layer->state[i] = v;
  }

  // now layer->state holds the new value
  return;
}

void BackpropLayer(NNLayer* l1, NNLayer* l2, double learning_rate){
  // l1 is current layer, l2 is next (closer to output) layer
  // l2 has already been backprop'd, now we update l1

  // equation for backprop:
  // t(l2.weights) * l2.err (x) deriv(l1.output)
  // (x) is Hadamard product
  // we can just store this in state, so the final equation is:
  int rows = l2->in_size;
  int cols = l2->out_size;


  // weights of l2 are o1 x o2 (o2 columns, o1 rows)
  // transpose should traverse rows, columns

  // each value becomes sum(t(l2.weights)[,i] * l2.err) * deriv(l1.state[i])
  // weight i,j updated as learning_rate * input[j] * error[i]
  double tmp, s, change;
  double *l2err = l2->state;
  double *l2w = l2->weights;
  for(int j=0; j<rows; j++){
    tmp = 0;
    s = l1->state[j];
    // multiply by row of transpose
    for(int i=0; i<cols; i++){
      tmp += l2err[i] * l2w[rows*i + j];
      //Rprintf("%f\n", learning_rate * s * l2err[i]);
      change = learning_rate * s * l2err[i];
      if(change < 0) change = fmax(change, -100);
      if(change > 0) change = fmin(change, 100);
      l2w[rows*i + j] -= change;

      if(j == 0)
        l2->bias[i] -= learning_rate*l2err[i];
    }

    // find derivative of output
    l1->state[j] = l1->activation->deriv(l1->state[j]);
    tmp *= l1->state[j];

    // store result
    l2->state[j] = tmp;
  }


  return;
}

void BackpropNetwork(double* expOutput, NN* nn){
  int nlayers = nn->num_layers-1;
  double lr = nn->learning_rate;

  // start with the last layer
  nn->loss_function->deriv(nn->layers[nlayers]->state, expOutput, nn->sizes[nlayers]);

  // iterate back to front
  // THIS CANNOT BE PARALLELIZED
  for(int i=nlayers-1; i>0; i--){
    BackpropLayer(nn->layers[i], nn->layers[i+1], lr);
  }

  return;
}

void BackpropNetworkBatch(double* avgLoss, NN* nn){
  int nlayers = nn->num_layers-1;
  double lr = nn->learning_rate;
  //int outsize = nn->sizes[nlayers];
  NNLayer* lastLayer = nn->layers[nlayers];
  for(int i=0; i<nn->sizes[nlayers]; i++){
    lastLayer->state[i] = avgLoss[i];
  }

  for(int i=nlayers-1; i>0; i--){
    BackpropLayer(nn->layers[i], nn->layers[i+1], lr);
  }

  /*
  for(int i=0; i<outsize; i++){
    nn->loss[i] = 0;
    nn->lossderiv[i] = 0;
  }
  */

  return;
}

void PredictSample(double* input, NN* nn){
  NNLayer* inputlayer = nn->layers[0];
  for(int i=0; i<nn->sizes[0]; i++)
    inputlayer->state[i] = input[i];

  double* tmp;
  for(int i=1; i<nn->num_layers; i++){
    tmp = nn->layers[i-1]->state;
    RunLayerForward(tmp, nn->layers[i]);
  }

  return;
}

void UpdateLoss(double* exp, NN* nn){
  int outsize = nn->sizes[nn->num_layers-1];
  double* lastState = nn->layers[nn->num_layers-1]->state;
  for(int i=0; i<outsize; i++){
    nn->loss[i] = lastState[i];
    nn->lossderiv[i] = lastState[i];
  }

  nn->loss_function->func(nn->loss, exp, outsize);
  nn->loss_function->deriv(nn->lossderiv, exp, outsize);

  return;
}

void TrainSample(double* input, double* output, NN* nn){
  PredictSample(input, nn);
  BackpropNetwork(output, nn);
  return;
}

/******************************/


/***** Testing *****/
void PrintLayer(NNLayer* layer){
  int rows = layer->in_size;
  int cols = layer->out_size;
  Rprintf("  -   Input: %d\n", rows);
  Rprintf("  -  Output: %d\n", cols);
  Rprintf("  -   AFunc: %s\n", layer->activation->label);
  Rprintf("  -   State: {");
  for(int i=0; i<cols; i++){
    Rprintf("%.2f", layer->state[i]);
    if(i+1 != cols)
      Rprintf(", ");
  }
  Rprintf("}\n");
  Rprintf("  -    Bias: {");
  for(int i=0; i<cols; i++){
    Rprintf("%.2f", layer->bias[i]);
    if(i+1 != cols)
      Rprintf(", ");
  }
  Rprintf("}\n");
  Rprintf("  - Weights:\n");
  for(int i=0; i<rows; i++){
    Rprintf("\t[");
    for(int j=0; j<cols; j++){
      Rprintf("%5.2f", layer->weights[j*rows + i]);
      if(j+1 != cols) Rprintf(", ");
    }
    Rprintf("]\n");
  }
  return;
}

void PrintNetwork(NN* nn){
  Rprintf("A network with %d layers.\n", nn->num_layers);
  Rprintf("Input size: %d\n", nn->input_size);
  Rprintf("Output size: %d\n", nn->sizes[nn->num_layers-1]);
  Rprintf("Learning rate: %f\n\n", nn->learning_rate);
  for(int i=0; i<nn->num_layers; i++){
    Rprintf("**** Layer %2d ****\n", i+1);
    PrintLayer(nn->layers[i]);
    Rprintf("******************\n");
  }
  return;
}

double testfxn(double* v){
  return(fSigmoid(v[0] + v[1] * v[2] + 0.5*v[3] + v[4]));
}

void testAfxn(AFunc af){
  for(double i=-10; i<10.1; i+=0.5){
    Rprintf("%6.2f %5.2f %5.2f\n", i, af.func(i), af.deriv(i));
  }
  return;
}

double genRandFloat(int max){
  return (rand() % 2 ? -1 : 1) * (float)rand() / (float)(RAND_MAX/max);
}

void testnetwork(){
  srand(time(NULL));
  //srand(11);
  int num_layers = 4;
  int layer_sizes[] = {10, 10, 10, 1};
  char a_functions[] = {'R', 'R', 'R', 'S'};
  int input_size = 5;
  int out_size = 1;
  double lr = 0.001;
  char lfxn = 'S';

  Rprintf("Allocating network...\n");
  NN* nn = allocNetwork(num_layers, layer_sizes, a_functions, input_size, lr, lfxn);

  int training_iter = 10000;
  int epoch = 10;
  int testing_iter = 100;
  int maxVal = 100;

  double ipt[input_size];
  double expOut[out_size];
  double loss[out_size];
  double lossderiv[out_size];
  double iptsum;
  //PrintNetwork(nn);
  Rprintf("\n************\n* Training  *\n************\n");
  for(int i=0; i<training_iter; i++){
    if((i+1) % epoch == 0){
      Rprintf("Average Loss: [");
      for(int j=0; j<out_size; j++){
        loss[j] /= (double)epoch;
        Rprintf("%.2f", loss[j]);
        if(j+1 < out_size) Rprintf(", ");
      }
      Rprintf("]\n");
      BackpropNetworkBatch(lossderiv, nn);
      for(int j=0; j<out_size; j++){
        loss[j] = 0;
        lossderiv[j] = 0;
      }
    }
    iptsum = 0;
    for(int j=0; j<input_size; j++){
      ipt[j] = genRandFloat(maxVal);
      iptsum += fabs(ipt[j]);
    }
    for(int j=0; j<input_size; j++)
      ipt[j] /= iptsum;

    expOut[0] = testfxn(ipt);
    PredictSample(ipt, nn);
    UpdateLoss(expOut, nn);
    for(int j=0; j<out_size; j++){
      if(isnan(nn->layers[nn->num_layers-1]->state[j])){
        PrintNetwork(nn);
        return;
      }
      loss[j] += nn->loss[j];
      lossderiv[j] += nn->lossderiv[j];
    }
  }

  //PrintNetwork(nn);
  for(int i=0; i<out_size; i++)
    loss[i] = 0;

  PrintNetwork(nn);

  Rprintf("\n\n***********\n* Testing *\n***********\n");
  for(int i=0; i<testing_iter; i++){
    for(int j=0; j<input_size; j++){
      ipt[j] = genRandFloat(maxVal);
    }
    expOut[0] = testfxn(ipt);
    PredictSample(ipt, nn);

    UpdateLoss(expOut, nn);
    for(int j=0; j<out_size; j++){
      loss[j] += nn->loss[j];
    }
  }
  Rprintf("Prediction Error: [");
  for(int j=0; j<out_size; j++){
      Rprintf("%.2f", loss[j]/testing_iter);
      if(j+1 != out_size){
        Rprintf(", ");
      }
    }
    Rprintf("]\n");

  Rprintf("\n\nDeallocating resources...\n");
  deallocNetwork(nn);
  return;
}
/*
int main(){
  testnetwork();
  //testAfxn(Sigmoid);
  return 0;
}
*/
