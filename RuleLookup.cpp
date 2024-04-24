#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calcASMR(NumericVector rSMR,NumericVector CMD, DataFrame Rules) {
  NumericVector rSMRLev = Rules["SMRLevel"];
  NumericVector CMDctf = Rules["CMD"];
  NumericVector aSMR = Rules["aSMR"];
  int len = CMD.length();
  IntegerVector v;
  NumericVector out(len);
  for(int i = 0; i < len; i++){
    v = Rcpp::seq(0, aSMR.length()-1);
    if(rSMR[i] < 5){
      v = v[rSMRLev == 0];
    }else if(rSMR[i] == 5){
      v = v[rSMRLev == 5];
    }else if(rSMR[i] == 6){
      v = v[rSMRLev == 6];
    }else{
      v = v[rSMRLev == 7];
    }
    int j;
    for(j = v[0]; j < v[v.length()]; j++){
      if(CMD[i] <= CMDctf[j]){
        break;
      }
    }
    out[i] = aSMR[j];
  }
  return(out);
}

