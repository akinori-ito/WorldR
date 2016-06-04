#include <Rcpp.h>
#include <math.h>
#include <world/dio.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector analyze(NumericVector wave, double frameshift) {
  DioOption option = {0};
  InitializeDioOption(&option);
  option.frame_period = frameshift;
  option.speed = 1;


}
