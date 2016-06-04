#include <Rcpp.h>
#include <math.h>
#include <world/d4c.h>
#include <world/dio.h>
#include <world/matlabfunctions.h>
#include <world/cheaptrick.h>
#include <world/stonemask.h>
#include <world/synthesis.h>


using namespace Rcpp;

// [[Rcpp::export]]
List worldAnalysis(NumericVector& wave, double frameshift, int fs) {

  DioOption d_option = {0};
  InitializeDioOption(&d_option);
  d_option.frame_period = frameshift;
  d_option.speed = 1;
  d_option.f0_floor = 71.0;
  d_option.allowed_range = 0.1;

  // F0 analysis
  int x_length = wave.size();
  double *x = new double[x_length];
  for (int i = 0; i < x_length; i++)
    x[i] = wave(i);

  int length = GetSamplesForDIO(fs,x_length,frameshift); //Frame length
  double *f0 = new double[length];
  double *time_axis = new double[length];
  double *refined = new double[length];

  Dio(x, x_length, fs, &d_option, time_axis,f0);
  StoneMask(x, x_length, fs, time_axis, f0,length, refined);
  NumericVector Rf0(length);
  NumericVector Rtime(length);
  for (int i = 0; i < length; i++) {
    Rf0(i) = refined[i];
    Rtime(i) = time_axis[i];
  }

  // Spectral envelope estimation
  CheapTrickOption c_option = {0};
  InitializeCheapTrickOption(&c_option);
  c_option.q1 = -0.15;
  c_option.f0_floor = 71.0;

  int fftsize = GetFFTSizeForCheapTrick(fs, &c_option);
  int specsize = fftsize/2+1;
  double **spectrogram = new double*[length];
  for (int i = 0; i < length; i++) {
    spectrogram[i] = new double[specsize];
  }

  CheapTrick(x, x_length, fs, time_axis, f0, length, &c_option, spectrogram);
  NumericMatrix spec(length,specsize);
  for (int t = 0; t < length; t++) {
      for (int f = 0; f < specsize; f++) {
          spec(t,f) = spectrogram[t][f];
      }
  }

  // Aperiodicity estimation
  D4COption d4option = {0};
  InitializeD4COption(&d4option);
  double **aperiodicity = new double*[length];
  for (int i = 0; i < length; i++) {
    aperiodicity[i] = new double[specsize];
  }
  D4C(x, x_length, fs, time_axis,
      f0, length, fftsize, &d4option, aperiodicity);
  NumericMatrix aper(length,specsize);
  for (int t = 0; t < length; t++) {
    for (int f = 0; f < specsize; f++) {
      aper(t,f) = aperiodicity[t][f];
    }
  }

  List ret = List::create(
    Named("frameshift") = frameshift,
    Named("samp.rate") = fs,
    Named("length") = length,
    Named("F0") = Rf0,
    Named("time_axis") = Rtime,
    Named("spec") = spec,
    Named("aperiodicity") = aper);

  delete[] x;
  delete[] f0;
  delete[] time_axis;
  delete[] refined;
  for (int t = 0; t < length; t++) {
    delete[] spectrogram[t];
    delete[] aperiodicity[t];
  }
  delete[] spectrogram;
  delete[] aperiodicity;

  return ret;
}
