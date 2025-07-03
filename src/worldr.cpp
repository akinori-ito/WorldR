#include <Rcpp.h>
#include <math.h>
#include <world/d4c.h>
#include <world/dio.h>
#include <world/matlabfunctions.h>
#include <world/cheaptrick.h>
#include <world/stonemask.h>
#include <world/synthesis.h>


using namespace Rcpp;
List spectral_analysis(double *x, int x_length,
                       double *f0, double *time_axis, double f0floor,
                       int fs, int frameshift,
                       int length, NumericVector &Rf0, NumericVector &Rtime)
{
  // Spectral envelope estimation
  CheapTrickOption c_option = {0};
  InitializeCheapTrickOption(fs,&c_option);
  c_option.q1 = -0.15;
  c_option.f0_floor = f0floor;

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

  for (int t = 0; t < length; t++) {
    delete[] spectrogram[t];
    delete[] aperiodicity[t];
  }
  delete[] spectrogram;
  delete[] aperiodicity;

  return ret;
}

// [[Rcpp::export]]
List worldAnalysis_(NumericVector& wave, double frameshift, int fs,
                    double f0floor, double allowed_range) {

  DioOption d_option = {0};
  InitializeDioOption(&d_option);
  d_option.frame_period = frameshift;
  d_option.speed = 1;
  d_option.f0_floor = f0floor;
  d_option.allowed_range = allowed_range;

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
  delete[] refined;
  List ret = spectral_analysis(x, x_length,
                           f0,time_axis,f0floor,
                           fs, frameshift, length,
                           Rf0,Rtime);
  delete[] x;
  delete[] f0;
  delete[] time_axis;

  return ret;
}

// [[Rcpp::export]]
List worldAnalysis_f0(NumericVector& wave, NumericVector& Rf0,
                      double frameshift, int fs,
                      double f0floor, double allowed_range) {

  int x_length = wave.size();
  double *x = new double[x_length];
  for (int i = 0; i < x_length; i++)
    x[i] = wave(i);

  int length = Rf0.size();
  double *f0 = new double[length];
  double *time_axis = new double[length];

  NumericVector Rtime(length);
  for (int i = 0; i < length; i++) {
    f0[i] = Rf0(i);
    Rtime(i) = time_axis[i] = frameshift*i/1000;
  }

  List ret = spectral_analysis(x, x_length,
                               f0,time_axis,f0floor,
                               fs, frameshift, length,
                               Rf0,Rtime);
  delete[] x;
  delete[] f0;
  delete[] time_axis;

  return ret;
}

// [[Rcpp::export]]
NumericVector worldF0Estimation(NumericVector& wave, double frameshift, int fs,
                    double f0floor, double allowed_range) {

  DioOption d_option = {0};
  InitializeDioOption(&d_option);
  d_option.frame_period = frameshift;
  d_option.speed = 1;
  d_option.f0_floor = f0floor;
  d_option.allowed_range = allowed_range;

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
  for (int i = 0; i < length; i++) {
    Rf0(i) = refined[i];
  }
  delete[] refined;
  delete[] x;
  delete[] f0;
  delete[] time_axis;

  return Rf0;
}


// [[Rcpp::export]]
NumericVector worldSynthesis_(List world) {
  double frameshift = as<double>(world["frameshift"]);
  int fs = as<int>(world["samp.rate"]);
  int length = as<int>(world["length"]);
  NumericVector Rf0 = world["F0"];
  NumericMatrix Rspec = world["spec"];
  NumericMatrix Raperio = world["aperiodicity"];

  double *F0 = new double[length];
  double **spectro = new double*[length];
  double **aperio = new double*[length];
  int specsize = Rspec.ncol();
  for (int i = 0; i < length; i++) {
    spectro[i] = new double[specsize];
    aperio[i] = new double[specsize];
    F0[i] = Rf0(i);
    for (int j = 0; j < specsize; j++) {
      spectro[i][j] = Rspec(i,j);
      aperio[i][j] = Raperio(i,j);
    }
  }

  int y_length = static_cast<int>((length-1)*frameshift/1000.0*fs)+1;
  double *y = new double[y_length];
  int fft_size = (specsize-1)*2;

  Synthesis(F0,length,spectro,aperio,fft_size,frameshift, fs,
            y_length, y);

  NumericVector Ry(y_length);
  for (int i = 0; i < y_length; i++)
    Ry(i) = y[i];

  delete[] F0;
  for (int i = 0; i < length; i++) {
    delete[] spectro[i];
    delete[] aperio[i];
  }
  delete[] spectro;
  delete[] aperio;
  delete[] y;
  return Ry;
}

