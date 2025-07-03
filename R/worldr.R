#library(tuneR)

#' @title Analyze speech signal using World
#' @description \code{world.analysis()} analyzes the input speech signal into F0, spectrum and aperiodicity.
#'\code{world.analysis()} analyzes the given \code{\link{Wave}} object into F0, spectrum and aperiodicity. If the Wave object has two or more channels, only the first channel (i.e. the left channel) will be used.
#'If \code{f0} is specified, \code{world.analysis()} uses the given F0 for spectral analysis. Otherwise, it calculates F0 using Dio and StoneMask algorithms.
#' @param wave A \code{\link{Wave}} object to analyze.
#' @param f0 A vector of F0 sequence. If NULL, the F0 is calculated in the function.
#' @param frameshift Value of frame shift given in ms.
#' @param f0floor The lower limit of F0.
#' @param allowed_range The allowed range.
#' @return The return value is a list, having the following members:
#' - samp.rate: The sampling rate (Hz).
#' - frameshift: The frame period (ms).
#' - length: The number of frames.
#' - F0: A vector of F0 values (Hz).
#' - spec: A matrix of spectrum. The [i,j]-th element denotes the j-th frequency element of the i-th frame.
#' - aperiodicity}{A matrix of aperiodicity, having the same dimension as the spectrum.
#' @export
world.analysis <- function(w,f0=NULL,frameshift=5.0,f0floor=71.0,allowed_range=0.1) {
  if (is.null(f0)) {
    r <- worldAnalysis_(w@left,frameshift,w@samp.rate,f0floor,allowed_range)
  } else {
    r <- worldAnalysis_f0(w@left,f0,frameshift,w@samp.rate,f0floor,allowed_range)

  }
  r$class <- "World"
  return(r)
}

#' @title Extract fundamental frequencies (F0) using World
#' @description `world.f0` extracts the fundamental frequency (F0) sequence using the World's standard F0 extraction algorithms, Dio and StoneMask.
#' If the Wave object has two or more channels, only the first channel (i.e. the left channel) will be used.
#' @param wave A Wave object to analyze.
#' @param f0 A vector of F0 sequence. If NULL, the F0 is calculated in the function.
#' @param frameshift Value of frame shift given in ms.
#' @param f0floor The lower limit of F0 (Hz).
#' @param allowed_range The allowed range
#' @return A vector of F0 sequence.
#' @export
world.f0 <- function(w,frameshift=5.0,f0floor=71.0,f0ceil=800.0,allowed_range=0.1,method="dio") {
  if (method == "dio") {
    return(worldF0Estimation_dio(w@left,frameshift,w@samp.rate,f0floor,allowed_range))
  } else if (method == "harvest") {
    return(worldF0Estimation_dio(w@left,frameshift,w@samp.rate,f0floor,f0ceil))
  } else {
    stop("world.f0: unknown method ",method)
  }
}

#' @title Synthesize speech signal using World
#' @description `world.synthesis` synthesizes the speech signal from the F0, spectrum and aperiodicity obtained by `world.analysis()`.
#' @param world A list calculated by `world.analysis()`.
#' @param normalize if TRUE, the output wave is normalized so that the maximum amplitude to be 32767.
#' @return An Wave object, having one channel, PCM format. The sampling rate of the output will be the same as that of the input signal to the `world.analysis()`.
#' @export
world.synthesis <- function(world,normalize=FALSE) {
  w <- worldSynthesis_(world)
  w[is.nan(w)] <- 0
  mx <- max(abs(w))
  if (normalize | mx > 32767) {
    # do normalize
    w <- w/mx*32767
  }
  return(Wave(floor(w),samp.rate=world$samp.rate,bit=16))
}

#' @title Shift F0 values of the speech signal
#' @description `world.shift.f0` changes the F0 obtained by `world.analysis()`.
#' @param world A list calculated by `world.analysis()`.
#' @param rate A float value to be multiplied to the F0.
#' @return A list in which the F0 values are changed.
#' @export
world.shift.f0 <- function(world,rate) {
  w <- world
  w$F0 <- w$F0*rate
  return(w)
}

#' @title Stretch/shring the tempo of the speech signal
#' @description `world.stretch.time` changes the tempo of the speech without changing the F0.
#' @param world A list calculated by `world.analysis()`.
#' @param rate A float value for stretching
#' @return A list in which the tempo is changed.
#' @export
world.stretch.time <- function(world,rate) {
  w <- list()
  len <- world$length
  specsize <- ncol(world$spec)
  newlen <- floor(len/rate)
  w$frameshift <- world$frameshift
  w$samp.rate <- world$samp.rate
  w$length <- newlen
  w$F0 <- signal::interp1(0:(len-1),world$F0,
                          (0:(newlen-1))/(newlen-1)*(len-1),
                           method="linear")
  w$spec <- matrix(nrow=newlen,ncol=specsize)
  w$aperiodicity <- matrix(nrow=newlen,ncol=specsize)
  for (f in 1:specsize) {
    w$spec[,f] <- signal::interp1(0:(len-1),world$spec[,f],
                                  (0:(newlen-1))/(newlen-1)*(len-1),
                                  method="linear")
    w$aperiodicity[,f] <- signal::interp1(0:(len-1),world$aperiodicity[,f],
                                          (0:(newlen-1))/(newlen-1)*(len-1),
                                          method="linear")
  }
  w$class <- "World"
  return(w)
}

#' @title Stretch/shring the spectrum of the speech signal
#' @description `world.stretch.spectrum` changes the spectrum of the speech without changing the F0.
#' @param world A list calculated by `world.analysis()`.
#' @param rate A float value for stretching
#' @return A list in which the spectrum is changed.
#' @export
world.stretch.spectrum <- function(world,rate) {
  w <- list()
  len <- world$length
  specsize <- ncol(world$spec)
  w$frameshift <- world$frameshift
  w$samp.rate <- world$samp.rate
  w$length <- world$length
  w$F0 <- world$F0
  w$spec <- matrix(nrow=len,ncol=specsize)
  w$aperiodicity <- matrix(nrow=len,ncol=specsize)
  orgspecsize <- floor(specsize/rate)
  if (orgspecsize <= specsize) {
    for (t in 1:len) {
      w$spec[t,] <- signal::interp1(0:(orgspecsize-1),world$spec[t,],
                                   (0:(specsize-1))/(specsize-1)*(orgspecsize-1),
                                   method="linear")
      w$aperiodicity[t,] <- signal::interp1(0:(orgspecsize-1),world$aperiodicity[t,],
                                   (0:(specsize-1))/(specsize-1)*(orgspecsize-1),
                                   method="linear")
    }
  }
  else {
    for (t in 1:len) {
      s <- runif(orgspecsize)*1.0e-5
      s[1:specsize] <- world$spec[t,]
      w$spec[t,] <- signal::interp1(0:(orgspecsize-1),s,
                                    (0:(specsize-1))/(specsize-1)*(orgspecsize-1),
                                    method="linear")
      s[1:specsize] <- world$aperiodicity[t,]
      w$aperiodicity[t,] <- signal::interp1(0:(orgspecsize-1),s,
                                            (0:(specsize-1))/(specsize-1)*(orgspecsize-1),
                                            method="linear")
    }
  }
  w$class <- "World"
  return(w)
}
