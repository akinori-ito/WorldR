#library(tuneR)

world.analysis <- function(w,f0=NULL,frameshift=5.0,f0floor=71.0,allowed_range=0.1) {
  if (is.null(f0)) {
    r <- worldAnalysis_(w@left,frameshift,w@samp.rate,f0floor,allowed_range)
  } else {
    r <- worldAnalysis_f0(w@left,f0,frameshift,w@samp.rate,f0floor,allowed_range)

  }
  r$class <- "World"
  return(r)
}

world.f0 <- function(w,frameshift=5.0,f0floor=71.0,allowed_range=0.1) {
  return(worldF0Estimation(w@left,frameshift,w@samp.rate,f0floor,allowed_range))
}

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

world.shift.f0 <- function(world,rate) {
  w <- world
  w$F0 <- w$F0*rate
  return(w)
}

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

world.test <- function() {
  library(tuneR)
  w<-readWave("src/World/test/vaiueo2d.wav")
  r <- world.analysis(w)
  return(r)
}
