\name{world.analysis}
\alias{world.analysis}
\docType{package}
\title{
Analyze speech signal using World
}
\description{
\code{world.analysis} analyzes the input speech signal into F0, spectrum and aperiodicity.

}
\usage{
world.analysis(wave,f0=NULL,frameshift=5.0,f0floor=71.0,allowed_range=0.1)
}
\arguments{
\item{wave}{A \code{\link{Wave}} object to analyze.}
\item{f0}{A vector of F0 sequence. If NULL, the F0 is calculated in the function.}
\item{frameshift}{Value of frame shift given in ms.}
\item{f0floor}{The lower limit of F0.}
\item{allowed_range}{}
}
\details{
\code{world.analysis()} analyzes the given \code{\link{Wave}} object into F0, spectrum and aperiodicity. If the Wave object has two or more channels, only the first channel (i.e. the left channel) will be used.

If \code{f0} is specified, \code{world.analysis()} uses the given F0 for spectral analysis. Otherwise, it calculates F0 using Dio and StoneMask algorithms.

The return value is a list, having the following members:
\describe{
\item{samp.rate}{The sampling rate (Hz).}
\item{frameshift}{The frame period (ms).}
\item{length}{The number of frames.}
\item{F0}{A vector of F0 values (Hz).}
\item{spec}{A matrix of spectrum. The [i,j]-th element denotes the j-th frequency element of the i-th frame.}
\item{aperiodicity}{A matrix of aperiodicity, having the same dimension as the spectrum.}
}
}
\seealso{
\code{\link{world.synthesis}}, \code{\link{world.f0}}
}
\examples{
\dontrun{

}
}
