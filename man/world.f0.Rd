\name{world.f0}
\alias{world.f0}
\docType{package}
\title{
Extract fundamental frequencies (F0) using World
}
\description{
\code{world.f0} extracts the fundamental frequency (F0) sequence using the World's standard F0 extraction algorithms, Dio and StoneMask.
}
\usage{
world.f0(wave,f0=NULL,frameshift=5.0,f0floor=71.0,allowed_range=0.1)
}
\arguments{
\item{wave}{A \code{\link{Wave}} object to analyze.}
\item{f0}{A vector of F0 sequence. If NULL, the F0 is calculated in the function.}
\item{frameshift}{Value of frame shift given in ms.}
\item{f0floor}{The lower limit of F0.}
\item{allowed_range}{}
}
\details{
\code{world.f0()} analyzes the given \code{\link{Wave}} object and extracts F0.
If the Wave object has two or more channels, only the first channel (i.e. the left channel) will be used.

The return value is a vector of F0 sequence (in Hz).
}
\seealso{
\code{\link{world.analysis}}
}
\examples{
\dontrun{

}
}
