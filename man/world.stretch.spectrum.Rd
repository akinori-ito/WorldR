\name{world.stretch.spectrum}
\alias{world.stretch.spectrum}
\docType{package}
\title{
Stretch/shring the spectrum of the speech signal
}
\description{
\code{world.stretch.time} changes the tempo of the speeh without changing the F0.

}
\usage{
world.stretch.spectrum(worldobj,rate)
}
\arguments{
\item{worldobj}{A list calculated by \code{world.analysis()}.}
\item{rate}{A float value for stretching the spectrum.}
}
\details{
\code{world.stretch.spectrum()} stretches/compresses the spectrum of the speech. The argument \code{rate} is the stretch rate.

The return value is a list of parameters analyzed by World.
}
\seealso{
\code{\link{world.analysis}}, \code{\link{world.synthesis}}, \code{\link{world.stretch.time}}, \code{\link{world.stretch.spectrum}}
}
\examples{
\dontrun{

}
}
