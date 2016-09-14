\name{world.stretch.time}
\alias{world.stretch.time}
\docType{package}
\title{
Change the tempo of the speech signal
}
\description{
\code{world.stretch.time} changes the tempo of the speeh without changing the F0.

}
\usage{
world.stretch.time(worldobj,rate)
}
\arguments{
\item{worldobj}{A list calculated by \code{world.analysis()}.}
\item{rate}{A float value for changing the tempo.}
}
\details{
\code{world.stretch.time()} stretches/compresses the speech signal. The argument \code{rate} is the rate. If \code{rate=2.0}, the speech signal becomes two times as quick as the original signal.

The return value is a list of parameters analyzed by World.
}
\seealso{
\code{\link{world.analysis}}, \code{\link{world.synthesis}}, \code{\link{world.stretch.time}}, \code{\link{world.stretch.spectrum}}
}
\examples{
\dontrun{

}
}
