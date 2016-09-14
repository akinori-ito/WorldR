\name{world.shift.f0}
\alias{world.shift.f0}
\docType{package}
\title{
Shift F0 values of the speech signal
}
\description{
\code{world.shift.f0} changes the F0 obtained by \code{world.analysis()}.

}
\usage{
world.synthesis(worldobj,rate)
}
\arguments{
\item{worldobj}{A list calculated by \code{world.analysis()}.}
\item{rate}{A float value to be multiplied to the F0.}
}
\details{
\code{world.shift.f0()} multiplies \code{rate} to the F0 vaues contained in \code{worldobj}. If \code{rate=2.0}, it doubles the F0.

The return value is a list in which the F0 values are changed.
}
\seealso{
\code{\link{world.analysis}}, \code{\link{world.synthesis}}, \code{\link{world.stretch.time}}, \code{\link{world.stretch.spectrum}}
}
\examples{
\dontrun{

}
}
