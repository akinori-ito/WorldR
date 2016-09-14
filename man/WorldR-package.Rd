\name{WorldR-package}
\alias{WorldR-package}
\alias{WorldR}
\docType{package}
\title{
WorldR: wrapper of World
}
\description{
WorldR is a wrapper of the World, the speech analysis and synthesis system.
}
\details{
It provides the following functionalities:
\itemize{
\item Analyze a speech signal (the Wave object) into F0, spectrum and aperiodicity (\code{\link{world.analysis}})
\item Synthesize a speech signal from the F0, spectrum and aperiodicity (\code{\link{world.synthesis}})
\item Extract F0 (\code{\link{world.f0}})
\item Change F0 of the speech without changing the tempo (\code{\link{world.shift.f0}})
\item Stretch and compress the speech signal without changing F0 (\code{\link{world.stretch.time}})
\item Change the vocal identity (virtually the length of the vocal tract) without changing F0 and tempo (\code{\link{world.stretch.spectrum}})
}
}
\author{
Akinori Ito

Maintainer: Akinori Ito <aito@fw.ipsj.or.jp>
}
\references{
 World: a high-quality speech analysis and synthesis system. http://ml.cs.yamanashi.ac.jp/world/english/index.html
}
\keyword{ package }
\seealso{
\code{\link{tuneR}}
}
\examples{
\dontrun{
library(tuneR)
library(WorldR)
wave <- readWave("wavefile.wav")
ana <- world.analysis(wave)
slow <- world.stretch.time(ana,0.5)
slowwave <- world.synthesis(slow)
writeWave(slowwave,"slowwave.wav")
  }
}
