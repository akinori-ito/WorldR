\name{world.synthesis}
\alias{world.synthesis}
\docType{package}
\title{
Synthesize speech signal using World
}
\description{
\code{world.synthesis} synthesizes the speech signal from the F0, spectrum and aperiodicity obtained by \code{world.analysis()}.

}
\usage{
world.synthesis(worldobj,normalize=FALSE)
}
\arguments{
\item{worldobj}{A list calculated by \code{world.analysis()}.}
\item{normalize}{if TRUE, the output wave is normalized so that the maximum amplitude to be 32767.}
}
\details{
\code{world.synthesis()} synthesizes the speech signal from the analysis result obtained by \code{\link{world.analysis}}.

The return value is an \code{\link{Wave}} object, having one channel, PCM format. The sampling rate of the output will be the same as that of the input signal to the \code{world.analysis()}.
}
\seealso{
\code{\link{world.analysis}}
}
\examples{
\dontrun{

}
}
