\name{combine}
\alias{combine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Combine MassArrayData objects }
\description{
  Function to join two MassArrayData objects by sequence positions and samples
}
\usage{
combine(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ MassArrayData object }
  \item{y}{ MassArrayData object }
  \item{\dots}{ Other arguments passed to combine not supported at this time.  }
}
\value{
  Returns a single MassArrayData object that contains a union of samples and amplicons and spectral data from both MassArrayData objects in input
}
\author{  Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
samples(MassArray.example.data)
combined.data <- combine(MassArray.example.data[2,], MassArray.example.data[1,])
samples(combined.data)
}
\seealso{ See Also \code{\linkS4class{MassArrayData}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
