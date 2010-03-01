\name{numCollisions}
\alias{numCollisions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Count peak collisions }
\description{
   Function to count the number of peak collisions (i.e. molecular weights are indistinguishable given the specified resolution)
}
\usage{
numCollisions(peaks, resolution = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{ Numerical list of molecular weights (in Da) corresponding to a set of peaks }
  \item{resolution}{ Resolution (in Da), used to specify the ability to distinguish two different molecular weights. For a resolution of \option{0.5} (default), two molecular weights are considered identical if they are less than \option{0.5} Da apart.}
}
\value{
  Returns a list of peak collision counts for each peak in the original list.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
numCollisions(1:5, 1.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }
\keyword{ manip }
