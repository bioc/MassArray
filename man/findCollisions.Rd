\name{findCollisions}
\alias{findCollisions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find peak collisions }
\description{
  Function to determine which subset(s) of peaks collide with eachother (i.e. molecular weights are indistinguishable given the specified resolution) 
}
\usage{
findCollisions(peaks, resolution = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{ Numerical list of molecular weights (in Da) corresponding to a set of peaks }
  \item{resolution}{ Resolution (in Da), used to specify the ability to distinguish two different molecular weights. For a resolution of \option{0.5} (default), two molecular weights are considered identical if they are less than \option{0.5} Da apart.}
}
\value{
  Returns a list of peak collisions for each peak in the original list, thus the data object returned is in the form of a nested list.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
findCollisions(1:5, 1.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }
\keyword{ manip }
