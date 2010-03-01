\name{findPeaks}
\alias{findPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find peaks }
\description{
  Function to determine which peak(s) in a list of peaks match a given molecular weight.
}
\usage{
findPeaks(MW, peaks, resolution = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MW}{ Molecular weight target (in Da) }
  \item{peaks}{ List of molecular weights corresponding to unique peaks }
  \item{resolution}{ Resolution (in Da), used to specify the ability to distinguish two different molecular weights. For a resolution of \option{1} (default), two molecular weights are considered identical if they are less than \option{1} Da apart.}
}
\value{
  Returns the index or indices of peak(s) within the input list that have a molecular weight which matches that specified as input 
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
findPeaks(3.1, 6:1, res=0)
findPeaks(3.1, 6:1, res=0.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }
\keyword{ manip }
