\name{findFragments}
\alias{findFragments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find fragments }
\description{
  Function to identify which fragment(s) in a list of fragments match a given molecular weight
}
\usage{
findFragments(MW, fragments, resolution = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MW}{ Molecular weight target (in Da) }
  \item{fragments}{ List of molecular weights corresponding to unique fragments }
  \item{resolution}{ Resolution (in Da), used to specify the ability to distinguish two different molecular weights. For a resolution of \option{1} (default), two molecular weights are considered identical if they are less than \option{1} Da apart.}
}
\value{
  Returns the index or indices of fragment(s) within the input list that have a molecular weight which matches that specified as input 
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
findFragments(3913,MassArray.example.data$fragments.T, resolution=0.1)
findFragments(3913,MassArray.example.data$fragments.T, resolution=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }
\keyword{ manip }
