\name{isAssayable}
\alias{isAssayable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Is assayable? }
\description{
  Function to determine whether or not a given molecular weight is assayable (i.e. within the usable mass window specified)
}
\usage{
isAssayable(MW, lower.threshold = 1500, upper.threshold = 7000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MW}{ Numerical input corresponding to molecular weight }
  \item{lower.threshold}{Lower limit (in Da) of usable mass window (default: \option{1500})}
  \item{upper.threshold}{Upper limit (in Da) of usable mass window (default: \option{7000})}
}
\value{
  Returns a logical corresponding to whether or not the molecular weight input falls within the usable mass window specified
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
isAssayable(5000)
isAssayable(1200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ logic }
