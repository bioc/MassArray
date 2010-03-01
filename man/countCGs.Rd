\name{countCGs}
\alias{countCGs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Count number of CGs}
\description{
  Function to count the number of CG dinucleotides in a given sequence (can include special characters for degenerate bases - i.e. 'Y' or 'R')
}
\usage{
countCGs(sequence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{ Nucleotide sequence input as a character string }
}
\value{
  Returns a numerical count of the number of CG dinucleotides in a given sequence, \option{NA} if sequence input is \option{NA}
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
countCGs("AAACGCGAAAAAAAYGAAA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ character }
