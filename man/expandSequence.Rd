\name{expandSequence}
\alias{expandSequence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Expand nucleotide sequence }
\description{
  Function to process shorthand form of a nucleotide sequence, where a given base pair followed by a number specifies a run of the indicated nucleotide for the specified length (ex: "A6TTCGA4")
}
\usage{
expandSequence(sequence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{ Nucleotide sequence input as a character string }
}
\value{
  Returns an expanded nucleotide sequence as a character string
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
expandSequence("A6TTCGA4")
expandSequence("C10C10")
expandSequence("AT1CG")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ character }
\keyword{ manip }
