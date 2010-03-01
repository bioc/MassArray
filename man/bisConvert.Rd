\name{bisConvert}
\alias{bisConvert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bisulphite conversion }
\description{
 Bisulphite convert nucleotide sequence input
}
\usage{
bisConvert(sequence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{ Nucleotide sequence in the form of a character string }
}
\value{
  Returns a character value corresponding to the bisulphite converted input sequence.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
bisConvert("AAATTCGGAACCC")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ character }
