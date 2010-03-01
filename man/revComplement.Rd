\name{revComplement}
\alias{revComplement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reverse complement }
\description{
  Function to find the reverse complement
}
\usage{
revComplement(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ sequence input to use for reverse complement. \code{x} can be a character string or a MassArrayData object. }
}
\value{
  Returns the reverse complement of a character string or MassArrayData object, depending upon input data type.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\seealso{ \code{\link{revComplement-methods}} }
\examples{
revComplement("AATCCGGGGGAA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ character }
