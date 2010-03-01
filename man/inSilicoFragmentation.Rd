\name{inSilicoFragmentation}
\alias{inSilicoFragmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ In silico fragmentation }
\description{
 Function to perform a complete in silico fragmentation of input sequence. Provides wrapper to a number of different functions, each of which determines additional information about each fragment.
}
\usage{
inSilicoFragmentation(sequence, fwd.tag = "", rev.tag = "", type = c("T", "C"), lower.threshold = 1500, upper.threshold = 7000, fwd.primer = 0, rev.primer = 0, multiple.conversion = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{ Nucleotide sequence input as a character string }
  \item{fwd.tag}{Nucleotide tag sequence 5' of the forward primer}
  \item{rev.tag}{T7-containing nucleotide tag sequence 5' of the reverse primer}
  \item{type}{ One of 'T' or 'C' indicating which cleavage reaction to use }
  \item{lower.threshold}{Lower limit (in Da) of usable mass window (default: \option{1500})}
  \item{upper.threshold}{Upper limit (in Da) of usable mass window (default: \option{7000})}
  \item{fwd.primer}{Length (in bp) of forward primer}
  \item{rev.primer}{Length (in bp) of reverse primer}
  \item{multiple.conversion}{Logical value specifying whether or not to include multiple CGs on the same conversion control fragment where possible (default is \code{FALSE}).}
}
\details{
  In silico fragmentation analysis includes RNAse A digestion, peak mapping and overlap detection, CG detection, assayability and conversion controls.
}
\value{
  Returns a list of \code{MassArrayFragment} objects, each with extensive contextual and other information
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
inSilicoFragmentation("GGGTTAGTCC")
}
\seealso{ See Also as \code{\linkS4class{MassArrayFragment}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip }
\keyword{character }
\keyword{ attribute }
