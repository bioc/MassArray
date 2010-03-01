\name{convControl}
\alias{convControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Conversion control }
\description{
  Function to identify all potential conversion controls in a given input sequence, for a given list of fragments
}
\usage{
convControl(sequence, fragmentation, multiple = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{ Nucleotide sequence input as a character string }
  \item{fragmentation}{ List of MassArrayFragment objects corresponding to the fragmentation pattern of the sequence input }
  \item{multiple}{ Logical value specifying whether or not to include multiple CGs on the same conversion control fragment where possible (default is \code{FALSE}). }
}
\details{
 Potential conversion controls are identified from the nucleotide sequence input through pattern recognition of fragments that contain non-CG cytosines.
 Any conversion controls that contain CG dinucleotides or have a molecular weight outside of the usable mass window are screened out. 
  Additionally, conversion controls that are located in identified primer sequence or have molecular weight overlap with other fragments are removed from consideration. 
  Lastly, if the consideration of the fragment as a conversion control will cause new molecular weight overlap(s) with one or more other fragments, the control is also removed from consideration.
}
\value{
  Returns a list of MassArrayFragment objects identical to the input, with the exception that conversion controls are labeled and updated accordingly.
}
\author{  Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
MassArray.example.data$fragments.T[[54]]
conversion.data <- convControl(MassArray.example.data$sequence, MassArray.example.data$fragments.T)
conversion.data[[54]]
}
\seealso{ See Also as \code{\linkS4class{MassArrayFragment}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ character }
