\name{samples}
\alias{samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{samples<-}
\title{ Operate on sample name(s) }
\description{
  Function to access (and/or assign) sample name information for a MassArrayData object
}
\usage{
samples(object)
samples(object) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{MassArrayData}}
  \item{value}{ List of character strings containing sample name information, one for each sample contained in the \code{MassArrayData} object }
}
\value{
  Returns a list of character strings containing sample name information for each sample in \code{MassArrayData} object. If updating a \code{MassArrayData} object, function returns the object with updated sample name information
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
samples(MassArray.example.data)
samples(MassArray.example.data)[2] <- "SECOND"
samples(MassArray.example.data)
}
\seealso{ See Also as \code{\linkS4class{MassArrayData}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
