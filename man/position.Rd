\name{position}
\alias{position}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{position<-}
\title{ Operate on positional information }
\description{
  Function to access (and/or assign) positional information for a MassArrayData object
}
\usage{
position(object)
position(object) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{MassArrayData}}
  \item{value}{ Character string containing positional information of the form "chrXX:XXXX-XXXX" }
}
\value{
  Returns a character string containing positional information of the form "chrXX:XXXX-XXXX" if accessing a \code{MassArrayData} object. If updating a \code{MassArrayData} object, function returns the object with updated positional information
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
position(MassArray.example.data)
position(MassArray.example.data) <- "chrB:2001-2374"
position(MassArray.example.data)
}
\seealso{ See Also \code{\linkS4class{MassArrayData}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
