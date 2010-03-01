\name{importEpiTyperData.new}
\alias{importEpiTyperData.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Import EpiTyper data (v.1.0.5) }
\description{
  Function to read and import an EpiTyper datafile (v.1.0.5) and store it as a list of MassArraySpectrum objects
}
\usage{
importEpiTyperData.new(data, MassArrayObject, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ location of EpiTyper datafile as a character string }
  \item{MassArrayObject}{ Pre-existent MassArrayData object in which store relevant sample and spectral information from datafile }
  \item{verbose}{ Logical specifying whether or not to display descriptive progress updates as datafile is processed }
}
\details{
  EpiTyper v.1.0.5 datafiles must only contain a single amplicon, thus the user must export peak data for one amplicon at a time.
}
\value{
  Returns a list of MassArraySpectrum objects each populated by spectral data
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\seealso{ See Also \code{\linkS4class{MassArraySpectrum}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
\keyword{ file }
