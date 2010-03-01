\name{calcPercentConversion}
\alias{calcPercentConversion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate percent conversion }
\description{
  Function to calculate percent methylation (wrapper for \code{calcMeth()} function) for each identified conversion control
}
\usage{
calcPercentConversion(fragments, peaks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fragments}{ List of MassArrayFragment objects }
  \item{peaks}{ List of MassArrayPeak objects comprising spectral data to be used for conversion control calculations }
}
\details{
  This function serves as a wrapper function for \code{calcMeth()}, such that percent methylation is calculated for all conversion controls within the input list of fragments.
}
\value{
  Returns a list of numerical values corresponding to percent methylation for each conversion control, with 0% indicating complete conversion and 100% indicating a complete failure of bisulphite conversion.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
conversion.data <- calcPercentConversion(MassArray.example.data$fragments.T, MassArray.example.data$samples[[1]]$peaks)
mean(conversion.data)
}
\seealso{ See Also \code{\link{calcMeth}}, \code{\link{convControl}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
\keyword{ arith }
