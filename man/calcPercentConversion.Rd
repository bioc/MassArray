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
  Returns a list of numerical values (from 0 to 1) corresponding to percent methylation for each conversion control, with 0% indicating complete conversion and 1 (100%) indicating a complete failure of bisulphite conversion. 
  Note that each element within the returned list will represent conversion control(s) for a single sample, while each element may contain multiple values with each value corresponding to data obtained from a single conversion control fragment.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
conversion.data <- calcPercentConversion(MassArray.example.data$fragments.T, MassArray.example.data$samples[[1]]$peaks)
mean(conversion.data)

# NOTE: conversion control data may already be contained within a MassArrayData object; these data can be accessed and visualized by the following (or alternative) commands
conversion.data <- unlist(lapply(lapply(MassArray.example.data$samples, slot, "quality.conversion"), median, na.rm=TRUE))
barplot(conversion.data)
}
\seealso{ See Also \code{\link{calcMeth}}, \code{\link{convControl}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
\keyword{ arith }
