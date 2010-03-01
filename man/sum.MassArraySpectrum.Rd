\name{sum.MassArraySpectrum}
\alias{sum.MassArraySpectrum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sum MassArraySpectrum objects }
\description{
  Function to collapse multiple MassArraySpectrum objects into a single MassArraySpectrum representing the sum of each
}
\usage{
\method{sum}{MassArraySpectrum}(x, ..., trim = 0, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ One or multiple MassArraySpectrum objects to include in sum }
  \item{\dots}{ Any additional MassArraySpectrum objects to include in sum }
  \item{trim}{ Numerical value between 0 and 0.5 specifying the proportion of spectra to remove from consideration on a per peak basis, such that the SNR of each peak is calculated as the trimmed mean of the same peak across all included spectra. }
  \item{na.rm}{ Logical value passed to \code{mean}, indicating whether \code{NA} values should be stripped before the computation proceeds.  }
}
\value{
  Returns a single \code{MassArraySpectrum} object that represents the union of all unique peaks from the component MassArraySpectrum objects, with SNR for each peak representing the average value of that peak across all spectra
}
\author{  Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
MassArray.example.data$samples[[1]]$peaks[[11]]$height
MassArray.example.data$samples[[1]] <- sum.MassArraySpectrum(MassArray.example.data$samples[[1]], MassArray.example.data$samples[[2]])
MassArray.example.data$samples[[1]]$peaks[[11]]$height
}
\seealso{ See Also as \code{\linkS4class{MassArraySpectrum}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
