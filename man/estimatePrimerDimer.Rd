\name{estimatePrimerDimer}
\alias{estimatePrimerDimer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate Primer Dimers }
\description{
  Function to estimate level of signal due to primer dimers in a given spectrum
}
\usage{
estimatePrimerDimer(fragments, peaks, method = c("ratio", "mann-whitney"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fragments}{ List of MassArrayFragment objects corresponding to the sample }
  \item{peaks}{ List of MassArrayPeak objects comprising spectral data for a complete assay }
  \item{method}{ Specifies which algorithm to use when estimating primer dimer levels (either \code{"ratio"} or \code{"mann-whitney"}) }
}
\details{
  Primer dimers are calculated by: 1) identifying fragments that occur within the expected primer sequence, 2) identifying which of these fragments is assayable, and 3) comparing the overall signal for primer peaks and peaks from the rest of the amplicon.
}
\value{
  Returns a list containing primer dimer ratios or significance estimates (i.e. p-values) depending on the analytical method specified (\code{"ratio"} or \code{"mann-whitney"}, respectively)
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
primer.data <- estimatePrimerDimer(MassArray.example.data$fragments.T, MassArray.example.data$samples[[1]]$peaks, "ratio")
mean(primer.data)
primer.data <- estimatePrimerDimer(MassArray.example.data$fragments.T, MassArray.example.data$samples[[1]]$peaks, "mann-whitney")
mean(primer.data)
}
\seealso{ See Also \code{\linkS4class{MassArrayData}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
