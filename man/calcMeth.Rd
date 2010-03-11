\name{calcMeth}
\alias{calcMeth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate percent methylation }
\description{
  Function to calculate percent methylation from a collection of peaks corresponding to a single fragment.
}
\usage{
calcMeth(peaks, peaklist, method = c("weighted", "proportion"), num.cg.fragments = rep(1, length(peaks)), non.cg.fragments = rep(0, length(peaks)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{ List of molecular weights, sorted from low to high, corresponding to the unmethylated and methylated peaks for a given fragment, respectively }
  \item{peaklist}{ List of MassArrayPeak objects comprising spectral data for a complete assay}
  \item{method}{ Specifies which algorithm to use when calculating percent methylation (either \code{"weighted"} or \code{"proportion"}) }
  \item{num.cg.fragments}{ Number of CG-containing fragments contributing to each of the input peaks, automatically defaulting to \code{1} }
  \item{non.cg.fragments}{ Number of fragments (without CGs) contributing to each of the input peaks, automatically defaulting to \code{0} }
}
\details{
  Note that the current release of this function performs as expected for the large majority of cases.  However, certain complex combinations of peak overlaps are not handled at this time.  
  This may affect data for a minority of points, particularly those containing multiple overlaps with alternative fragments.  Please ensure more in-depth review of such loci.
}
\value{
  Returns a numerical values corresponding to percent methylation, with 0% indicating no methylation and 100% indicating complete methylation.
}
\references{   Coolen, M.W., et al. (2007) Genomic profiling of CpG methylation and allelic specificity using quantitative high-throughput mass spectrometry: critical evaluation and improvements, \emph{Nucleic Acids Research}, \strong{35}(18), e119. }
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
frags <- MassArray.example.data$fragments.T[[6]]$"MW"
calcMeth(frags, MassArray.example.data$samples[[1]]$peaks, "weighted")
calcMeth(frags, MassArray.example.data$samples[[1]]$peaks, "proportion")
}
\seealso{ See Also \code{\linkS4class{MassArrayPeak}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ arith}
\keyword{ manip }
