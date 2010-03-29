\name{calcMeth}
\alias{calcMeth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate percent methylation }
\description{
  Function to calculate percent methylation from a collection of peaks corresponding to a single fragment.
}
\usage{
calcMeth(SNR.list, fragments = rep(1, length(SNR.list)), non.cg.fragments = numeric(0), method = c("weighted", "proportion"), prune.non.cg.peaks = TRUE, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNR.list}{ List of signal-to-noise ratios, sorted from low to high MWs, corresponding to the unmethylated and methylated peaks for a given set of fragments }
  \item{fragments}{ List of all fragments contributing to each of the input peaks, automatically defaulting to a single fragment }
  \item{non.cg.fragments}{ List of fragments (without CGs) contributing to any of the input peaks, automatically defaulting to \code{numeric(0)} }
  \item{method}{ Specifies which algorithm to use when calculating percent methylation (either \code{"weighted"} or \code{"proportion"}) }
  \item{prune.non.cg.peaks}{ Boolean value determining whether or not to remove non-CG-containing fragments prior to analysis or whether to include them in the calculating model (note that setting this option to \code{FALSE} could result in a considerable increase in analytical time); option automatically defaults to \code{TRUE} }
  \item{na.rm}{ Boolean value determing whether or not to return an error on input of any unspecified data (\code{NA}), automatically defaulting to \code{FALSE} } 
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
peaks <- findPeaks(frags, unlist(lapply(MassArray.example.data$samples[[1]]$peaks, slot, "MW.actual")))
SNRs <- unlist(lapply(MassArray.example.data$samples[[1]]$peaks[peaks], slot, "SNR"))
frag.list <- list(1:3, 1:3, 1:3, 1:3)
calcMeth(SNRs, fragments=frag.list, method="weighted")
calcMeth(SNRs, fragments=frag.list, method="proportion")
}
\seealso{ See Also \code{\linkS4class{MassArrayPeak}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ arith}
\keyword{ manip }
