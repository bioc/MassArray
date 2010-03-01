\name{evaluateSNPs}
\alias{evaluateSNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Evaluate SNPs }
\description{
  Function to analyze a MassArrayData object for all potential single nucleotide polymorphisms (SNPs) indicated by new and/or missing peaks in the the spectral data for one or more samples
}
\usage{
evaluateSNPs(x, verbose = TRUE, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ MassArrayData object containing spectral data for one or more samples }
  \item{verbose}{ Logical specifying whether or not to display descriptive progress updates as SNPs are analyzed }
  \item{plot}{ Logical specifying whether or not to display graphical representation of fragmentation profiles (default is \code{TRUE})}
}
\details{
  This function performs an exhaustive search for all potential SNPs (single base pair substitutions or deletions) that may give rise to new and/or missing peaks.  Graphical output is displayed 
  by default, and extensive data describing putative SNPs is also returned.
  
  Note that the graphical output does not contain a built-in legend at this time, but the plot may be interepreted as follows: 
  In the uppermost panel the T-cleavage fragmentation profile is shown for a given amplicon (C-cleavage reactions occupy a split screen whenever relevant). CG dinucleotides (filled circles) are numbered and colored in blue. Other fragments are colored 
  according to their ability to be assayed:  fragment molecular weight outside the testable mass window (gray), fragment molecular weight overlapping with another fragment (red), fragment 
  containing a potential conversion control (green), or fragment uniquely assayable but containing no CGs (black). Putative SNPs are shown directly below their location within the amplicon 
  fragmentation profile. Each row represents analysis from a single sample. Small, gray symbols represent potential SNPs that do not have sufficient evidence (presence of a new peak with 
  corresponding absence of an expected peak). Larger black symbols indicate a potential SNP with both new peaks and missing expected peaks. Triangles indicate base pair substitution, while 
  circles indicate single base pair deletion.
}
\value{
	Returns a list of potential SNPs for each identified new peak in the spectral data.  Note that each new peak may be explained by any number of potential SNPs; the list returned only includes the most reliable hits, but the redundant nature of the data necessitates returning a nested list, such that each new peak is associated with the following list elements:
  \item{SNP}{Contains a list of SNPs, each of which takes the form "\code{position}:\code{base}" where \code{position} is the base pair location within the amplicon sequence, and \code{base} is the mutated character}
  \item{SNR}{Contains a numerical list of signal-to-noise ratios corresponding to the expected original peak for the fragment mapping to the identified SNP position}
  \item{fragment}{Contains a numerical list of fragment IDs which map the SNP position to a specific fragment}
  \item{SNP.quality}{Contains a numerical list (values ranging from \code{0} to \code{2}, with \code{0} being a highly unlikely SNP and \code{2} being a SNP with increased likelihood.  This number is calculated as a function of new peak SNR and expected peak SNR.}
  \item{samples}{Contains a list of samples whose spectral data contained the given new peak}
  \item{count}{Specifies the number of unique SNP and sample pairs, exactly equivalent to the length of \code{SNP}, \code{SNR}, \code{fragment}, \code{SNP.quality}, or \code{samples}}
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
SNP.data <- evaluateSNPs(MassArray.example.data[2,])
}
\seealso{ See Also \code{\link{identifySNPs}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
\keyword{ manip }
