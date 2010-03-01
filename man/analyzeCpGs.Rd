\name{analyzeCpGs}
\alias{analyzeCpGs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Analyze CG methylation }
\description{
  Function to determine percent methylation for all CGs from input fragmentation
}
\usage{
analyzeCpGs(fragments, peaks, method = c("weighted", "proportion"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fragments}{ List of MassArrayFragment objects }
  \item{peaks}{ List of MassArrayPeak objects comprising spectral data for a given assay}
  \item{method}{ Specifies which algorithm to use when calculating percent methylation (either "weighted" or "proportion") }
}
\details{
  Wrapper function for \code{calcMeth()}, takes fragmentation pattern and spectral data as input and applies percent methylation calculation for all CG-containing, non conversion control fragments
}
\value{
  Returns a list of numerical values corresponding to percent methylation for each CG dinucleotide, with 0% indicating no methylation and 100% indicating complete methylation.
}
\author{  Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
cpg.data <- analyzeCpGs(MassArray.example.data$fragments.T, MassArray.example.data$samples[[1]]$peaks, method="weighted")
barplot(cpg.data, xlab="CpG (Number)", ylim=c(0,1), ylab="Methylation (Percent)")
}
\seealso{ See Also \code{\link{calcMeth}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
\keyword{ arith }
