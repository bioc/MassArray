\name{identifySNPs}
\alias{identifySNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identify SNPs }
\description{
  Function to identify potential single nucleotide polymorphisms (SNPs) which allow mapping of a novel peak sequence to the expected amplicon sequence
}
\usage{
identifySNPs(peak.sequence, sequence, rxn = c("T", "C"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak.sequence}{ Nucleotide sequence (can also be base composition - ex: "A6G2C1T3") as a character string }
  \item{sequence}{ Nucleotide sequence for wildtype/expected amplicon as a character string }
  \item{rxn}{ One of "T" or "C" specifying which cleavage reaction to use for SNP analysis }
}
\details{
  The algorithm steps through the \code{sequence}, substituting one nucleotide at a time with the other three base pairs or a blank character (deletion), in order to determine a base compositional match to the input \code{peak.sequence} which represents a peak not found in the native \code{sequence}.
}
\value{
  Returns a list of potential SNP matches for the input \code{peak.sequence}. Each element of the list contains multiple items as follows:
  \item{sequence}{ corresponds to \code{peak.sequence}}
  \item{position}{ corresponds to the matched position within \code{sequence}}
  \item{base}{ corresponds to the altered nucleotide (i.e. "A", "T", "C", "G", or "")}
  \item{type}{ corresponds to the class of SNP (i.e. "substitution" or "deletion")}
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
## SINGLE SUBSTITUTION
identifySNPs("AAGT","AATTTT")
## MULTIPLE SUBSTITUTION POSSIBILITIES
identifySNPs("A1G1T1","AATTTT")
## DELETION
identifySNPs("AAT","AGATTTT")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ character }
\keyword{ manip }
