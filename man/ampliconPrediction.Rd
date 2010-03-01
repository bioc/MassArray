\name{ampliconPrediction}
\alias{ampliconPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Amplicon prediction }
\description{
  Function to predict amplicon fragmentation pattern and details for T&C reactions on the plus and minus strands
}
\usage{
ampliconPrediction(sequence, lower.threshold = 1500, upper.threshold = 7000, fwd.tag = "AGGAAGAGAG", rev.tag = "AGCCTTCTCCC", plot = TRUE, table = TRUE, lwd = 1, cex = 1, multiple.conversion = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{ Nucleotide sequence input as a character string }
  \item{lower.threshold}{Lower limit (in Da) of usable mass window (default: \option{1500})}
  \item{upper.threshold}{Upper limit (in Da) of usable mass window (default: \option{7000})}
  \item{fwd.tag}{Nucleotide tag sequence 5' of the forward primer}
  \item{rev.tag}{T7-containing nucleotide tag sequence 5' of the reverse primer}
  \item{plot}{ Logical specifying whether or not to display graphical representation of fragmentation profiles (default is \code{TRUE})}
  \item{table}{ Logical specifying whether or not to return tabular representation of fragmentation profiles (default is \code{TRUE}) }
  \item{lwd}{The line width used for fragmentation display, a positive number, defaulting to \code{1}}
  \item{cex}{A numerical value (defaulting to \code{1}) giving the amount by which plotting text and symbols should be magnified relative to the default}
  \item{multiple.conversion}{Logical value specifying whether or not to include multiple CGs on the same conversion control fragment where possible (default is \code{FALSE}).}
}
\details{
  Plotted fragmentation patterns contain a number of detailed features including: CG positions, molecular weight overlaps, conversion controls, fragment assayability, and more.

  Note that the graphical output does not contain a built-in legend at this time, but the plot may be interepreted as follows: 
  Putative fragmentation patterns are shown for T and C-cleavage reactions on both the plus and minus strands of an input amplicon, with the T-forward, T-reverse, C-forward, and C-reverse 
  shown in descending order. CG dinucleotides (filled circles) are numbered and colored in blue. Other fragments are colored according to their ability to be assayed:  fragment molecular 
  weight outside the testable mass window (gray), fragment molecular weight overlapping with another fragment (red), fragment containing a potential conversion control (green), or 
  fragment uniquely assayable but containing no CGs (black). Linked arrowheads denote molecular weight overlaps between multiple CG-containing fragments. Yellow highlights represent 
  tagged or primer sequences, while lavender highlights denote user-specified "required" sites.
}
\value{
  If \code{table} is \code{TRUE}, returns a list containing the following items:
  \item{summary }{A summary matrix of logical values specifying whether or not each CG is assayable by a given combination of cleavage reaction and DNA strand}
  \item{counts}{A numerical tally of the quantity of CGs that are assayable by each assay}
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu})  }
\examples{
ampliconPrediction("TGGAACACCCAGCAAAGATCAAGCAGGAAAGGGCGCACGCAGCCTTCGTTGCTAACCTCCTCTGGACTCTGGTACCCCAGGCACCGCGAATGCTCCCCACCTCAGCCCCCTGACCTTTACCATCGCTGAAGCGGGCGTCGCTGATGTCTGCGGCGAGCCTGCCGACCAGCCCAGCTGCCCAGAGGAGCAGCCAGGCAAGGGCGCTGGCAGCCAGGACGCCGGAGCCCGACGCCCGAGAGGGGCGCGCGGAGCAAGCTGCGGTCACGGGAGGAACCTGAGCACGCAGAGCGTACCCCCACCTTCCACGGTGACCCGGACAGAACGCTCCTTGCGCTCCCACCCTAGGACCCCCTGTAACTCCAGGTTCCTGAGA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{attribute }
\keyword{ character }
\keyword{ aplot }
