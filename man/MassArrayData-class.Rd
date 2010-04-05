\name{MassArrayData-class}
\docType{class}
\alias{MassArrayData-class}
\alias{$,MassArrayData-method}
\alias{$<-,MassArrayData-method}
\alias{[,MassArrayData-method}
\alias{initialize,MassArrayData-method}

\title{Class "MassArrayData" }
\description{A data structure containing MassArray data and associated information for a single amplicon}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MassArrayData", sequence, file, verbose, fwd.tag, rev.tag, fwd.primer, rev.primer, strand, lower.threshold, upper.threshold, header, skip, sep, comment.char, fill, method, position, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{sequence}:}{Nucleotide sequence for unconverted amplicon}
    \item{\code{chr}:}{Chromosomal position of amplicon}
    \item{\code{start}:}{Chromosomal position of amplicon}
    \item{\code{end}:}{Chromosomal position of amplicon}
    \item{\code{strand}:}{DNA strand used for amplicon (can be '+' or '-')}
    \item{\code{fwd.tag}:}{Nucleotide tag sequence 5' of the forward primer}
    \item{\code{rev.tag}:}{T7-containing nucleotide tag sequence 5' of the reverse primer}
    \item{\code{fwd.primer}:}{Length (in bp) of forward primer}
    \item{\code{rev.primer}:}{Length (in bp) of reverse primer}
    \item{\code{lower.threshold}:}{Lower limit (in Da) of usable mass window (default: \option{1500})}
    \item{\code{upper.threshold}:}{Upper limit (in Da) of usable mass window (default: \option{7000})}
    \item{\code{fragments.T}:}{List containing objects of class \code{MassArrayFragment}, corresponding to the T-cleavage reaction for the amplicon on the specified \code{strand}}
    \item{\code{fragments.C}:}{List containing objects of class \code{MassArrayFragment}, corresponding to the C-cleavage reaction for the amplicon on the specified \code{strand}}
    \item{\code{samples}:}{List containing object of class \code{MassArraySpectrum}, each corresponding to spectral data from a single sample}
    \item{\code{groups}:}{List of the group name to which each sample belongs}
    \item{\code{CpG.data}:}{Matrix containing analyzed methylation data, where each row is a sample and each column is a CG dinucleotide site}
    \item{\code{CpG.data.combined}:}{Matrix containing methylation data combined from multiple objects (or collapsed from within a single object), where each row is a sample and each column is a CG dinucleotide site}
  }
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "MassArrayData")}: ... }
    \item{\$<-}{\code{signature(x = "MassArrayData")}: ... }
    \item{[}{\code{signature(x = "MassArrayData")}: ... }
    \item{initialize}{\code{signature(.Object = "MassArrayData")}: ... }
	 }
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
showClass("MassArrayData")
}
\keyword{classes}
