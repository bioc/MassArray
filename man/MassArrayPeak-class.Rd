\name{MassArrayPeak-class}
\docType{class}
\alias{MassArrayPeak-class}
\alias{$,MassArrayPeak-method}
\alias{$<-,MassArrayPeak-method}
\alias{initialize,MassArrayPeak-method}

\title{Class "MassArrayPeak" }
\description{A data structure containing information and data for a single peak from a single spectrum}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MassArrayPeak", ID, MW.theoretical, MW.actual, probability, SNR, height, sample.intensity, ref.intensity, sequence, adduct, type, charge, collisions, components, missing, new, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ID}:}{Peak ID inicating indexed position within a potentially large list of peaks}
    \item{\code{MW.theoretical}:}{Expected molecular weight of peak based on nucleotide sequence}
    \item{\code{MW.actual}:}{Observed molecular weight from experimental data}
    \item{\code{probability}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{SNR}:}{Signal-to-noise ratio}
    \item{\code{height}:}{Raw peak height}
    \item{\code{sample.intensity}:}{Raw sample intensity}
    \item{\code{ref.intensity}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sequence}:}{Nucleotide composition or sequence(s) corresponding to peak}
    \item{\code{adduct}:}{One of 'Na', 'K', or '' indicating whether or not peak represents a salt adduct of another expected peak}
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{charge}:}{Degree of ionization of fragment (default is '1' indicating a single positive charge per fragment)}
    \item{\code{collisions}:}{Number of peaks that share the same molecular weight as the current peak}
    \item{\code{components}:}{Number of fragments expected to give rise to a peak of this molecular weight }
    \item{\code{missing}:}{Logical indicating whether or not the expected peak is missing from the spectral data }
    \item{\code{new}:}{Logical indicating whether or not the observed peak is unexpected given the amplicon sequence}
  }
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "MassArrayPeak")}: ... }
    \item{\$<-}{\code{signature(x = "MassArrayPeak")}: ... }
    \item{initialize}{\code{signature(.Object = "MassArrayPeak")}: ... }
	 }
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\seealso{
	See Also \code{\linkS4class{MassArrayData}}
}
\examples{
showClass("MassArrayPeak")
}
\keyword{classes}
