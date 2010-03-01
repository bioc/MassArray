\name{MassArraySpectrum-class}
\docType{class}
\alias{MassArraySpectrum-class}
\alias{$,MassArraySpectrum-method}
\alias{$<-,MassArraySpectrum-method}
\alias{initialize,MassArraySpectrum-method}

\title{Class "MassArraySpectrum" }
\description{A data structure containing MassArray spectral data for a single sample}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MassArraySpectrum", sample, rxn, strand, peaks, quality.conversion, quality.spectra, quality.primerdimer, quality.contaminant, quality.adducts, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{sample}:}{Sample name}
    \item{\code{rxn}:}{Cleavage reaction (either 'T' or 'C')}
    \item{\code{strand}:}{DNA strand for amplicon (either '+' or '-')}
    \item{\code{peaks}:}{List containing objects of class \code{MassArrayPeak}}
    \item{\code{quality.conversion}:}{Overall level(s) of remnant unconverted cytosines, as measured by one or more conversion controls}
    \item{\code{quality.spectra}:}{(currently not supported)}
    \item{\code{quality.primerdimer}:}{(currently not supported)}
    \item{\code{quality.contaminant}:}{(currently not supported)}
    \item{\code{quality.adducts}:}{Overall ratio(s) of Na and/or K adduct peak heights to expected peak heights}
  }
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "MassArraySpectrum")}: ... }
    \item{\$<-}{\code{signature(x = "MassArraySpectrum")}: ... }
    \item{initialize}{\code{signature(.Object = "MassArraySpectrum")}: ... }
	 }
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\seealso{
	See Also \code{\linkS4class{MassArrayData}}
}
\examples{
showClass("MassArraySpectrum")
}
\keyword{classes}
