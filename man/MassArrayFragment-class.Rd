\name{MassArrayFragment-class}
\docType{class}
\alias{MassArrayFragment-class}
\alias{$,MassArrayFragment-method}
\alias{$<-,MassArrayFragment-method}
\alias{initialize,MassArrayFragment-method}

\title{Class "MassArrayFragment" }
\description{A data structure containing information for a single fragment of an amplicon}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MassArrayFragment", ID, sequence, assay.name, name, position, type, direction, extra, bisulfite.converted, assayable, primer, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ID}:}{Unique integer indexing the fragment's position within a potential list of multiple fragments}
    \item{\code{assay.name}:}{(currently not supported)}
    \item{\code{name}:}{(currently not supported)}
    \item{\code{sequence}:}{Bisulphite converted nucleotide sequence of fragment}
    \item{\code{position}:}{Relative position of fragment within the amplicon}
    \item{\code{length}:}{Length (in bp) of fragment sequence}
    \item{\code{CpGs}:}{Number of CG dinucleotides contained within the fragment}
    \item{\code{MW}:}{Predicted molecular weight(s) of fragment, including methylated and unmethylated mass, adducts, etc.}
    \item{\code{collisions}:}{Number of fragments that share the same molecular weight as the current fragment}
    \item{\code{collision.IDs}:}{IDs of other fragments that share the same molecular weight as the current fragment}
    \item{\code{CG.collisions}:}{Number of CG-containing fragments that share the same molecular weight as the current fragment}
    \item{\code{CG.collision.IDs}:}{IDs of other CG-containing fragments that share the same molecular weight as the current fragment}
    \item{\code{type}:}{Specifies either 'T' or 'C' cleavage reaction}
    \item{\code{direction}:}{DNA strand used for fragment sequence (can be '+' or '-')}
    \item{\code{extra}:}{One of "5PPP-3P", "5OH", or "5PPP-3OH" (default)}
    \item{\code{bisulfite.converted}:}{Logical indicating whether the fragment sequence represents bisulfite converted sequence}
    \item{\code{assayable}:}{Logical indicating whether or not the fragment molecular weight is within the usable mass window}
    \item{\code{conversion.control}:}{Logical indicating whether or not the fragment is designated as a potential conversion control}
    \item{\code{required}:}{Logical indicating whether or not the fragment is designated as 'required' by the user}
    \item{\code{ignored}:}{Logical indicating whether or not the fragment is to be ignored}
    \item{\code{primer}:}{Logical indicating whether or not the fragment overlaps with primer or tagged sequence}
  }
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "MassArrayFragment")}: ... }
    \item{\$<-}{\code{signature(x = "MassArrayFragment")}: ... }
    \item{initialize}{\code{signature(.Object = "MassArrayFragment")}: ... }
	 }
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\seealso{
See Also \code{\linkS4class{MassArrayData}}
}
\examples{
showClass("MassArrayFragment")
}
\keyword{classes}
