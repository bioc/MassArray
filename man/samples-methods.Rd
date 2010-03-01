\name{samples-methods}
\docType{methods}
\alias{samples-methods}
\alias{samples,MassArrayData-method}
\alias{samples<-,MassArrayData,missing-method}
\alias{samples<-,MassArrayData,character-method}
\title{ Operate on sample name(s) (methods)}
\description{
 Methods to access (and/or assign) sample name information for a MassArrayData object
}
\section{Methods}{
\describe{

\item{object = "MassArrayData"}{ Access sample name information for MassArrayData object }

\item{object = "MassArrayData", value = "missing"}{ Handle empty function call, simply return the MassArrayData object }

\item{object = "MassArrayData", value = "character"}{ Assign sample name of MassArrayData object to \code{value} }
}}
\seealso{ \code{\link{samples}} }
\keyword{methods}
