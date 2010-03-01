\name{position-methods}
\docType{methods}
\alias{position-methods}
\alias{position,MassArrayData-method}
\alias{position<-,MassArrayData,missing-method}
\alias{position<-,MassArrayData,character-method}
\title{ Operate on positional information (methods)}
\description{
 Methods to access (and/or assign) positional information for a MassArrayData object
}
\section{Methods}{
\describe{

\item{object = "MassArrayData"}{ Access positional information for MassArrayData object }

\item{object = "MassArrayData", value = "missing"}{ Handle empty function call, simply return the MassArrayData object }

\item{object = "MassArrayData", value = "character"}{ Assign position of MassArrayData object to \code{value} }
}}
\seealso{ \code{\link{position}} }
\keyword{methods}
