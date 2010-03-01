\name{combine-methods}
\docType{methods}
\alias{combine-methods}
\alias{combine,MassArrayData,MassArrayData-method}
\alias{combine,MassArrayData,missing-method}
\title{ Combine MassArrayData objects (methods)}
\description{
Methods for joining two MassArrayData objects by sequence positions and samples, or simply operating on a single MassArrayData object to combine samples, depending on input
}
\section{Methods}{
\describe{

\item{x = "MassArrayData", y = "MassArrayData"}{ Combine two MassArrayData objects by position and then by sample }

\item{x = "MassArrayData", y = "missing"}{ Combine duplicate samples within the same MassArrayData object }
}}
\seealso{ See Also \code{\link{combine}} }
\keyword{methods}
\keyword{ manip}
