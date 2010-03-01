\name{createWiggle}
\alias{createWiggle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create wiggle track }
\description{
 Function to create and write a wiggle track (UCSC Genome Browser format) to flat file from methylation data contained in a MassArrayData object
}
\usage{
createWiggle(x, file = "", append = FALSE, colors = NULL, na.rm = FALSE, sep = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \code{MassArrayData} object containing methylation data for at least one sample. }
  \item{file}{ location of file to write wiggle track information; if \option{""}, wiggle track prints to the standard output connection: see \code{\link{cat}}. }
  \item{append}{ logical; if \option{TRUE}, the output is appended to an existent wiggle track file. If \option{FALSE} (default), a new file with a new header is created and any existing file of the same name is destroyed. }
  \item{colors}{ vector of colors, indicates which colors to use for which wiggle track }
  \item{na.rm}{ logical; if \option{TRUE} (default), missing values are removed from data. If \option{FALSE} any missing values cause an error }
  \item{sep}{ a string used to separate columns. Using \option{sep = "\t"} (default) gives tab-delimited output. }
}
\references{ UCSC Genome Browser, \url{http://genome.ucsc.edu/goldenPath/help/customTrack.html}:
Kent, W.J., Sugnet, C. W., Furey, T. S., Roskin, K.M., Pringle, T. H., Zahler, A. M., and Haussler, D. The Human Genome Browser at UCSC. \emph{Genome Res}. \bold{12}(6), 996-1006 (2002). }
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
createWiggle(MassArray.example.data)
}
\seealso{ \code{\link{write}}, \code{\link{cat}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
