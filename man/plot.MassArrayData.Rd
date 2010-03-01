\name{plot.MassArrayData}
\alias{plot.MassArrayData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot MassArrayData }
\description{
  Function to generate graphical output for methylation data in a \code{MassArrayData} object
}
\usage{
\method{plot}{MassArrayData}(x, ..., collapse = TRUE, bars = TRUE, scale = TRUE, sequence = TRUE, labels = TRUE, colors = TRUE, main = position(x), width = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of class \code{MassArrayData} }
  \item{\dots}{ Other arguments to plot, currently not supported at this time }
  \item{collapse}{ Logical specifying whether or not to combine samples by unique \code{group} (see \code{\linkS4class{MassArrayData}}). If \code{TRUE}, each methylation values are averaged across all samples in each unique group. If \code{FALSE}, all samples are retained and plotted individually }
  \item{bars}{ Logical specifying whether or not to display error bars. If \code{TRUE} (and \code{collapse} is \code{TRUE}), the median absolute deviation is calculated for each group and plotted as an error bar for each methylation value. If \code{FALSE}, no error bars are displayed}
  \item{scale}{ Logical specifying whether or not to keep the x axis to scale.  If \code{TRUE}, methylation values are plotted as a function of relative position within the amplicon sequence. If \code{FALSE}, positional information is ignored and methylation values are evenly spaced across the plot window.}
  \item{sequence}{ Logical specifying whether or not to display the nucleotide sequence for the amplicon }
  \item{labels}{ Logical specifying whether or not to display data labels }
  \item{colors}{ Logical specifying whether or not to plot in color.  If \code{TRUE}, colors are used.  If \code{FALSE}, plotting occurs in black and white and grayscale. }
  \item{main}{ Label/title for overall plot (default is \code{""} }
  \item{width}{ Numerical value specifying the display width to use for each methylation value; number corresponds to the number of base pairs to include in both directions from the methylation position (default is \code{1.5})}
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
plot(MassArray.example.data,collapse=FALSE,bars=FALSE,scale=FALSE)
}
\seealso{ See Also \code{\linkS4class{MassArrayData}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
