\name{calcPercentAdduct}
\alias{calcPercentAdduct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate percent salt adducts }
\description{
  Function to calculate ratio of salt adduct peak heights to reference/unmodified peaks
}
\usage{
calcPercentAdduct(peaks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{ List of MassArrayPeak objects comprising complete spectral data }
}
\details{
  Salt adducts (either Na or K) are identified and compared to each of their unmodified reference peaks
}
\value{
  Returns a list of numerical values corresponding to the ratios of salt adduct peak heights to their unmodified reference peaks
}
\author{  Reid F. Thompson (\email{rthompso@aecom.yu.edu}), John M. Greally (\email{jgreally@aecom.yu.edu}) }
\examples{
data(MassArray.example.data)
adduct.ratios <- calcPercentAdduct(MassArray.example.data$samples[[1]]$peaks)
median(adduct.ratios)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ arith }
\keyword{ attribute }
