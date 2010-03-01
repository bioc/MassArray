\name{MassArray.example.data}
\Rdversion{1.1}
\alias{MassArray.example.data}
\docType{data}
\title{
MassArray Data object
}
\description{
This data contains MassArray spectral information for two samples.
}
\usage{MassArray.example.data}
\format{
  The format is:
Formal class 'MassArrayData' [package "MassArray"] with 17 slots
  ..@ sequence         : chr "CCAGGTCCAAAGGTTCAGACCAGTCTGAACCTGTCCAGGGGCACTCCATATTTTCCTACCTGTCCCTCTTTGCTTGTAAAAACAAATTAAACAGGGATCCCAGCAACTTCGGGGCATGTGTGTAACT"| __truncated__
  ..@ chr              : chr(0) 
  ..@ start            : int(0) 
  ..@ end              : int(0) 
  ..@ strand           : chr "+"
  ..@ fwd.tag          : chr "AGGAAGAGAG"
  ..@ rev.tag          : chr "AGCCTTCTCCC"
  ..@ fwd.primer       : num 29
  ..@ rev.primer       : num 27
  ..@ lower.threshold  : num 1500
  ..@ upper.threshold  : num 9000
  ..@ fragments.T      :List of 89
  .. ..$ :Formal class 'MassArrayFragment' [package "MassArray"] with 21 slots
  .. .. .. ..@ ID                 : int 1
  .. .. .. ..@ assay.name         : chr ""
  .. .. .. ..@ name               : chr ""
  .. .. .. ..@ sequence           : chr "GGGAGAAGGCT"
  .. .. .. ..@ position           : int 385
  .. .. .. ..@ length             : int 11
  .. .. .. ..@ CpGs               : int 0
  .. .. .. ..@ MW                 : num 3913
  .. .. .. ..@ collisions         : int 0
  .. .. .. ..@ collision.IDs      :List of 1
  .. .. .. .. ..$ : int(0) 
  .. .. .. ..@ CG.collisions      : int 0
  .. .. .. ..@ CG.collision.IDs   : list()
  .. .. .. ..@ type               : chr "T"
  .. .. .. ..@ direction          : chr "+"
  .. .. .. ..@ extra              : chr "5PPP-3P"
  .. .. .. ..@ bisulfite.converted: logi TRUE
  .. .. .. ..@ assayable          : logi TRUE
  .. .. .. ..@ conversion.control : logi FALSE
  .. .. .. ..@ required           : logi FALSE
  .. .. .. ..@ ignored            : logi FALSE
  .. .. .. ..@ primer             : logi TRUE
  ..@ samples          :List of 2
  .. ..$ :Formal class 'MassArraySpectrum' [package "MassArray"] with 9 slots
  .. .. .. ..@ sample             : chr "A"
  .. .. .. ..@ rxn                : chr "T"
  .. .. .. ..@ strand             : chr "+"
  .. .. .. ..@ peaks              :List of 184
  .. .. .. .. ..$ :Formal class 'MassArrayPeak' [package "MassArray"] with 16 slots
  .. .. .. .. .. .. ..@ ID              : int 1
  .. .. .. .. .. .. ..@ MW.theoretical  : num 1111
  .. .. .. .. .. .. ..@ MW.actual       : num NA
  .. .. .. .. .. .. ..@ probability     : num 0
  .. .. .. .. .. .. ..@ SNR             : num 0
  .. .. .. .. .. .. ..@ height          : num NA
  .. .. .. .. .. .. ..@ sample.intensity: num NA
  .. .. .. .. .. .. ..@ ref.intensity   : num 0.1
  .. .. .. .. .. .. ..@ sequence        : chr "ACACAAT"
  .. .. .. .. .. .. ..@ adduct          : chr ""
  .. .. .. .. .. .. ..@ type            : chr "Modified"
  .. .. .. .. .. .. ..@ charge          : int 1
  .. .. .. .. .. .. ..@ collisions      : int 0
  .. .. .. .. .. .. ..@ components      : int 0
  .. .. .. .. .. .. ..@ missing         : logi TRUE
  .. .. .. .. .. .. ..@ new             : logi FALSE
  .. .. .. ..@ quality.conversion : num [1:4] 0.0529 0 0 0
  .. .. .. ..@ quality.spectra    : num NA
  .. .. .. ..@ quality.primerdimer: num [1:7] 8.51 15.83 3.28 1.04 0 ...
  .. .. .. ..@ quality.contaminant: num NA
  .. .. .. ..@ quality.adducts    : num [1:114] 1 1 0.97 0.231 0.412 ...
  ..@ groups           : chr(0) 
  ..@ CpG.data         : num [1:2, 1:18] 0.0322 0.0449 0.1468 0.3641 0.1468 ...
  ..@ CpG.data.combined: num [1:2, 1:18] 0.0322 0.0449 0.1468 0.3641 0.1468 ...
}
\source{
Thompson et al. 2009
}
\examples{
data(MassArray.example.data)
}
\keyword{datasets}
