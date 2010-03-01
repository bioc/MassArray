setClass("MassArraySpectrum",
	representation(
		sample = "character",
		rxn = "character",
		strand = "character",
		peaks = "list",
		quality.conversion = "numeric",
		quality.spectra = "numeric",
		quality.primerdimer = "numeric",
		quality.contaminant = "numeric",
		quality.adducts = "numeric"
	),
	prototype(
		sample = "",
		rxn = "",
		strand = "+",
		peaks = list(),
		quality.conversion = as.numeric(0),
		quality.spectra = as.numeric(0), 
		quality.primerdimer = as.numeric(0),
		quality.contaminant = as.numeric(0),
		quality.adducts = as.numeric(0)
	)
)


setMethod("initialize",
	"MassArraySpectrum",
	function(.Object,
		sample,
		rxn = c("T", "C", "CT", "TC"),
		strand = c("+", "-"),
		peaks = list(),
		quality.conversion = NA,
		quality.spectra = NA,
		quality.primerdimer = NA,
		quality.contaminant = NA,
		quality.adducts = NA,
		...	
	) {
		.Object@sample <- as.character(sample)
		.Object@rxn <- match.arg(rxn)
		.Object@strand <- match.arg(strand)
		.Object@peaks <- peaks
		.Object@quality.conversion <- as.numeric(quality.conversion)
		.Object@quality.spectra <- as.numeric(quality.spectra)
		.Object@quality.primerdimer <- as.numeric(quality.primerdimer)
		.Object@quality.contaminant <- as.numeric(quality.contaminant)
		.Object@quality.adducts <- as.numeric(quality.adducts)
		return(.Object)
	}
)


setValidity("MassArraySpectrum",
	function(object) {
		return(TRUE)
	}
)


setMethod("$", "MassArraySpectrum",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)


setMethod("$<-", "MassArraySpectrum",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'MassArraySpectrum'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)


