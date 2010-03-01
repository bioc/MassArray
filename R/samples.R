setGeneric("samples", 
	function(object) {
		standardGeneric("samples")
	}
)

setMethod("samples", 
	signature = c("MassArrayData"),
	function(object) {
		return(as.character(unlist(lapply(object@samples, slot, "sample"))))
	}
)

setGeneric("samples<-", 
	function(object, value) {
		standardGeneric("samples<-")
	}
)

setMethod("samples<-", 
	signature = c("MassArrayData", "missing"),
	function(object, value) {
		return(object)	
	}
)

setMethod("samples<-", 
	signature = c("MassArrayData", "character"),
	function(object, value) {
		if (length(value) != length(object@samples)) {
			stop("number of samples does not match number of spectra in 'MassArrayData' object")
		}
		if (length(value) != dim(object@CpG.data)[1]) {
			stop("number of samples does not match dimensions of methylation data in 'MassArrayData' object")
		}
		for (i in 1:length(value)) {
			object@samples[[i]]$sample <- as.character(value[i])
		}
		rownames(object@CpG.data) <- value
		return(object)	
	}
)
