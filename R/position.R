setGeneric("position", 
	function(object) {
		standardGeneric("position")
	}
)

setMethod("position", 
	signature = c("MassArrayData"),
	function(object) {
		if (is.null(object@chr) | is.null(object@start) | is.null(object@end)) return("")
		if (identical(object@chr, character(0)) | identical(object@start, integer(0)) | identical(object@end, integer(0))) return("")
		if ((object@chr == "") | is.na(object@start) | is.na(object@end)) return("")
		return(paste("chr", object@chr, ":", object@start, "-", object@end, sep=""))
	}
)

setGeneric("position<-", 
	function(object, value) {
		standardGeneric("position<-")
	}
)

setMethod("position<-", 
	signature = c("MassArrayData", "missing"),
	function(object, value) {
		return(object)	
	}
)

setMethod("position<-", 
	signature = c("MassArrayData", "character"),
	function(object, value) {
		if (length(grep("^chr[1-9A-Z][0-9]?[:][1-9][0-9]*-[1-9]*[0-9]*$", value)) != 1) {
			if (value != "") {
				warning("positional information incorrect (must be of form chrXX:XXXX-XXXX)")
			}
			return(object)
		}
		fields <- unlist(strsplit(value, "[:-]"))
		if (is.na(fields[3])) {
			fields <- c(fields, as.numeric(fields[2]) + object@end - object@start)
		}
		if (fields[3] < fields[2]) {
			warning("positional information out of bounds (end < start)")
			return(object)
		}
		if (as.integer(fields[3]) - as.integer(fields[2]) + 1 != nchar(object@sequence)) {
			warning("positional information does not match sequence length")
			return(object)
		}
		object@chr <- sub("chr", "", fields[1])
		object@start <- as.integer(fields[2])
		object@end <- as.integer(fields[3])
		return(object)	
	}
)
