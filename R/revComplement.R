setGeneric("revComplement", 
	function(x) {
		standardGeneric("revComplement")
	}
)

setMethod("revComplement", 
	signature = c("character"),
	function(x) {
		x <- toupper(x)
		x <- chartr("ATCGYR()<>", "TAGCRY)(><", x)
		x <- paste(rev(unlist(strsplit(x, ""))), sep="", collapse="")
		return(as.character(x))
	}
)

setMethod("revComplement", 
	signature = c("MassArrayData"),
	function(x) {
		if (!validObject(x)) {
			stop("input (x) is not a valid object of 'MassArrayData'")
		}
		x@sequence <- revComplement(x@sequence)
		N <- nchar(paste(x@fwd.tag, x@sequence, x@rev.tag, sep=""))
		x@CpG.data <- x@CpG.data[, dim(x@CpG.data)[2]:1]
		x@strand <- c("+", "-")[which(c("+", "-") != x@strand)]
		x@fragments.T <- rev(x@fragments.T)
		x@fragments.T <- lapply(x@fragments.T, 
			function(fragment) {
				fragment$position <- as.integer(N - (fragment$position + fragment$length - 1) + 1)
				return(fragment)
			}
		)
		x@fragments.C <- rev(x@fragments.C)
		x@fragments.C <- lapply(x@fragments.C, 
			function(fragment) {
				fragment$position <- as.integer(N - (fragment$position + fragment$length - 1) + 1)
				return(fragment)
			}
		)
		return(x)
	}
)