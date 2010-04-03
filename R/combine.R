setGeneric("combine", 
	function(x, y, ...) {
		standardGeneric("combine")
	}
)

setMethod("combine",
	signature = c("MassArrayData", "missing"),
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("input (x) is not a valid object of 'MassArrayData'")
		}
		sample.names <- samples(x)
		sample.names.unique <- unique(sample.names)
		N <- dim(x@CpG.data)[2]
		x@CpG.data.combined <- matrix(NA, nrow=length(sample.names.unique), ncol=N, dimnames=list(sample.names.unique, 1:N))
		for (i in 1:length(sample.names.unique)) {
			j <- which(sample.names == sample.names.unique[i])
			if (length(j) <= 1) {
				x@CpG.data.combined[i, ] <- x@CpG.data[j, ]
			}
			else {
				x@CpG.data.combined[i, ] <- apply(x@CpG.data[j, ], 2, mean, na.rm=TRUE)
			}
		}
		return(x)
	}
)


setMethod("combine",
	signature = c("MassArrayData", "MassArrayData"),
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("input (x) is not a valid object of 'MassArrayData'")
		}
		if (!validObject(y)) {
			stop("input (y) is not a valid object of 'MassArrayData'")
		}
		if (!identical(x@chr, y@chr)) {
			stop("sequences are non-adjacent (chromosome mismatch)")
		}
		if ((length(x@start) < 1) | (length(x@end) < 1) | (length(y@start) < 1) | (length(y@end) < 1)) {
			stop("positional information not specified -- see position() function")
		}
		original <- x@strand
		## NEED TO UPDATE THIS STRAND SPECIFICITY FOR MASSARRAYDATA OBJECTS!!!
		if (identical(x@strand, "-")) {
			x <- revComplement(x)
		}
		if (identical(y@strand, "-")) {
			y <- revComplement(y)
		}
		bounds <- c(max(x@start, y@start), min(x@end, y@end))
		overlap <- (bounds[1] <= bounds[2])
		bounds <- bounds[order(bounds)]
		x.cgs <- c(gregexpr("(CG|YG|CR)", x@sequence)[[1]] + x@start - 1, bounds)
		y.cgs <- c(gregexpr("(CG|YG|CR)", y@sequence)[[1]] + y@start - 1, bounds)
		cgs <- unique(c(x.cgs, y.cgs, bounds))
		x.cgs <- x.cgs[order(x.cgs)]
		y.cgs <- y.cgs[order(y.cgs)]
		cgs <- cgs[order(cgs)]
		x.A <- which(x.cgs == bounds[1])[1]
		x.B <- which(x.cgs == bounds[2])[1]
		y.A <- which(y.cgs == bounds[1])[1]
		y.B <- which(y.cgs == bounds[2])[1]
		A <- which(cgs == bounds[1])[1]
		if ((x.B - x.A) != (y.B - y.A)) {
			stop("sequences are non-adjacent (overlap CG mismatch)")
		}
		## UPDATE FRAGMENTATION PROFILE POSITIONS TO BE RELATIVE TO NEWLY COMBINED SEQUENCE INFORMATION
		if (x@start < y@start) {
#			num.preceding.fragments <- length(which(unlist(lapply(x@fragments.T, slot, "position")) < y@start))
			for (i in 1:length(y@fragments.T)) {
				y@fragments.T[[i]]$assay.name <- as.character(position(y))
				y@fragments.T[[i]]$position <- y@fragments.T[[i]]$position + y@start - x@start
				y@fragments.T[[i]]$ID <- as.integer(y@fragments.T[[i]]$ID + length(x@fragments.T))
				y@fragments.T[[i]]$collision.IDs <- relist(unlist(y@fragments.T[[i]]$collision.IDs) + length(x@fragments.T), y@fragments.T[[i]]$collision.IDs)
			}				
			for (i in 1:length(y@fragments.C)) {
				y@fragments.C[[i]]$assay.name <- as.character(position(y))
				y@fragments.C[[i]]$position <- y@fragments.C[[i]]$position + y@start - x@start
				y@fragments.C[[i]]$ID <- as.integer(y@fragments.C[[i]]$ID + length(x@fragments.C))
				y@fragments.C[[i]]$collision.IDs <- relist(unlist(y@fragments.C[[i]]$collision.IDs) + length(x@fragments.C), y@fragments.C[[i]]$collision.IDs)
			}
		}
		else {
			for (i in 1:length(x@fragments.T)) {
				x@fragments.T[[i]]$assay.name <- as.character(position(x))
				x@fragments.T[[i]]$position <- x@fragments.T[[i]]$position + x@start - y@start
			}				
			for (i in 1:length(x@fragments.C)) {
				x@fragments.C[[i]]$assay.name <- as.character(position(x))
				x@fragments.C[[i]]$position <- x@fragments.C[[i]]$position + x@start - y@start
			}
			for (i in 1:length(y@fragments.T)) {
				y@fragments.T[[i]]$assay.name <- as.character(position(y))
				y@fragments.T[[i]]$ID <- as.integer(y@fragments.T[[i]]$ID + length(x@fragments.T))
				y@fragments.T[[i]]$collision.IDs <- relist(unlist(y@fragments.T[[i]]$collision.IDs) + length(x@fragments.T), y@fragments.T[[i]]$collision.IDs)
			}				
			for (i in 1:length(y@fragments.C)) {
				y@fragments.C[[i]]$assay.name <- as.character(position(y))
				y@fragments.C[[i]]$ID <- as.integer(y@fragments.C[[i]]$ID + length(x@fragments.C))
				y@fragments.C[[i]]$collision.IDs <- relist(unlist(y@fragments.C[[i]]$collision.IDs) + length(x@fragments.C), y@fragments.C[[i]]$collision.IDs)
			}
		}
		x@fragments.T <- c(x@fragments.T, y@fragments.T)
		x@fragments.C <- c(x@fragments.C, y@fragments.C)
		if (overlap) {
			if (substr(x@sequence, bounds[1] - x@start + 1, bounds[2] - x@start + 1) != substr(y@sequence, bounds[1] - y@start + 1, bounds[2] - y@start + 1)) {
				stop("sequences are non-adjacent (overlap sequence mismatch)")
			}
			x@sequence <- paste(substr(x@sequence, 0, bounds[1] - x@start), 
					substr(y@sequence, 0, bounds[1]-y@start), 
					substr(x@sequence, bounds[1] - x@start + 1, bounds[2] - x@start + 1), 
					substr(x@sequence, bounds[2] - x@start + 2, nchar(x@sequence) + 1), 
					substr(y@sequence, bounds[2] - y@start + 2, nchar(y@sequence) + 1), 
					sep=""
				)
		}
		else {
			x@sequence <- paste(substr(x@sequence, 0, bounds[1] - x@start + 1), 
					substr(y@sequence, 0, bounds[1] - y@start + 1), 
					paste(rep(".", bounds[2] - bounds[1] - 1), collapse=""), 
					substr(x@sequence, bounds[2] - x@start + 1, nchar(x@sequence) + 1), 
					substr(y@sequence, bounds[2] - y@start + 1, nchar(y@sequence) + 1), 
					sep=""
				)
		}
		overlap <- x.B-x.A-1
		x.pre <- x.A-1
		y.pre <- y.A-1
		x.post <- length(x.cgs) - x.B
		y.post <- length(y.cgs) - y.B
		## COMBINE CpG METHYLATION DATA MATRICES
		N <- as.integer(x.pre + y.pre + overlap + x.post + y.post)
		CpG.data <- matrix(NA, nrow=length(c(x@samples, y@samples)), ncol=N, dimnames=list(c(samples(x), samples(y)), 1:N))
		if (length(x@samples) > 0) {
			CpG.data[1:length(x@samples), (A-x.A+1):(A-x.A+length(x.cgs)-2)] <- as.matrix(x@CpG.data)
		}
		if (length(y@samples) > 0) {
			CpG.data[1:length(y@samples)+length(x@samples), (A-y.A+1):(A-y.A+length(y.cgs)-2)] <- as.matrix(y@CpG.data)
		}
		## COMBINE SPECTRA
		x@samples <- c(x@samples, y@samples)
		x@groups <- c(x@groups, y@groups)
		samples.combined <- groups <- c()
		sample.names <- samples(x)
		sample.names.unique <- unique(sample.names)
		x@CpG.data <- CpG.data
		x@CpG.data.combined <- matrix(NA, nrow=length(sample.names.unique), ncol=N, dimnames=list(sample.names.unique, 1:N))
		for (i in 1:length(sample.names.unique)) {
			j <- which(sample.names == sample.names.unique[i])
			if (length(j) <= 1) {
				x@CpG.data.combined[i, ] <- CpG.data[j, ]
			}
			else {
				x@CpG.data.combined[i, ] <- apply(CpG.data[j, ], 2, mean, na.rm=TRUE)
			}
			if (length(unique(x@groups[j])) == 1) {
				groups <- c(groups, unique(x@groups[j]))
			}
			else {
				groups <- c(groups, x@groups[j][1])
				warning("multiple groups specified for a single sample, only the first will be used")
			}
			samples.combined <- c(samples.combined, list(sum.MassArraySpectrum(x@samples[j])))
		}
		x@CpG.data <- x@CpG.data.combined
		x@samples <- samples.combined
		x@groups <- as.character(groups)
		rownames(x@CpG.data.combined) <- sample.names.unique
		x@start <- min(x@start, y@start)
		x@end <- max(x@end, y@end)
		if (original == "-") {
			x <- revComplement(x)
		}
		
		return(x)
		
		
#		x@data <- matrix(NA, nrow=length(x@samples), ncol=x@N, dimnames=list(x@samples, 1:x@N))
#		conversion <- c(x@conversion, y@conversion)
#		x@conversion <- as.numeric(rep(NA, length(x@samples)))
#		for (i in 1:length(x@samples)) {
#			x@conversion[i] <- mean(conversion[j], na.rm=TRUE)
#			x@rxn[i] <- rxn[j][1]
#		}
	}
)
