setClass("MassArrayData",
	representation(
		sequence = "character",
		chr = "character",
		start = "integer",
		end = "integer",
		strand = "character",
		fwd.tag = "character",
		rev.tag = "character",
		fwd.primer = "numeric",
		rev.primer = "numeric",
		lower.threshold = "numeric", 
		upper.threshold = "numeric", 
		fragments.T = "list",
		fragments.C = "list",
		samples = "list",
		groups = "character",
		CpG.data = "matrix",
		CpG.data.combined = "matrix"
	),
	prototype(
		sequence = character(),
		chr = character(),
		start = integer(),
		end = integer(),
		fwd.tag = character(),
		rev.tag = character(),
		fwd.primer = numeric(),
		rev.primer = numeric(),
		lower.threshold = numeric(),
		upper.threshold = numeric(),
		fragments.T = list(),
		fragments.C = list(),
		samples = list(),
		groups = character(),
		CpG.data = matrix(),
		CpG.data.combined = matrix()
	)
)


setMethod("initialize",
	"MassArrayData",
	function (.Object,
		sequence, 
		file = "", 
		verbose=TRUE, 
		fwd.tag="AGGAAGAGAG", 
		rev.tag="AGCCTTCTCCC", 
		fwd.primer=0, 
		rev.primer=0, 
		strand=c("+", "-"), 
		lower.threshold=1500, 
		upper.threshold=9000, 
		header=TRUE, 
		skip=1, 
		sep="\t", 
		comment.char="#", 
		fill=TRUE, 
		method=c("weighted", "proportion"),
		position="",
		multiple.conversion=FALSE,
		...
	) {	
		.Object@strand <- match.arg(strand)
		method <- match.arg(method)
		.Object@sequence <- gsub("[^ACGTNXYR.()><]", "", toupper(sequence))
		.Object@fwd.tag <- gsub("[^ACTG]", "", toupper(fwd.tag))
		.Object@rev.tag <- gsub("[^ACTG]", "", toupper(rev.tag))
		.Object@lower.threshold <- as.numeric(lower.threshold)
		.Object@upper.threshold <- as.numeric(upper.threshold)
		switch(.Object@strand,
			## STANDARD REV-T7 RXN
			"+" = sequence <- .Object@sequence,
			## PERFORM ASSAY ON OPPOSITE STRAND, USE FWD-T7
			"-" = sequence <- revComplement(.Object@sequence)
		)
		required <- gregexpr("[(][ACGTNNXYR.]+[)]", .Object@sequence)[[1]]
		.Object@sequence <- gsub("[()]", "", .Object@sequence)
		attr(required, "match.length") <- attr(required, "match.length") - 2
		required <- required + nchar(.Object@fwd.tag)
		required <- required - 2 * (1:length(required) - 1)
		for (primer.flags in gregexpr("[<>]", .Object@sequence)[[1]]) {
			adjust.length <- (required <= primer.flags & required + attr(required, "match.length") - 1 >= primer.flags)
			attr(required[adjust.length], "match.length") <- attr(required[adjust.length], "match.length") - 1
			adjust.position <- (required >= primer.flags)
			required[adjust.position] <- required[adjust.position] - 1
		}
		if (fwd.primer <= 0) {
			.Object@fwd.primer <- max(0, attr(regexpr("^[^>]*(?=[>])", .Object@sequence, perl=TRUE), "match.length"))
		}
		else {
			.Object@fwd.primer <- fwd.primer	
		}
		if (.Object@fwd.primer <= 0) {
			warning("FWD primer not specified! Please consult the MassArray package manual for details on specifying primers")
		}
		if (rev.primer <= 0) {
			.Object@rev.primer <- max(0, attr(regexpr("(?<=[<])[^<]*$", .Object@sequence, perl=TRUE), "match.length"))
		}
		else {
			.Object@rev.primer <- rev.primer	
		}
		if (.Object@rev.primer <= 0) {
			warning("REV primer not specified! Please consult the MassArray package manual for details on specifying primers")
		}
		.Object@sequence <- gsub("[<>]", "", .Object@sequence)
		position(.Object) <- position
		if (verbose) cat("Performing inSilico Fragmentation: T, ")
		.Object@fragments.T <- inSilicoFragmentation(.Object@sequence, fwd.tag, rev.tag, type="T", .Object@lower.threshold, .Object@upper.threshold, .Object@fwd.primer, .Object@rev.primer, multiple.conversion)
		fragment.starts <- unlist(lapply(.Object@fragments.T, slot, "position"))
		fragment.ends <- fragment.starts + unlist(lapply(.Object@fragments.T, slot, "length")) - 1
		for (i in required) {
			fragments.i <- which((fragment.ends >= i) & (fragment.starts <= i + attr(required, "match.length")[which(required == i)] - 1))
			for (j in fragments.i) {
				.Object@fragments.T[[j]]$required <- TRUE
			}
		}
		if (verbose) cat("C ... ")
		.Object@fragments.C <- inSilicoFragmentation(.Object@sequence, .Object@fwd.tag, .Object@rev.tag, type="C", .Object@lower.threshold, .Object@upper.threshold, .Object@fwd.primer, .Object@rev.primer, multiple.conversion)
		fragment.starts <- unlist(lapply(.Object@fragments.C, slot, "position"))
		fragment.ends <- fragment.starts + unlist(lapply(.Object@fragments.C, slot, "length")) - 1
		for (i in required) {
			fragments.i <- which((fragment.ends >= i) & (fragment.starts <= i + attr(required, "match.length")[which(required == i)] - 1))
			for (j in fragments.i) {
				.Object@fragments.C[[j]]$required <- TRUE
			}
		}

		if (verbose) cat("FINISHED\n")
		if (!file.exists(file)) {
			if (verbose) warning("Cannot open data file:", file)
			CpG.num <- gregexpr("CG", paste(.Object@fwd.tag, .Object@sequence, .Object@rev.tag, sep=""))[[1]]
			if (CpG.num[1] < 0) {
				CpG.num <- 0	
			}
			else {
				CpG.num <- length(CpG.num)
			}
			.Object@CpG.data <- .Object@CpG.data.combined <- matrix(NA, nrow=length(.Object@samples), ncol=CpG.num)
			return(.Object)
		}
		if (verbose) cat("Importing matched peaks file (", file, "):\n", sep="")
	    data <- read.table(file, header=header, skip=skip, sep=sep, comment.char=comment.char, fill=fill)
		## EPITYPER DATA (MATCHED PEAKS) CAN COME IN TWO FORMATS (OLD & NEW) => HANDLE DIFFERENTLY
	    if ("Amplicon" %in% colnames(data)) {
			.Object@samples <- importEpiTyperData.new(data, .Object, verbose)
	    }
	    else {
			.Object@samples <- importEpiTyperData(data, .Object, verbose)
	    }
		if (length(.Object@samples) < 1) return(.Object)
		if (verbose) cat("\tAnalyzing conversion control(s) ... ")
		## CALCULATE CONVERSION CONTROLS FOR EACH SAMPLE
		for (i in 1:length(.Object@samples)) {
			switch(.Object@samples[[i]]$rxn,
				"T" = {
					.Object@samples[[i]]@quality.conversion <- calcPercentConversion(.Object@fragments.T, .Object@samples[[i]]@peaks)
				},
				"C" = {
					.Object@samples[[i]]@quality.conversion <- calcPercentConversion(.Object@fragments.C, .Object@samples[[i]]@peaks)
				}
			)
		}
		## ESTIMATE PRIMER DIMER LEVELS FOR EACH SAMPLE
		if (verbose) cat("FINISHED\n\tEstimating primer dimer level(s) ... ")
		if (.Object@fwd.primer * .Object@rev.primer > 0) {
			for (i in 1:length(.Object@samples)) {
				sample.i <- .Object@samples[[i]]
				switch(sample.i$rxn,
					"T" = primerdimer <- estimatePrimerDimer(.Object@fragments.T, sample.i$peaks),
					"C" = primerdimer <- estimatePrimerDimer(.Object@fragments.C, sample.i$peaks)
				)
				.Object@samples[[i]]$quality.primerdimer <- primerdimer
			}
			if (verbose) cat("FINISHED\n\tEstimating adduct level(s) ... ")
		}
		else {
			if (verbose) cat("ERROR (no primers specified)\n\tEstimating adduct level(s) ... ")	
		}
		## ESTIMATE PREVALENCE OF ADDUCTS IN EACH SAMPLE
		for (i in 1:length(.Object@samples)) {
			sample.i <- .Object@samples[[i]]
			switch(sample.i$rxn,
				   "T" = adduct <- calcPercentAdduct(sample.i$peaks),
				   "C" = adduct <- calcPercentAdduct(sample.i$peaks)
				   )
			.Object@samples[[i]]$quality.adducts <- adduct
		}
		## CALCULATE CpG METHYLATION CALLS FOR CPGS
		if (verbose) cat("FINISHED\n\tAnalyzing CpG methylation ...")
		CpG.num <- gregexpr("CG", paste(.Object@fwd.tag, .Object@sequence, .Object@rev.tag, sep=""))[[1]]
		if (CpG.num[1] < 0) {
			CpG.num <- 0	
		}
		else {
			CpG.num <- length(CpG.num)
		}
		.Object@CpG.data <- matrix(NA, nrow=length(.Object@samples), ncol=CpG.num)
		if (length(.Object@samples) > 0) {
			for (i in 1:length(.Object@samples)) {
				if (verbose) cat(".")
				sample.i <- .Object@samples[[i]]
				switch(sample.i$rxn,
					"T" = CpG.data.i <- rev(analyzeCpGs(.Object@fragments.T, sample.i$peaks, method)),
					"C" = CpG.data.i <- rev(analyzeCpGs(.Object@fragments.C, sample.i$peaks, method))
				)
				.Object@CpG.data[i, ] <- CpG.data.i
			}
		}
		.Object@CpG.data.combined <- .Object@CpG.data

		if (verbose) cat(" FINISHED\n")
		return(.Object)
	}
)


setValidity("MassArrayData",
	function(object) {
		## TEST WHETHER CpG COUNT MATCHES IN BOTH T&C RXN
		CpG.num <- length(gregexpr("CG", paste(object@fwd.tag, object@sequence, object@rev.tag, sep=""))[[1]])
		CpGs.T <- which(unlist(lapply(object@fragments.T, slot, "CpGs")) > 0 & !unlist(lapply(object@fragments.T, slot, "conversion.control")))
		CpGs.T <- CpGs.T[which(!duplicated(unlist(lapply(object@fragments.T[CpGs.T], slot, "position"))))]
		CpG.num.T <- sum(unlist(lapply(object@fragments.T[CpGs.T], slot, "CpGs")))
		if (!identical(CpG.num, CpG.num.T)) return(FALSE)
		CpGs.C <- which(unlist(lapply(object@fragments.C, slot, "CpGs")) > 0 & !unlist(lapply(object@fragments.C, slot, "conversion.control")))
		CpGs.C <- CpGs.C[which(!duplicated(unlist(lapply(object@fragments.C[CpGs.C], slot, "position"))))]
		CpG.num.C <- sum(unlist(lapply(object@fragments.C[CpGs.C], slot, "CpGs")))
		if (!identical(CpG.num, CpG.num.C)) return(FALSE)
		## TEST WHETHER SAMPLE COUNT MATCHES CpG METHYLATION DATA
		if (!identical(length(object@samples), dim(object@CpG.data)[1])) return(FALSE)
		## TEST IF POSITIONAL INFORMATION CONTAINS VALID BOUNDARIES (I.E. START <= END)
		if ((length(object@start) > 0) & (length(object@end > 0))) {
			if (object@start > object@end) return(FALSE)
		}
		return(TRUE)
	}
)

setMethod("$", "MassArrayData",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)


setMethod("$<-", "MassArrayData",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'MassArrayData'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)


setMethod("[", "MassArrayData",
	function (x, i, j, ..., drop=NA) {
		if (!validObject(x)) {
			stop("not a valid object of 'MassArrayData'")
		}
		if (!missing("i")) {
			if (is.logical(i)) {
				if (length(i) != length(x@samples)) {
					stop("(subscript 'i') logical subscript does not match dimensions of MassArrayData object")
				}
				i <- which(i)
			}
			if (is.character(i)) {
				matched <- match(samples(x), i)
				i <- which(!is.na(matched))
			}
			i <- as.integer(i)
			if (is.integer(i)) {
				if (length(i) < 1) {
					x@samples <- list()
					x@CpG.data <- x@CpG.data.combined <- matrix(nrow=0, ncol=dim(x@CpG.data)[2])
					x@groups <- character()
					return(x)	
				}
				if (max(i) > length(x@samples) | min(i) < 1) {
					stop("(subscript 'i') subscript out of bounds")
				}
				samples <- list()
				for (sample.i in i) {
					samples <- c(samples, x@samples[[sample.i]])
				}
				x@samples <- samples
				x@CpG.data <- matrix(x@CpG.data[i, ], ncol=dim(x@CpG.data)[2], dimnames=list(samples(x), 1:dim(x@CpG.data)[2]))
				x@CpG.data.combined <- matrix(x@CpG.data.combined[which(!is.na(match(unique(samples(x)), samples(x)))), ], ncol=dim(x@CpG.data)[2], dimnames=list(unique(samples(x)), 1:dim(x@CpG.data)[2]))
				x@groups <- x@groups[i]
			}
		}
		if (!missing("j")) {
			if (is.character(j)) {
				stop("(subscript 'j') subscript must contain numeric or logical input")
			}
			if (is.logical(j)) {
				if (length(j) != dim(x@CpG.data)[2]) {
					stop("(subscript 'j') logical subscript does not match dimensions of MassArrayData object")
				}
				j <- which(j)
			}
			j <- as.integer(j)
			if (is.integer(j)) {
				if (max(j) > dim(x@CpG.data)[2] | min(j) < 1) {
					stop("(subscript 'j') subscript out of bounds")
				}
				x@CpG.data <- matrix(x@CpG.data[, j], ncol=length(j), dimnames=list(samples(x), j))
				x@CpG.data.combined <- matrix(x@CpG.data.combined[, j], ncol=length(j), dimnames=list(samples(x), j))
				## REMOVE CGs FROM SEQUENCE AND FRAGMENTATION LISTS
				CG.positions <- unlist(gregexpr("(CG|YG|CR)", x@sequence))
				for (CG in CG.positions[setdiff(1:length(CG.positions), j)]) {
					substr(x@sequence, CG, CG) <- "T"
				}
				x@fragments.T <- inSilicoFragmentation(x@sequence, x@fwd.tag, x@rev.tag, type="T", x@lower.threshold, x@upper.threshold, x@fwd.primer, x@rev.primer)
				x@fragments.C <- inSilicoFragmentation(x@sequence, x@fwd.tag, x@rev.tag, type="C", x@lower.threshold, x@upper.threshold, x@fwd.primer, x@rev.primer)
			}
		}
		return(x)
	}
)

