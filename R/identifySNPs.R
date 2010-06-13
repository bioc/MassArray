identifySNPs <- function (peak.sequence, sequence, rxn=c("T", "C")) {
	rxn <- match.arg(rxn)
	if (is.null(peak.sequence)) return()
	if (length(peak.sequence) < 1) return()
	if (nchar(peak.sequence) < 1) return()
	sequence <- gsub("[?]", "", toupper(sequence))
	## DEFINE ACCESSORY FUNCTIONS
	## BASE COMPOSITION WILL ANALYZE A SEQUENCE AND GENERATE INDIVIDUAL NUCLEOTIDE COUNTS
	baseComposition <- function(sequence) {
		if (length(grep("R", sequence)) == 1) {
			return(c(baseComposition(sub("R", "A", sequence)), baseComposition(sub("R", "G", sequence))))
		}
		if (length(grep("Y", sequence)) == 1) {
			return(c(baseComposition(sub("Y", "C", sequence)), baseComposition(sub("Y", "T", sequence))))
		}
		if (length(grep("N", sequence)) == 1) {
			return(c(baseComposition(sub("N", "A", sequence)), baseComposition(sub("N", "T", sequence)), baseComposition(sub("N", "C", sequence)), baseComposition(sub("N", "G", sequence))))
		}
		base.counts <- list("A" = 0,
			"T" = 0,
			"C" = 0,
			"G" = 0
		)
		for (base in names(base.counts)) {
			base.count <- gregexpr(base, sequence)[[1]]
			if (base.count[1] > 0) base.counts[base] <- length(base.count)
		}
		return(base.counts)	
	}
	## LOCAL FRAGMENTATION WILL EXTRACT ALL POSSIBLE LOCAL FRAGMENTS AROUND THE "?" CHARACTER
	localFragmentation <- function(sequence, ignore) {
		switch(rxn,
			"T" = fragment <- gregexpr(paste("[T][^T]*[?][^T]*[T]|^[^T]*[?][^T]*[T]|[T][^T]*[?][^T]*$", sep=""), sequence)[[1]],
			"C" = fragment <- gregexpr(paste("[C][^C]*[?][^C]*[C]|^[^C]*[?][^C]*[C]|[C][^C]*[?][^C]*$", sep=""), sequence)[[1]]
		)
		if (fragment[1] < 0) return(FALSE)
		if (length(fragment) != 1) return(FALSE)
		sequence.sub <- substr(sequence, fragment, fragment + attr(fragment, "match.length") - 1)
		if (length(grep("R", sequence.sub)) == 1) {
			if (is.character(base <- localFragmentation(sub("R", "A", sequence.sub), ignore))) return(base)
			if (is.character(base <- localFragmentation(sub("R", "G", sequence.sub), ignore))) return(base)
		}
		if (length(grep("Y", sequence.sub)) == 1) {
			if (is.character(base <- localFragmentation(sub("Y", "C", sequence.sub), ignore))) return(base)
			if (is.character(base <- localFragmentation(sub("Y", "T", sequence.sub), ignore))) return(base)
		}
		if (length(grep("N", sequence.sub)) == 1) {
			if (is.character(base <- localFragmentation(sub("N", "A", sequence.sub), ignore))) return(base)
			if (is.character(base <- localFragmentation(sub("N", "T", sequence.sub), ignore))) return(base)
			if (is.character(base <- localFragmentation(sub("N", "C", sequence.sub), ignore))) return(base)
			if (is.character(base <- localFragmentation(sub("N", "G", sequence.sub), ignore))) return(base)
		}
		for (base in setdiff(c("A", "T", "C", "G", ""), ignore)) {
			sequence.base.sub <- sub("[?]", base, sequence.sub)
			switch(rxn,
				"T" = sub.fragments <- gregexpr("[ACGNXYR.]*T", sequence.base.sub)[[1]],
				"C" = sub.fragments <- gregexpr("[ATGNXYR.]*C", sequence.base.sub)[[1]]
			)
			sub.fragments <- substr(rep(sequence.base.sub, length(sub.fragments)), 
				sub.fragments, 
				as.integer(sub.fragments) + attr(sub.fragments, "match.length") - 1
			)
			sub.fragments <- lapply(sub.fragments, baseComposition)
			SNP.matched <- any(unlist(lapply(sub.fragments, identical, base.counts)))
			if (SNP.matched) return(as.character(base))
		}
		return(FALSE)
	}
	peak.sequence <- expandSequence(peak.sequence)
	## DETERMINE BASE COUNTS FOR PEAK SEQUENCE
	base.counts <- baseComposition(peak.sequence)
	
	N <- nchar(sequence)
	SNPs <- c()
	for (i in 1:N) {
		sequence.i <- paste(substr(sequence, 0, i - 1), "?", substr(sequence, i + 1, N), sep="")
		if (is.character(base <- localFragmentation(sequence.i, substr(sequence, i, i)))) {
			if (base == "") {
				SNPs <- c(SNPs, list(list(sequence=peak.sequence, position=i, base=base, type="deletion")))
			}
			else {
				SNPs <- c(SNPs, list(list(sequence=peak.sequence, position=i, base=base, type="substitution")))
			}
		}
	}
	return(SNPs)
}