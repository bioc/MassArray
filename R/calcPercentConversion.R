calcPercentAdduct <- function(peaks) {
	sequences <- lapply(peaks, slot, "sequence")
	## IDENTIFY ADDUCTS AND THEIR MATCHING REFERENCE PEAKS
	adduct.peaks <- which((unlist(lapply(peaks, slot, "adduct")) != "") & (lapply(sequences, length) == 1))
	adduct.peaks <- unique(unlist(lapply(peaks[adduct.peaks], slot, "sequence")))
	## IDENTIFY MATCHING REFERENCE PEAKS
	adduct.ratios <- unlist(lapply(adduct.peaks,
		function (adduct) {
			## WHICH PEAKS MATCH THE SEQUENCE FOR THIS ADDUCT PEAK?
			matched.by.sequence <- which(unlist(lapply(lapply(sequences, "%in%", adduct), any)))
			## IDENTIFY ADDUCTS
			adducts <- matched.by.sequence[which(unlist(lapply(peaks[matched.by.sequence], slot, "adduct")) != "")]
			adducts <- sum(unlist(lapply(peaks[adducts], slot, "SNR")))
			## FILTER OUT ADDUCT, DOUBLY CHARGED, ABORTIVE CYCLING, AND OTHER MODIFIED PEAKS
			non.adducts <- matched.by.sequence[which(unlist(lapply(peaks[matched.by.sequence], slot, "type")) != "Modified")]
			non.adducts <- sum(unlist(lapply(peaks[non.adducts], slot, "SNR")))
			## CALCULATE RATIO OF ADDUCT PEAK HEIGHTS TO THEIR REFERENCES
			return(adducts / (adducts + non.adducts))
		}
	))

	return(median(adduct.ratios))
}


calcPercentConversion <- function(fragments, peaks) {
	controls <- which(unlist(lapply(fragments, slot, "conversion.control")))
	controls.meth <- c()
	for (i in controls) {
		controls.meth <- c(controls.meth, calcMeth(fragments[[i]]$MW, peaks))
	}
	return(as.numeric(controls.meth))
}


calcMeth <- function(peaks, peaklist, method=c("weighted", "proportion"), num.cg.fragments=rep(1, length(peaks))) {
	method <- match.arg(method)
	peaks <- peaks[order(peaks)]
	which.peaks.matched <- findPeaks(peaks, unlist(lapply(peaklist, slot, "MW.actual")))
	which.peaks.matched.theory <- findPeaks(peaks, unlist(lapply(peaklist, slot, "MW.theoretical")))
	which.peaks.matched[which(is.na(which.peaks.matched))] <- which.peaks.matched.theory[which(is.na(which.peaks.matched))]
	numerator <- 0
	denominator <- 0
	if (any(is.na(which.peaks.matched))) return(NA)
	N <- length(peaks)
	## CALCULATE NET PEAK INTENSITY FOR ALL FRAGMENTS BEING MEASURED
	SNR.sum <- sum(unlist(lapply(peaklist[which.peaks.matched], slot, "SNR")), na.rm=TRUE)
	for (i in 1:N) {
		SNR.i <- peaklist[[which.peaks.matched[i]]]$SNR
		## HOW MANY FRAGMENTS DETERMINE THE SIGNAL FOR THIS PEAK?
		num.fragments <- max(1, peaklist[[which.peaks.matched[i]]]$components)
		## SCALE SNR DATA TO ACCOUNT FOR PEAK COLLISIONS WITH CpG-CONTAINING AND/OR OTHER FRAGMENTS
		## NOTE: THE FOLLOWING SCALING ASSUMES THAT TOTAL PEAK INTENSITY IS THE SAME FOR ALL FRAGMENTS
		## NOTE: THE FORMULA REPRESENTS THE ALREADY SIMPLIFIED SOLUTION OF MULTIPLE LINEAR EQUATIONS
		## NOTE: EX... SNR1=A1+B1; SNR2=A2; A1+A2=Anet=Bnet=B1
		## NOTE: FOR NOW, THIS IS SUFFICIENT... BUT MAY EXTEND IN FUTURE TO ACCOUNT FOR ADD'L FRAGMENT DEGENERACY (WILL REQUIRE RECURSIVE FUNCTIONS)
		SNR.i <- SNR.i - SNR.sum * (num.fragments - num.cg.fragments[i]) / num.fragments
		## IF OTHER PEAKS OVERLAPPING => REMOVE FRACTION FROM CONSIDERATION!!!
		if (is.na(SNR.i)) next
		denominator <- denominator + SNR.i
		switch(method,
			"weighted" = numerator <- numerator + (i - 1) * SNR.i / (N - 1),
			"proportion" = numerator <- numerator + min(1, i - 1) * SNR.i
		)
	}
	if (denominator <= 0) return(NA)
	return(min(1, numerator / denominator))
}


analyzeCpGs <- function(fragments, peaks, method=c("weighted", "proportion")) {
	## HOW MANY CGs TO ANALYZE?
	CpG.num <- sum(unlist(lapply(fragments, slot, "CpGs"))[!unlist(lapply(fragments, slot, "conversion.control"))])
	## STORE CG DATA AS NUMERIC VECTOR
	CpG.data <- rep(NA, CpG.num)
	## IDENTIFY FRAGMENTS THAT HAVE CGs AND ARE NOT CONVERSION CONTROLS
	CpGs <- which((unlist(lapply(fragments, slot, "CpGs")) > 0) & (!unlist(lapply(fragments, slot, "conversion.control"))))
	CpG.fragment.map <- c()
	for (i in CpGs) {
		CpG.fragment.map <- c(CpG.fragment.map, rep(i, fragments[[i]]$CpGs))
	}
	for (i in 1:CpG.num) {
		## FRAGMENT COLLISIONS WITH CpG-CONTAINING FRAGMENTS
		collisions.i <- fragments[[CpG.fragment.map[i]]]$collision.IDs
		## COUNT NUMBER OF THESE THAT CONTAIN CpGs THEMSELVES
		cg.collision.counts <- unlist(lapply(collisions.i,
			function(x) {
				if (is.null(x) | identical(x, integer(0))) return(1)
				return(length(which((unlist(lapply(fragments[x], slot, "CpGs")) > 0)
					& !(unlist(lapply(fragments[x], slot, "conversion.control"))))) + 1)
			}
		))
		CpG.data[i] <- calcMeth(fragments[[CpG.fragment.map[i]]]$MW, peaks, method, cg.collision.counts)
	}
	return(CpG.data)
}
