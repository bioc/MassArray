estimatePrimerDimer <- function (fragments, peaks, method=c("ratio", "mann-whitney")) {
	## FIRST IDENTIFY FRAGMENTS OVERLAPPING (OR NOT) PRIMERS
	is.primer <- unlist(lapply(fragments, slot, "primer"))
	fragments.primer <- which(is.primer)
	fragments.noprimer <- which(!is.primer)
	## THEN EXTRACT THE MOLECULAR WEIGHTS (ONLY CONSIDER MWs THAT ARE ASSAYABLE) . . .
	MWs.primer <- unlist(lapply(fragments[fragments.primer], slot, "MW"))
	MWs.noprimer <- unlist(lapply(fragments[fragments.noprimer], slot, "MW"))
	MWs.primer <- MWs.primer[isAssayable(MWs.primer)]
	MWs.noprimer <- MWs.noprimer[isAssayable(MWs.noprimer)]
	## . . . AND MATCH THE APPROPRIATE PEAKS
	peaks.primer <- findPeaks(MWs.primer, unlist(lapply(peaks, slot, "MW.theoretical")))
	peaks.noprimer <- findPeaks(MWs.noprimer, unlist(lapply(peaks, slot, "MW.theoretical")))
	## DO NOT CONSIDER MISSING OR OTHERWISE UNMATCHED PEAKS
	peaks.primer <- peaks.primer[which(!is.na(peaks.primer))]
	peaks.noprimer <- peaks.noprimer[which(!is.na(peaks.noprimer))]
	## CALCULATE SCALED SNR FOR EACH PEAK
	data.primer <- unlist(lapply(peaks[peaks.primer], slot, "SNR")) / unlist(lapply(peaks[peaks.primer], slot, "ref.intensity"))
	data.noprimer <- unlist(lapply(peaks[peaks.noprimer], slot, "SNR")) / unlist(lapply(peaks[peaks.noprimer], slot, "ref.intensity"))
#	SNR.p <- mean(unlist(lapply(peaks[peaks.primer], slot, "SNR")),na.rm=TRUE)
#	SNR.np <- mean(unlist(lapply(peaks[peaks.noprimer], slot, "SNR")),na.rm=TRUE)
#	R.p <- mean(unlist(lapply(peaks[peaks.primer], slot, "ref.intensity")),na.rm=TRUE)
#	R.np <- mean(unlist(lapply(peaks[peaks.noprimer], slot, "ref.intensity")),na.rm=TRUE)

	## IF INSUFFICIENT DATA TO CALCULATE PRIMER DIMER
	if ((length(data.primer) < 1) | (length(data.noprimer) < 1)) return(NA)

	method <- match.arg(method)
	## SIMPLE RATIO OF OVERALL PEAK HEIGHTS
	if (method == "ratio") {
		return(pmax(0, data.primer / median(data.noprimer, na.rm=TRUE) - 1))
	}
	## MANN-WHITNEY U-TEST:  DO POPULATIONS OF PEAKS FROM PRIMER AND NON-PRIMER HAVE DIFFERENT MEANS?
	else if (method == "mann-whitney") {
		test <- wilcox.test(data.primer, data.noprimer, paired=FALSE, conf.int=TRUE, exact=FALSE)
		return(test$p.value)
	}
}