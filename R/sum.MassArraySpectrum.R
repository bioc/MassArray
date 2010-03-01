sum.MassArraySpectrum <- function (x, ..., trim=0, na.rm=TRUE) useMethod("sum")

sum.MassArraySpectrum <- function (x, ..., trim=0, na.rm=TRUE) {
	spectra <- c(x, ...)
	samples <- unlist(lapply(spectra, slot, "sample"))
	rxn <- unique(unlist(lapply(spectra, slot, "rxn")))
	if (length(rxn) > 1) {
		warning("Spectral data contains both T&C reactions")
		rxn <- "TC"
	}
	strand <- unlist(lapply(spectra, slot, "strand"))
#	quality.conversion <- unlist(lapply(lapply(spectra, slot, "quality.conversion"), mean, trim=trim, na.rm=na.rm))
	quality.conversion <- unlist(lapply(spectra, slot, "quality.conversion"))
	quality.spectra <- mean(unlist(lapply(spectra, slot, "quality.spectra")), trim=trim, na.rm=na.rm)
	peaks <- list()
	for (i in spectra) {
		peaks.i <- i$peaks
		MW <- unlist(lapply(peaks, "[", "MW"))
		## BECAUSE USING THEORETICAL MW VALUES, PEAKS SHOULD BE EXACT MATCHES => CAN EMPLOY match() FUNCTION INSTEAD OF findPeaks()
		matched.peaks <- match(unlist(lapply(peaks.i, slot, "MW.theoretical")), MW)
		defined.peaks <- which(!is.na(matched.peaks))
		for (j in defined.peaks) {
			peaks[[matched.peaks[j]]]$peaks <- c(peaks[[matched.peaks[j]]]$peaks, peaks.i[j])
		}
		additional.peaks <- which(is.na(matched.peaks))
		for (j in additional.peaks) {
			peaks <- c(peaks, list(list(MW=peaks.i[[j]]$MW.theoretical, peaks=peaks.i[j])))
		}
	}
	MW <- unlist(lapply(peaks, "[", "MW"))
	peaks <- lapply(peaks, "[", "peaks")
	for (i in 1:length(peaks)) {
		if (peaks[[i]]$peaks[[1]]$adduct == "") adduct <- NULL
		else adduct <- peaks[[i]]$peaks[[1]]$adduct
		MW.actual <- mean(unlist(lapply(peaks[[i]]$peaks, slot, "MW.actual")), trim=trim, na.rm=na.rm)
		if (is.na(MW.actual)) MW.actual <- MW[i]
		peaks[[i]] <- new("MassArrayPeak", ID=as.integer(i),
						MW.theoretical=MW[i],
						MW.actual=MW.actual,
						probability=mean(unlist(lapply(peaks[[i]]$peaks, slot, "probability")), trim=trim, na.rm=na.rm),
						SNR=sum(unlist(lapply(peaks[[i]]$peaks, slot, "SNR")), na.rm=na.rm),
						height=sum(unlist(lapply(peaks[[i]]$peaks, slot, "height")), na.rm=na.rm),
						intensity=sum(unlist(lapply(peaks[[i]]$peaks, slot, "intensity")), na.rm=na.rm),
						adduct=adduct, 
						type=peaks[[i]]$peaks[[1]]$type, 
						missing=all(unlist(lapply(peaks[[i]]$peaks, slot, "missing"))), 
						new=any(unlist(lapply(peaks[[i]]$peaks, slot, "new"))),
						collisions=sum(unlist(lapply(peaks[[i]]$peaks, slot, "collisions")), na.rm=na.rm),
						components=ceiling(mean(unlist(lapply(peaks[[i]]$peaks, slot, "components")), trim=trim, na.rm=na.rm))
					)
	}
	if (length(samples) <= 1) sample.name <- samples[1]
	else if (length(unique(samples)) <= 1) sample.name <- unique(samples)
	else sample.name <- paste("mean(", samples[1], " ... ", samples[length(samples)], ")", sep="")
	return(new("MassArraySpectrum", sample=sample.name , rxn=paste(rxn, collapse=""), strand=strand[1], peaks=peaks, quality.conversion=quality.conversion, quality.spectra=quality.spectra))
}
