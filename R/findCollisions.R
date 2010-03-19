findCollisions <- function (peaks, resolution=0.5) {
	collision <- function (MW, peaks, resolution) {
		peaks <- which(abs(peaks - MW) < resolution)
		if (length(peaks) <= 1) {
			return(0)	
		}
		else {
			return(peaks)
		}
	}
	return(lapply(peaks, collision, peaks, resolution))
}

numCollisions <- function (peaks, resolution=0.5) {
	collision <- function (MW, peaks, resolution) {
		peaks <- which(abs(peaks - MW) < resolution)
		return(length(peaks) - 1)
	}	
	return(unlist(lapply(peaks, collision, peaks, resolution)))
}

findPeaks <- function (MW, peaks, resolution=1) {
	peaks.found <- c()
	if ((length(MW) < 1) | (length(peaks) < 1)) return()
	for (i in 1:length(MW)) {
		peak.i <- which(abs(peaks - MW[i]) < resolution)
		if (length(peak.i) < 1) peak.i <- NA
		peaks.found <- c(peaks.found, peak.i)
	}
	return(peaks.found)
}

findFragments <- function (MW, fragments, resolution=1) {
	if (length(MW) < 1) return()
	MWs <- lapply(fragments, slot, "MW")
	MW.matches <- rep(NA, length(unlist(MWs)))
	if (length(MW.matches) < 1) return()
	MW.matches[findPeaks(MW, unlist(MWs), resolution)] <- 1
	MW.matches <- relist(!is.na(MW.matches), MWs)
	fragments.found <- which(unlist(lapply(MW.matches, any)))
	if (length(fragments.found) < 1) return()
	return(fragments[fragments.found])
}