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
	if (length(MW) < 1) return()
	for (i in 1:length(MW)) {
		peak.i <- which(abs(peaks - MW[i]) < resolution)
		if (length(peak.i) < 1) peak.i <- NA
		peaks.found <- c(peaks.found, peak.i)
	}
	return(peaks.found)
}

findFragments <- function (MW, fragments, resolution=1) {
	if (length(MW) < 1) return()
	peaks <- relist(!is.na(findPeaks(MW, unlist(lapply(fragments, slot, "MW")), resolution)), lapply(fragments, slot, "MW"))
	fragments.found <- which(unlist(lapply(peaks, all)))
	if (length(fragments.found) < 1) return()
	return(fragments[fragments.found])
}