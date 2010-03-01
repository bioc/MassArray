isAssayable <- function (MW, lower.threshold=1500, upper.threshold=7000) {
	return(as.logical((MW >= lower.threshold & MW <= upper.threshold)))
}
