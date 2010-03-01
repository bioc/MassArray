bisConvert <- function (sequence) {
	sequence <- toupper(sequence)
	sequence <- gsub("CG", "YG", sequence)
	sequence <- gsub("C", "T", sequence)
	return(sequence)
}
