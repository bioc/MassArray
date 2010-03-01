rnaDigest <- function (sequence, type=c("T", "C")) {
	sequence <- toupper(sequence)
	sequence.length <- nchar(sequence)
	type <- match.arg(type)
	switch(type,
		"T" = fragments <- gregexpr("[ACGNXYR.]*T", sequence)[[1]],
		"C" = fragments <- gregexpr("C[ATGNXYR.]*", sequence)[[1]]
	)
	num.fragments <- length(fragments)
	positions <- as.integer(fragments)
	fragments <- substr(rep(sequence, num.fragments), as.integer(fragments), as.integer(fragments) + attr(fragments, "match.length") - 1)
	N <- nchar(paste(fragments, sep="", collapse=""))
	if (N < sequence.length) {
		if (type == "T") {
			fragments <- c(fragments, substr(sequence, N + 1, sequence.length))
			positions <- c(positions, N + 1)
		}
		else {
			fragments <- c(substr(sequence, 1, sequence.length - N), fragments)
			positions <- c(1, positions)
		}
		num.fragments <- num.fragments + 1
	}
	digests <- c()
	for (i in 1:num.fragments) {
		digests <- c(digests,
			new("MassArrayFragment", ID=i, sequence=fragments[i], position=positions[i], type=type, bisulfite.converted=TRUE)
		)
	}
	return(digests)
}
