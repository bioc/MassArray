countCGs <- function (sequence) {
	sequence <- toupper(sequence)
	return(unlist(lapply(gregexpr("(YG|CR|CG)", sequence), 
		function (CGs) {
			if (is.na(CGs[1])) {
				return(NA)	
			}
			## NO MATCHES TO "CG"
			if (CGs[1] == -1) {
				return(0)	
			}
			## ONE OR MORE MATCHES TO "CG"
			else {
				return(length(CGs))	
			}
		}
	)))
}

