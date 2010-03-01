expandSequence <- function (sequence) {
	sequence <- toupper(gsub(" ", "", sequence))
	## SEQUENCE INPUT CAN CONTAIN BASE SHORT HAND (I.E. "A3T2G1" IS INTERPRETED AS "AAATTG")
	base.repeats <- gregexpr("[ATCG][0-9]+", sequence)
	sequence <- unlist(mapply(
		function (bases, sequence) {
			if (bases[1] < 0) return(sequence)
			sequence.old <- sequence
			for (base in bases) {
				sequence <- sub("[ATCG][0-9]+", paste(rep(substr(sequence.old, base, base), as.integer(substr(sequence.old, base + 1, base + attr(bases, "match.length")[which(base == bases)] - 1))), collapse=""), sequence)
			}	
			return(sequence)	
		},
		base.repeats,
		sequence
	))
	return(sequence)
}

calcMW <- function (sequence, extra=c("5OH", "5PPP-3P", "5PPP-3OH"), adduct=c("", "Na", "K"), rxn=c("T", "C")) {
	## ONLY DESIGNED TO CALCULATE MW FOR ONE SEQUENCE AT A TIME
	if (length(sequence) > 1) {
		warning("multiple sequence input (", length(sequence), " where 1 expected)")
		sequence <- sequence[1]
	}
	sequence <- expandSequence(sequence)
	extra <- match.arg(extra)
	adduct <- match.arg(adduct)
	rxn <- match.arg(rxn)
	## MOLECULAR WEIGHT MEASUREMENTS FOR EACH RNA BASE
	massA <- 329.2098
	massG <- 345.2091
	massC <- 289.1851
	## I DO NOT!!!!! UNDERSTAND THIS ADJUSTMENT, BUT DERIVED IT FROM THE SEQUENOM PREDICTED MASSES . . . EACH ADDITIONAL "T" IN SEQUENCE IS WEIGHTED WITH SLIGHTLY DIFFERENT MOLECULAR WEIGHT DEPENDING ON IF T OR C RXN
	switch(rxn,
		## NOTE: "T" IN RNA IS ACTUALLY A URACIL
		"T" = massT <- 306.169,
		"C" = massT <- 304.1967
	)

	molecularWeight <- function (sequence) {
		## NOTE: "R" INDICATES "A" OR "G"
		if (length(grep("R", sequence)) == 1) {
			return(unique(c(
				molecularWeight(sub("R", "A", sequence)), 
				molecularWeight(sub("R", "G", sequence))))
			)	
		}
		## NOTE: "Y" INDICATES "C" OR "T"
		if (length(grep("Y", sequence)) == 1) {
			return(unique(c(
				molecularWeight(sub("Y", "T", sequence)), 
				molecularWeight(sub("Y", "C", sequence))))
			)	
		}
		## NOTE: "N" INDICATES "A", "T", "C", OR "G"
		if (length(grep("N", sequence)) == 1) {
			return(unique(c(
				molecularWeight(sub("N", "A", sequence)), 
				molecularWeight(sub("N", "T", sequence)), 
				molecularWeight(sub("N", "C", sequence)), 
				molecularWeight(sub("N", "G", sequence))))
			)	
		}
		## ADJUSTMENT FOR INITIATING 5' PHOSPHATE AND 3' HYDROXYL GROUP
		switch(extra,
			"5OH" = MW <- 19.02327,
			"5PPP-3P" = MW <- 258.96227,
			"5PPP-3OH" = MW <- 178.9823
		)
		switch(adduct,
			"Na" = MW <- MW + 21.9818,
			"K" = MW <- MW + 38.0903
		)
		#### AGAIN, I DO NOT YET UNDERSTAND WHY!!!!!  BUT I KNOW THAT THIS CORRECTION IS NECESSARY IN ORDER TO GET C RXN TO AGREE WITH SEQUENOM!!!!
		switch(rxn,
			"T" = MW <- MW,
			"C" = MW <- MW + 15.9993
		)
		## ADD MOLECULAR WEIGHTS FOR EACH BASE IN SEQUENCE
		for (i in 1:nchar(sequence)) {
			base <- substr(sequence, i, i)
			switch(base,
				"A" = MW <- MW + massA,
				"G" = MW <- MW + massG,
				"C" = MW <- MW + massC,
				"T" = MW <- MW + massT,
				return(NA)
			)
		}
		return(round(MW, 3))
	}
	return(molecularWeight(sequence))
}
