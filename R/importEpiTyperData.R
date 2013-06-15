importEpiTyperData <- function (data, MassArrayObject, verbose=TRUE) {
	samples <- list()
	assay.rows <- c(grep("^T Methyl ", as.character(data[, 1])) - 1, dim(data)[1] + 1)
	if (length(assay.rows) <= 1) {
		if (verbose) warning("File does not contain an identifiable assay")
		return(samples)
	}
	else if (length(assay.rows) > 2) {
	   	if (verbose) warning("File contains multiple assays (please export data as a single assay per file)")
	   	return(samples)
	}
	if (verbose) cat("\tReading assay (", as.character(data[assay.rows[1], 1]), ")", sep="")
	assay.data <- data[assay.rows[1]:(assay.rows[2] - 1), ]
	T.row.divider <- grep("^T Methyl ", as.character(assay.data[, 1]))
	C.row.divider <- grep("^C Methyl ", as.character(assay.data[, 1]))
	T.data <- assay.data[T.row.divider:(C.row.divider - 1), ]
	C.data <- assay.data[C.row.divider:dim(assay.data)[1], ]
	T.sample.rows <- grep(" \\[Plate.*\\]$", T.data[, 1])
	T.samples <- sub(" \\[Plate.*\\]$", "", as.character(T.data[T.sample.rows, 1]))
	T.sample.rows <- c(T.sample.rows, dim(T.data)[1] + 1) + 1
	C.sample.rows <- grep(" \\[Plate.*\\]$", as.character(C.data[, 1]))
	C.samples <- sub(" \\[Plate.*\\]$", "", as.character(C.data[C.sample.rows, 1]))
	C.sample.rows <- c(C.sample.rows, dim(C.data)[1] + 1) + 1
	if (verbose) cat(", ", length(T.samples), "T+", length(C.samples), "C rxns (EpiTyper v1.0):\n", sep="")
	for (rxn in c("T", "C")) {
		switch(rxn,
			"T" = num.samples <- length(T.samples),
			"C" = num.samples <- length(C.samples)
		)
		if (num.samples < 1) next
		if (verbose) cat("\t\t", rxn, " reaction:\n", sep="")
		for (i in 1:num.samples) {
			switch(rxn,
				"T" = {
					sample.data <- T.data[T.sample.rows[i]:(T.sample.rows[i + 1] - 2), ]
					sample.name <- T.samples[i]
					sample.fragments <- MassArrayObject@fragments.T
				},
				"C" = {
					sample.data <- C.data[C.sample.rows[i]:(C.sample.rows[i + 1] - 2), ]
					sample.name <- C.samples[i]
					sample.fragments <- MassArrayObject@fragments.C
				}
			)	
			if (verbose) cat("\t\t\tReading sample (", sample.name, ") ... ", sep="")
			sample.fragments.MW <- unlist(lapply(sample.fragments, slot, "MW"))
			peaks <- c()
			for (j in 1:dim(sample.data)[1]) {
				switch(as.character(sample.data[j, grep("Ref[._ ]type", colnames(sample.data), ignore.case=TRUE)]),
					"NAAD" = adduct <- "Na",
					"K-AD" = adduct <- "K",
					adduct <- NULL
				)
				new <- FALSE
				switch(as.character(sample.data[j, grep("Ref[._ ]type", colnames(sample.data), ignore.case=TRUE)]),
					"NAAD" = type <- "Modified",
					"K-AD" = type <- "Modified",
					"UNKN" = {
						if (length(grep("^(Doubly charged|DOUBLY CHARGED)", as.character(sample.data[j, "Description"]))) > 0) type <- "Modified"
						else type <- "Unknown"
					},
					"NEW " = {
						type <- "Unknown"
						new <- TRUE
					},
					"MAIN" = type <- "Expected",
					"METH" = type <- "Expected",
					"NOME" = type <- "Expected",
					"ACYC" = type <- "Modified",
					"MAIN/ANCH" = type <- "Expected",
					"NOME/ANCH" = type <- "Expected",
					type <- "Unknown"
				)
				if (length(grep("^Missing", as.character(sample.data[j, 1]))) > 0) missing <- TRUE
				else missing <- FALSE
				peak.count <- length(which(!is.na(findPeaks(sample.fragments.MW, sample.data[j, grep("Sample.mass", colnames(sample.data), ignore.case=TRUE)]))))
				## SPLIT PEAK DESCRIPTIONS INTO THEIR POTENTIAL AND/OR ACTUAL COMPONENTS
				description <- unlist(strsplit(as.character(sample.data[j, "Description"]), "[;,]"))
				## REMOVE EXTRANEOUS WHITE SPACES
				description <- gsub("[ ]", "", description)
				## REMOVE NUMERICAL LABELS (HOLD POSITIONAL AND/OR SCORING INFORMATION WHICH IS NOT NEEDED HERE)
				description <- sub("([(][0-9.]*[)].*|[@].*)", "", description)
				charge <- length(grep("^(Doubly charged|DOUBLY CHARGED)", description)) + 1
				## REMOVE DOUBLY CHARGED OR ADDUCT LABELS
				description <- sub(".*[Oo][Ff]", "", description)
				## EXTRACT SEQUENCE FROM 5OH-...-3P AND OTHER SIMILAR FRAGMENTS
				description <- unique(sub(".*[-]([ATCG0-9]+)[-].*", "\\1", description))
				peaks <- c(peaks, new("MassArrayPeak", ID=as.integer(j),
					MW.theoretical=sample.data[j, grep("Reference.mass", colnames(sample.data), ignore.case=TRUE)],
					MW.actual=sample.data[j, grep("Sample.mass", colnames(sample.data), ignore.case=TRUE)],
					probability=sample.data[j, "Probability"],
					SNR=sample.data[j, "SNR"],
					height=sample.data[j, "Height"],
					sample.intensity=sample.data[j, grep("Sample.intensity", colnames(sample.data), ignore.case=TRUE)],
					ref.intensity=sample.data[j, grep("Reference.intensity", colnames(sample.data), ignore.case=TRUE)],
					sequence=as.character(description),
					adduct=adduct, type=type, charge=charge, missing=missing, new=new,
					collisions=nchar(gsub("[^;]", "", as.character(sample.data[j, "Description"]))),
					components=peak.count)
				)
			}
			## EXTRACT SEQUENCE INFORMATION FROM REFERENCE PEAK FOR ADDUCTS AND DOUBLY CHARGED FRAGMENTS
			for (i in 1:length(peaks)) {
				if (length(grep("^[0-9]+[.]?[0-9]*$", peaks[[i]]$sequence)) < 1) next
				peaks[[i]]$sequence <- as.character(peaks[[findPeaks(as.numeric(peaks[[i]]$sequence), unlist(lapply(peaks, slot, "MW.theoretical")))]]$sequence)
			}
			## INTRODUCE ARTIFICIAL CONVERSION CONTROLS INTO PEAK SPECTRA (IF NOT ALREADY PRESENT)
			peaks.expected <- lapply(sample.fragments, slot, "MW")
			peaks.expected <- peaks.expected[unlist(lapply(sample.fragments, slot, "conversion.control"))]
			peaks.expected <- unique(unlist(peaks.expected))
			peaks.observed <- unlist(unique(lapply(peaks, slot, "MW.actual")))
			peaks.found <- findPeaks(peaks.expected, peaks.observed, resolution=1)
			if (is.null(peaks.found)) {
				control.peaks.missing <- NULL
			}
			else {
				control.peaks.missing <- which(is.na(peaks.found))
			}
			if (length(control.peaks.missing) > 0) {
				for (k in 1:length(control.peaks.missing)) {
					peaks <- c(peaks, new("MassArrayPeak", ID=j+k,
						MW.theoretical=peaks.expected[control.peaks.missing[k]],
						MW.actual=peaks.expected[control.peaks.missing[k]],
						probability=NA,
						SNR=0,
						height=NA,
						sample.intensity=NA,
						ref.intensity=1,
						adduct=NULL, type="Control", missing=TRUE, new=FALSE,
						collisions=0, components=1))
				}
			}
			samples <- c(samples, new("MassArraySpectrum", sample=sample.name,
					rxn=rxn, strand=MassArrayObject@strand, peaks=peaks)
				)
			if (verbose) cat("FINISHED\n")
		}
	}
	return(samples)
}



importEpiTyperData.new <- function (data, MassArrayObject, verbose=TRUE) {
	samples <- list()
	if (length(unique(as.character(data[, "Amplicon"]))) < 1) {
		if (verbose) warning("File does not contain an identifiable assay")
		return(samples)
	}
	else if (length(unique(as.character(data[, "Amplicon"]))) > 1) {
	    if (verbose) warning("File contains multiple assays (please export data as a single assay per file)")
	    return(samples)
	}
	if (verbose) cat("\tReading assay (", unique(as.character(data[, "Amplicon"])), ")", sep="")
	T.data <- data[grep("^T Methyl ", as.character(data[, "Reaction"])), ]
	C.data <- data[grep("^C Methyl ", as.character(data[, "Reaction"])), ]
	T.samples <- unique(as.character(T.data[, "Sample"]))
	C.samples <- unique(as.character(C.data[, "Sample"]))
	if (verbose) cat(", ", length(T.samples), "T+", length(C.samples), "C rxns (EpiTyper v1.0.5):\n", sep="")
	for (rxn in c("T", "C")) {
		switch(rxn,
			"T" = num.samples <- length(T.samples),
			"C" = num.samples <- length(C.samples)
		)
		if (num.samples < 1) next
		if (verbose) cat("\t\t", rxn, " reaction:\n", sep="")
		for (i in 1:num.samples) {
			switch(rxn,
				"T" = {
					sample.data <- T.data[grep(T.samples[i], as.character(T.data[, "Sample"]), fixed=TRUE), ]
					sample.name <- sub(" \\[Plate.*\\]$", "", T.samples[i])
					sample.fragments <- MassArrayObject@fragments.T
				},
				"C" = {
					sample.data <- C.data[grep(C.samples[i], as.character(C.data[, "Sample"]), fixed=TRUE), ]
					sample.name <- sub(" \\[Plate.*\\]$", "", C.samples[i])
					sample.fragments <- MassArrayObject@fragments.C
				}
			)
			if (verbose) cat("\t\t\tReading sample (", sample.name, ") ... ", sep="")
			sample.fragments.MW <- unlist(lapply(sample.fragments, slot, "MW"))
			peaks <- c()
			for (j in 1:dim(sample.data)[1]) {
				switch(as.character(sample.data[j, "Ref.type"]),
					"NAAD" = adduct <- "Na",
					"K-AD" = adduct <- "K",
					adduct <- NULL
				)
				new <- FALSE
				switch(as.character(sample.data[j, "Ref.type"]),
					"NAAD" = type <- "Modified",
					"K-AD" = type <- "Modified",
					"UNKN" = {
						if (length(grep("^(Doubly charged|DOUBLY CHARGED)", as.character(sample.data[j, "Description"]))) > 0) type <- "Modified"
						else type <- "Unknown"
					},
					"NEW " = {
						type <- "Unknown"
						new <- TRUE
					},
					"MAIN" = type <- "Expected",
					"METH" = type <- "Expected",
					"NOME" = type <- "Expected",
					"ACYC" = type <- "Modified",
					"MAIN/ANCH" = type <- "Expected",
					"NOME/ANCH" = type <- "Expected",
					type <- "Unknown"
				)
				if (length(grep("^Missing", as.character(sample.data[j, "Peak"]))) > 0) missing <- TRUE
				else missing <- FALSE
				peak.count <- length(which(!is.na(findPeaks(sample.fragments.MW, sample.data[j, grep("Sample.mass", colnames(sample.data), ignore.case=TRUE)]))))
				## SPLIT PEAK DESCRIPTIONS INTO THEIR POTENTIAL AND/OR ACTUAL COMPONENTS
				description <- unlist(strsplit(as.character(sample.data[j, "Description"]), "[;,]"))
				## REMOVE EXTRANEOUS WHITE SPACES
				description <- gsub("[ ]", "", description)
				## REMOVE NUMERICAL LABELS (HOLD POSITIONAL AND/OR SCORING INFORMATION WHICH IS NOT NEEDED HERE)
				description <- sub("([(][0-9.]*[)].*|[@].*)", "", description)
				## IDENTIFY AND THEN REMOVE DOUBLY CHARGED OR ADDUCT LABELS
				charge <- length(grep("^(Doubly charged|DOUBLY CHARGED)", description)) + 1
				description <- sub(".*[Oo][Ff]", "", description)
				## EXTRACT SEQUENCE FROM 5OH-...-3P AND OTHER SIMILAR FRAGMENTS
				description <- unique(sub(".*[-]([ATCG0-9]+)[-].*", "\\1", description))
				peaks <- c(peaks, new("MassArrayPeak", ID=as.integer(j),
					MW.theoretical=sample.data[j, grep("Reference.mass", colnames(sample.data), ignore.case=TRUE)],
					MW.actual=sample.data[j, grep("Sample.mass", colnames(sample.data), ignore.case=TRUE)],
					probability=sample.data[j, "Probability"],
					SNR=sample.data[j, "SNR"],
					height=sample.data[j, "Height"],
					sample.intensity=sample.data[j, grep("Sample.intensity", colnames(sample.data), ignore.case=TRUE)],
					ref.intensity=sample.data[j, grep("Reference.intensity", colnames(sample.data), ignore.case=TRUE)],
					sequence=as.character(description),
					adduct=adduct, type=type, charge=charge, missing=missing, new=new,
					collisions=nchar(gsub("[^;]", "", as.character(sample.data[j, "Description"]))),
					components=peak.count)
				)
			}
			## EXTRACT SEQUENCE INFORMATION FROM REFERENCE PEAK FOR ADDUCTS AND DOUBLY CHARGED FRAGMENTS
			for (i in 1:length(peaks)) {
				if (length(grep("^[0-9]+[.]?[0-9]*$", peaks[[i]]$sequence)) < 1) next
				peaks[[i]]$sequence <- as.character(peaks[[findPeaks(as.numeric(peaks[[i]]$sequence), unlist(lapply(peaks, slot, "MW.theoretical")))]]$sequence)
			}
			## INTRODUCE ARTIFICIAL CONVERSION CONTROLS INTO PEAK SPECTRA (IF NOT ALREADY PRESENT)
			peaks.expected <- lapply(sample.fragments, slot, "MW")
			peaks.expected <- peaks.expected[unlist(lapply(sample.fragments, slot, "conversion.control"))]
			peaks.expected <- unique(unlist(peaks.expected))
			peaks.observed <- unlist(unique(lapply(peaks, slot, "MW.actual")))
			peaks.found <- findPeaks(peaks.expected, peaks.observed, resolution=1)
			if (is.null(peaks.found)) {
				control.peaks.missing <- NULL
			}
			else {
				control.peaks.missing <- which(is.na(peaks.found))
			}
			if (length(control.peaks.missing) > 0) {
				for (k in 1:length(control.peaks.missing)) {
					peaks <- c(peaks, new("MassArrayPeak", ID=j+k,
						MW.theoretical=peaks.expected[control.peaks.missing[k]],
						MW.actual=peaks.expected[control.peaks.missing[k]],
						probability=NA,
						SNR=0,
						height=NA,
						sample.intensity=NA,
						ref.intensity=1,
						adduct=NULL, type="Control", missing=TRUE, new=FALSE,
						collisions=0, components=1))
				}
			}
			samples <- c(samples, new("MassArraySpectrum", sample=sample.name,
					rxn=rxn, strand=MassArrayObject@strand, peaks=peaks)
				)
			if (verbose) cat("FINISHED\n")
		}
	}
	return(samples)
}