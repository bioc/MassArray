createWiggle <- function(x, file="", append=FALSE, colors=NULL, na.rm=FALSE, sep=" ") {
	N <- length(x$samples)
	if (N < 1) {
		stop("'MassArrayData' object must contain at least one sample")
	}
	position <- position(x)
	chr <- sub("^(chr.*)[:].*$", "\\1", position)
	start <- as.numeric(sub("^.*:([0-9]*)[-].*$", "\\1", position))
	end <- as.numeric(sub("^.*[-]([0-9]*)$", "\\1", position))
	if ((chr == "") | is.na(start) | is.na(end)) {
		stop("'MassArrayData' object does not contain positional information (i.e. 'chr1:27718536-27718918')")
	}
	if (any(is.na(x$CpG.data)) & !na.rm) {
		stop("'MassArrayData' object contains missing methylation values")
	}
	CpG.positions <- gregexpr("[CY][GR]", x$sequence)[[1]]
	sp <- getOption("scipen")
	options(scipen=10)
	if (!append) {
		cat("browser full altGraph\n", file=file)
	}
	hex.map <- matrix(data=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),dimnames=list(c(0,1,2,3,4,5,6,7,8,9,"A","B","C","D","E","F"),10))
	hex.to.10 <- function (hex) {
		digits <- unlist(strsplit(hex,split=""))
		base.10 <- 0
		for (i in 1:length(digits)) {
			base.10 <- base.10 + hex.map[digits[i],]*(16**(length(digits)-i))
		}
		return(base.10)
	}
	if (is.null(colors)) {
		colors <- rainbow(n=N, start=0, end=2/3)
	}
	if (length(colors) != N) {
		colors <- rainbow(n=N, start=0, end=2/3)
	}
	for (i in 1:length(colors)) {
		colors[i] <- paste(hex.to.10(substr(colors[i], 2, 3)),
							hex.to.10(substr(colors[i], 4, 5)),
							hex.to.10(substr(colors[i], 6, 7)), sep=",")
	}
	for (i in 1:N) {
		cat("track type=wiggle_0 name=\"", samples(x)[i], "\" description=\"MassArrayData\" visibility=full autoScale=off color=", colors[i], " altColor=", colors[i], " priority=", i, " yLineOnOff=on yLineMark=0\n", file=file, append=TRUE, sep="")
		if (any(is.na(x$CpG.data[i, ]))) {
			not.na <- which(!is.na(x$CpG.data[i, ]))
			apply(cbind(chr[not.na], start[not.na] + CpG.positions[not.na] - 1, start[not.na] + CpG.positions[not.na], x$CpG.data[i, not.na]), 1, cat, "\n", sep=sep, file=file, append=TRUE)
			rm(not.na)
		}
		else {
			apply(cbind(chr, start + CpG.positions - 1, start + CpG.positions, x$CpG.data[i, ]), 1, cat, "\n", sep=sep, file=file, append=TRUE)
		}
	}
	
	
	
	options(scipen=sp)	
}

