ampliconPrediction <- function(sequence, 
		lower.threshold=1500, 
		upper.threshold=7000,
		fwd.tag="AGGAAGAGAG", 
		rev.tag="AGCCTTCTCCC",
		plot=TRUE,
		table=TRUE,
		lwd=1,
		cex=1,
		multiple.conversion=FALSE) {
	sequence <- gsub("[^ACGTNXYR.()><]", "", toupper(sequence))
	fwd.tag <- gsub("[^ACTG]", "", toupper(fwd.tag))
	rev.tag <- gsub("[^ACTG]", "", toupper(rev.tag))
	lower.threshold <- as.numeric(lower.threshold)
	upper.threshold <- as.numeric(upper.threshold)
	N <- gregexpr("[CY][.()><]*[GR]", paste(fwd.tag, sequence, rev.tag, sep=""))[[1]]
	if (N[1] < 0) {
		stop("Input sequence contains no CG dinucleotides")
	}
	else {
		N <- length(N)
	}
	CpGs <- matrix(TRUE, nrow=N, ncol=6, dimnames=list(1:N, c("required", "summary", "T+", "C+", "T-", "C-")))
	if (plot) {
		plot(0, type="n", 
			xlim=c(1, 1000), 
			ylim=c(1, 5), xlab="Nucleotide position (bp)", ylab="",
			main="inSilico Assay Prediction", axes=FALSE)
		usr <- par("usr")
		scale <- function (x, min1, max1, min2=usr[1], max2=usr[2]) {
			return(min2 + (x - min1) * (max2 - min2) / (max1 - min1))
		}
	}
	## IN SILICO PREDICTION FOR BOTH FWD AND REV REACTIONS (+/- STRANDS, RESPECTIVELY)
	for (strand in c("+", "-")) {
		switch(strand,
			## STANDARD REV-T7 RXN
			"+" = sequence.i <- sequence,
			## PERFORM ASSAY ON OPPOSITE STRAND, USE FWD-T7
			"-" = {
			   sequence.i <- revComplement(sequence)
			   if (dim(CpGs)[1] > 1) {
				CpGs <- CpGs[N:1, ]
			   }
			}
		)
		required <- gregexpr("[(][ACGTNNXYR.]+[)]", sequence.i)[[1]]
		sequence.i <- gsub("[()]", "", sequence.i)
		attr(required, "match.length") <- attr(required, "match.length") - 2
		required <- required + nchar(fwd.tag)
		required <- required - 2 * (1:length(required) - 1)
		for (primer.flags in gregexpr("[<>]", sequence.i)[[1]]) {
			adjust.length <- (required <= primer.flags & required + attr(required, "match.length") - 1 >= primer.flags)
			attr(required[adjust.length], "match.length") <- attr(required[adjust.length], "match.length") - 1
			adjust.position <- (required >= primer.flags)
			required[adjust.position] <- required[adjust.position] - 1
		}
		fwd.primer <- max(0, attr(regexpr("^[^>]*(?=[>])", sequence.i, perl=TRUE), "match.length"))
		rev.primer <- max(0, attr(regexpr("(?<=[<])[^<]*$", sequence.i, perl=TRUE), "match.length"))
		sequence.i <- gsub("[<>]", "", sequence.i)
		N.seq <- nchar(paste(fwd.tag, sequence.i, rev.tag, sep=""))
		## IN SILICO PREDICTION FOR BOTH T&C REACTIONS
		for (rxn in c("T", "C")) {
			fragments <- inSilicoFragmentation(sequence.i, fwd.tag, rev.tag, type=rxn, lower.threshold, upper.threshold, fwd.primer, rev.primer, multiple.conversion)
			fragment.starts <- unlist(lapply(fragments, slot, "position"))
			fragment.ends <- fragment.starts + unlist(lapply(fragments, slot, "length")) - 1
			primer.fragments <- which((fragment.ends <= nchar(fwd.tag) + fwd.primer) | (fragment.starts >= nchar(sequence.i) + nchar(fwd.tag) - rev.primer + 1))
			for (i in primer.fragments) {
				fragments[[i]]$primer <- TRUE
			}
			for (i in required) {
				fragments.i <- which((fragment.ends >= i) & (fragment.starts <= i + attr(required, "match.length")[which(required == i)] - 1))
				for (j in fragments.i) {
					fragments[[j]]$required <- TRUE
				}
			}
			## IDENTIFY FRAGMENTS THAT HAVE CGs AND ARE NOT CONVERSION CONTROLS
			CpG.positions <- gregexpr("[CY][GR]", paste(fwd.tag, sequence.i, rev.tag, sep=""))[[1]]
			CpG.fragments <- unlist(lapply(CpG.positions, 
				function(x) {
					positions <- unlist(lapply(fragments, slot, "position"))
					lengths <- unlist(lapply(fragments, slot, "length"))
					return(which(x >= positions & x < positions + lengths))
				}
			))
			for (i in 1:N) {
				CpGs[i, "required"] <- CpGs[i, "required"] & fragments[[CpG.fragments[i]]]$required
#				if (fragments[[CpG.fragments[i]]]$CpGs > 1) {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], fragments[[CpG.fragments[i]]]$CpGs, sep="")
#				}
#				else {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], ".", sep="")
#				}
#				collisions <- fragments[[CpG.fragments[i]]]$collision.IDs
#				CpG.collisions <- unlist(lapply(collisions,
#					function(x) {
#						if (is.null(x) | identical(x, integer(0))) return(0)
#						return(length(which((unlist(lapply(fragments[x], slot, "CpGs")) > 0)
#							& !(unlist(lapply(fragments[x], slot, "conversion.control"))))))
#					}
#				))
#				if (any(CpG.collisions > 0)) {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], "D", sep="")
#				}
#				else {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], ".", sep="")
#				}
#				collisions <- fragments[[CpG.fragments[i]]]$collisions
#				if (any(collisions - CpG.collisions > 0)) {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], "O", sep="")
#				}
#				else {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], ".", sep="")
#				}
#				if (fragments[[CpG.fragments[i]]]$primer) {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], "P", sep="")
#				}
#				else {
#					CpGs[i, paste(rxn, strand, sep="")] <- paste(CpGs[i, paste(rxn, strand, sep="")], ".", sep="")
#				}				
				CpGs[i, paste(rxn, strand, sep="")] <- fragments[[CpG.fragments[i]]]$assayable
			}		
			## DISPLAY FRAGMENTATION PROFILES GRAPHICALLY
			if (plot) {
				fragment.margin <- 5.5-(which(paste(rxn, strand, sep="") == c("T+", "C+", "T-", "C-")))
				sequence.start <- scale(nchar(fwd.tag), 1, N.seq)
				sequence.end <- scale(N.seq - nchar(rev.tag) + 1, 1, N.seq)
				## PLOT SEQUENCE TAGS AND LABELS
				rect(scale(1, 1, N.seq), 
					scale(fragment.margin - 0.125, 1, 5, usr[3], usr[4]), 
					sequence.start, 
					scale(fragment.margin + 0.125, 1, 5, usr[3], usr[4]), 
					col="yellow", border="yellow"
				)
				rect(sequence.end, 
					scale(fragment.margin - 0.125, 1, 5, usr[3], usr[4]), 
					scale(N.seq, 1, N.seq), 
					scale(fragment.margin + 0.125, 1, 5, usr[3], usr[4]), 
					col="yellow", border="yellow"
				)
				if (strand == "-") {
					CpG.positions <- rev(N.seq - (CpG.positions + 1) + 1)
				}
				collision.count <- c()
				## PLOT EACH FRAGMENT, ONE AT A TIME
				for (fragment.i in fragments) {
					if (strand == "-") {
						fragment.i$position <- as.integer(N.seq - (fragment.i$position + fragment.i$length - 1) + 1)
					}
					CpGs.i <- which(CpG.positions >= fragment.i$position & CpG.positions < fragment.i$position + fragment.i$length)
					## DISPLAY FRAGMENTATION PROFILE
					segment.start <- scale(fragment.i$position, 1, N.seq)
					segment.end <- scale(fragment.i$position + fragment.i$length - 0.01, 1, N.seq)
					## HIGHLIGHT REQUIRED FRAGMENTS (USER-SPECIFIED)
					if (fragment.i$required) {
						rect(segment.start,
							scale(fragment.margin - 0.325, 1, 5, usr[3], usr[4]), 
							segment.end,
							scale(fragment.margin + 0.2, 1, 5, usr[3], usr[4]),
							col="#FFEEFF", border="#FFEEFF"
						)	
					}
					## HIGHLIGHT FRAGMENTS THAT ARE LOCATED WITHIN PRIMERS (DO NOT REPEAT HIGHLIGHTING OF TAG SEQUENCES)
					if (fragment.i$primer & fragment.i$position < N.seq - nchar(rev.tag) + 1 & fragment.i$position + fragment.i$length - 1 > nchar(fwd.tag)) {
#						rect(max(segment.start, scale(nchar(fwd.tag) + 1, 1, N.seq)), 
						rect(max(segment.start, scale(nchar(fwd.tag), 1, N.seq)), 
							scale(fragment.margin - 0.125, 1, 5, usr[3], usr[4]), 
#							min(segment.end, scale(N.seq - nchar(rev.tag), 1, N.seq)), 
							min(segment.end, scale(N.seq - nchar(rev.tag) + 1, 1, N.seq)), 
							scale(fragment.margin + 0.125, 1, 5, usr[3], usr[4]), 
							col="#FFFFCC", border="#FFFFCC"
						)
					}
					## IS THIS FRAGMENT ASSAYABLE (I.E. BETWEEN THE LOWER/UPPER MW THRESHOLDS)?
					if (fragment.i$assayable) {
						## SELECT APPROPRIATE COLORS FOR EACH FRAGMENT
						if (fragment.i$conversion.control) {
							color.i <- "green"
						}
						else if (fragment.i$CpGs > 0) {
							color.i <- "blue"	
						}
						else {
							color.i <- "black"
						}
						if (sum(fragment.i$collisions) > 0) {
							color.i <- "red" #FF6600
						}
						## PLOT FRAGMENT DELIMITERS
						segments(segment.start, 
							scale(fragment.margin - 0.125, 1, 5, usr[3], usr[4]), 
							segment.start, 
							scale(fragment.margin + 0.125, 1, 5, usr[3], usr[4]), 
							col=color.i, lwd=lwd/2
						)
						segments(segment.end, 
							scale(fragment.margin - 0.075, 1, 5, usr[3], usr[4]), 
							segment.end, 
							scale(fragment.margin + 0.075, 1, 5, usr[3], usr[4]), 
							col=color.i, lwd=lwd/2
						)
					}
					## ELSE... FRAGMENT IS NOT ASSAYABLE
					else {
						color.i <- "gray"
					}
					## FINISH PLOTTING FRAGMENTATION (HORIZONTAL SEGMENTS)
					segments(segment.start, 
						scale(fragment.margin, 1, 5, usr[3], usr[4]), 
						segment.end, 
						scale(fragment.margin, 1, 5, usr[3], usr[4]), 
						col=color.i, lwd=lwd/2
					)
					## PLOT CGs ALONG FRAGMENTATION PROFILE
					if ((fragment.i$CpGs > 0) & !fragment.i$conversion.control) {
						CpG.positions.i <- scale(CpG.positions[CpGs.i], 1, N.seq)
						points(CpG.positions.i, 
							rep(scale(fragment.margin, 1, 5, usr[3], usr[4]), length(CpGs.i)), 
							pch=19, col=color.i, cex=0.5*cex
						)
						text(CpG.positions.i, 
							rep(scale(fragment.margin + 0.225, 1, 5, usr[3], usr[4]), length(CpGs.i)), 
							labels=c(CpGs.i[1], rep("-", length(CpGs.i) - 1)), 
							cex=0.4*cex, col=color.i
						)
					}
				}
				## PLOT SEQUENCE TAGS AND LABELS
				segments(c(sequence.start, sequence.end), 
					scale(fragment.margin - 0.325, 1, 5, usr[3], usr[4]), 
					c(sequence.start, sequence.end), 
					scale(fragment.margin + 0.225, 1, 5, usr[3], usr[4]), 
					col="black", lwd=lwd
				)
				text(c(sequence.start, sequence.end), 
					scale(fragment.margin + 0.325, 1, 5, usr[3], usr[4]), 
					labels=paste(c(paste(rxn, "(", strand, ") 1", sep=""), (N.seq - nchar(fwd.tag) - nchar(rev.tag))), "bp"), 
					cex=0.75*cex
				)
				## DISPLAY WHICH CpG-CONTAINING FRAGMENTS COLLIDE WITH WHICH OTHER CpG-CONTAINING FRAGMENTS
				collision.count <- c()
				for (fragment.i in fragments) {
					if ((fragment.i$CpGs < 1) | (!fragment.i$assayable) | (length(unlist(fragment.i$CG.collision.IDs)) < 1)) next
					if (strand == "-") {
						fragment.i$position <- as.integer(N.seq - (fragment.i$position + fragment.i$length - 1) + 1)
					}
					CGs.i <- which(CpG.positions >= fragment.i$position & CpG.positions < fragment.i$position + fragment.i$length)
					points(scale(CpG.positions[CGs.i], 1, N.seq), 
						rep(scale(fragment.margin - 0.15, 1, 5, usr[3], usr[4]), length(CGs.i)), 
						pch=17, col="gray", cex=0.4*cex
					)
					CG.collisions <- unlist(fragment.i$CG.collision.IDs)
					if (strand == "-") CG.collisions <- N - CG.collisions + 1
					CG.collisions <- unique(c(CG.collisions, CGs.i))
					if (all(!collision.count %in% CG.collisions)) {
						collision.count <- c(collision.count, CGs.i[1])
						segments(scale(min(CpG.positions[CG.collisions]), 1, N.seq),
							scale(fragment.margin - 0.175 - 0.05*length(collision.count), 1, 5, usr[3], usr[4]),
							scale(max(CpG.positions[CG.collisions]), 1, N.seq),
							scale(fragment.margin - 0.175 - 0.05*length(collision.count), 1, 5, usr[3], usr[4]),
							lty="solid", col="gray", lwd=lwd/2
						)
						segments(scale(CpG.positions[CG.collisions], 1, N.seq),
							rep(scale(fragment.margin - 0.15, 1, 5, usr[3], usr[4]), length(CG.collisions)),
							scale(CpG.positions[CG.collisions], 1, N.seq),
							scale(fragment.margin - 0.175 - 0.05*length(collision.count), 1, 5, usr[3], usr[4]),
							lty="solid", col="gray", lwd=lwd/2
						)
					}
				}
			}
		}
	}
	if (plot) {
		for (i in 1:N.seq) {
			text(scale(i, 1, N.seq), 
				 scale(1.05, 1, 5, usr[3], usr[4]), 
				 labels=substr(paste(fwd.tag, revComplement(sequence.i), rev.tag, sep=""), i, i), 
				 cex=0.2*cex
				 )
		}
	}
	if (N > 1) {
		CpGs <- CpGs[N:1, ]
	}
	CpG.counts <- matrix(0, nrow=2, ncol=7, dimnames=list(c("all", "required"), c("summary", "+", "T+", "C+", "-", "T-", "C-")))
	for (i in 1:N) {
		CpGs[i, "summary"] <- any(CpGs[i, c("T+", "C+", "T-", "C-")])
		if (CpGs[i, "summary"]) CpG.counts[1, "summary"] <- CpG.counts[1, "summary"] + 1
		if (CpGs[i, "T+"]) CpG.counts[1, "T+"] <- CpG.counts[1, "T+"] + 1
		if (CpGs[i, "C+"]) CpG.counts[1, "C+"] <- CpG.counts[1, "C+"] + 1
		if (any(CpGs[i, c("T+", "C+")])) CpG.counts[1, "+"] <- CpG.counts[1, "+"] + 1
		if (CpGs[i, "T-"]) CpG.counts[1, "T-"] <- CpG.counts[1, "T-"] + 1
		if (CpGs[i, "C-"]) CpG.counts[1, "C-"] <- CpG.counts[1, "C-"] + 1
		if (any(CpGs[i, c("T-", "C-")])) CpG.counts[1, "-"] <- CpG.counts[1, "-"] + 1
		if (CpGs[i, "required"] & CpGs[i, "summary"]) CpG.counts[2, "summary"] <- CpG.counts[2, "summary"] + 1
		if (CpGs[i, "required"] & CpGs[i, "T+"]) CpG.counts[2, "T+"] <- CpG.counts[2, "T+"] + 1
		if (CpGs[i, "required"] & CpGs[i, "C+"]) CpG.counts[2, "C+"] <- CpG.counts[2, "C+"] + 1
		if (CpGs[i, "required"] & any(CpGs[i, c("T+", "C+")])) CpG.counts[2, "+"] <- CpG.counts[2, "+"] + 1
		if (CpGs[i, "required"] & CpGs[i, "T-"]) CpG.counts[2, "T-"] <- CpG.counts[2, "T-"] + 1
		if (CpGs[i, "required"] & CpGs[i, "C-"]) CpG.counts[2, "C-"] <- CpG.counts[2, "C-"] + 1
		if (CpGs[i, "required"] & any(CpGs[i, c("T-", "C-")])) CpG.counts[2, "-"] <- CpG.counts[2, "-"] + 1
	}
	return(list(summary=CpGs, counts=CpG.counts))
}