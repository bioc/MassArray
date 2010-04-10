evaluateSNPs <- function(x, verbose=TRUE, plot=TRUE) {

	## A FEW IMPORTANT NOTES:
	## 1) calculations are SUPER-SLOW RIGHT NOW . . . limiting step is the N-scanning SNP detection . . .
	##		current implementation loops through N-scan as many times as there are novel sequences!
	##		better implementation will be to do N-scanning once and through each N-scan test which sequences from set match!!!
	##		will likely see dramatic performance improvement! . . . just make sure output for new version same as old!
	## 2) does not factor in conversion controls with SNP new peak analysis . . . determine if fragment could be conversion
	##		control AND if SNP identity would be like the converted sequence => if so, then IGNORE as a SNP => will consider
	##		when analyze for conversion controls

	if ((class(x) != "MassArrayData") | !validObject(x)) {
		stop("not a valid object of 'MassArrayData'")
	}
	N <- length(x$samples)
	if (N < 1) stop("no samples to analyze")
	sequence.tagged <- revComplement(paste(x$fwd.tag, bisConvert(x$sequence), x$rev.tag, sep=""))
	N.seq <- nchar(sequence.tagged)
	rxns <- unique(unlist(lapply(x$samples, slot, "rxn")))
	SNP.list <- c()
	if (plot) {
		plot(0, type="n", 
			xlim=c(1, 1000), 
			ylim=c(1, 5), xlab="Nucleotide position (bp)", ylab="",
			main="inSilico SNP Prediction", axes=FALSE)
		usr <- par("usr")
		scale <- function (x, min1, max1, min2=usr[1], max2=usr[2]) {
			return(min2 + (x - min1) * (max2 - min2) / (max1 - min1))
		}
	}
	for (rxn in rxns) {
		which.rxn <- grep(rxn, unlist(lapply(x$samples, slot, "rxn")))
		if (length(which.rxn) < 1) next
		switch(rxn,
			"T" = fragments <- x$fragments.T,
			"C" = fragments <- x$fragments.C
		)
		fragment.starts <- unlist(lapply(fragments, slot, "position"))
		fragment.ends <- fragment.starts + unlist(lapply(fragments, slot, "length")) - 1
		fragments.total <- length(fragments)
		## IDENTIFY FRAGMENTS THAT HAVE CGs AND ARE NOT CONVERSION CONTROLS
		CpG.positions <- gregexpr("[CY][GR]", paste(x$fwd.tag, x$sequence, x$rev.tag, sep=""))[[1]]
		CpG.fragments <- unlist(lapply(CpG.positions, 
			function(x) {
				positions <- unlist(lapply(fragments, slot, "position"))
				lengths <- unlist(lapply(fragments, slot, "length"))
				return(which(x >= positions & x < positions + lengths))
			}
		))
		if (plot) {
			## DISPLAY FRAGMENTATION PROFILES GRAPHICALLY
			fragment.margin <- 4.5-2*(which(rxn == rxns) - 1)
			sequence.start <- scale(nchar(x$fwd.tag), 1, N.seq)
			sequence.end <- scale(N.seq - nchar(x$rev.tag) + 1, 1, N.seq)
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
			segments(c(sequence.start, sequence.end), 
				scale(fragment.margin - 0.325, 1, 5, usr[3], usr[4]), 
				c(sequence.start, sequence.end), 
				scale(fragment.margin + 0.225, 1, 5, usr[3], usr[4]), 
				col="black", lwd=1
			)
			text(c(sequence.start, sequence.end), 
				scale(fragment.margin + 0.325, 1, 5, usr[3], usr[4]), 
				labels=paste(c(paste(rxn, ": 1", sep=""), (N.seq - nchar(x$fwd.tag) - nchar(x$rev.tag))), "bp"), 
				cex=0.75
			)
#			if (strand == "-") {
#				CpG.positions <- rev(N.seq - (CpG.positions + 1) + 1)
#			}
			## PLOT EACH FRAGMENT, ONE AT A TIME
			for (fragment.i in fragments) {
#				if (strand == "-") {
#					fragment.i$position <- as.integer(N.seq - (fragment.i$position + fragment.i$length - 1) + 1)
#				}
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
				if (fragment.i$primer & fragment.i$position < N.seq - nchar(x$rev.tag) + 1 & fragment.i$position + fragment.i$length - 1 > nchar(x$fwd.tag)) {
					rect(max(segment.start, scale(nchar(x$fwd.tag) + 1, 1, N.seq)), 
						scale(fragment.margin - 0.125, 1, 5, usr[3], usr[4]), 
						min(segment.end, scale(N.seq - nchar(x$rev.tag), 1, N.seq)), 
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
						col=color.i, lwd=0.5
					)
					segments(segment.end, 
						scale(fragment.margin - 0.075, 1, 5, usr[3], usr[4]), 
						segment.end, 
						scale(fragment.margin + 0.075, 1, 5, usr[3], usr[4]), 
						col=color.i, lwd=0.5
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
					col=color.i, lwd=0.5
				)
				## PLOT CGs ALONG FRAGMENTATION PROFILE
				if ((fragment.i$CpGs > 0) & !fragment.i$conversion.control) {
					CpG.positions.i <- scale(CpG.positions[CpGs.i], 1, N.seq)
					points(CpG.positions.i, 
						rep(scale(fragment.margin, 1, 5, usr[3], usr[4]), length(CpGs.i)), 
						pch=19, col=color.i, cex=0.5
					)
					text(CpG.positions.i, 
						rep(scale(fragment.margin + 0.225, 1, 5, usr[3], usr[4]), length(CpGs.i)), 
						labels=c(CpGs.i[1], rep(".", length(CpGs.i) - 1)), 
						cex=0.4, col=color.i
					)
				}
			}
		}
		
		new.peaks <- unlist(lapply(x$samples[which.rxn], slot, "peaks"), recursive=FALSE)[which(unlist(lapply(unlist(lapply(x$samples[which.rxn], slot, "peaks"), recursive=FALSE), slot, "new")))]
		if (verbose) cat("Examining ", rxn, " reaction for SNPs: ", length(which.rxn), " samples, ", length(new.peaks), " 'NEW' peaks ... ", sep="")
		new.peaks.sequence <- unique(unlist(lapply(new.peaks, slot, "sequence")))
		## IGNORE SEQUENCES THAT ARE NOT WITHIN THE ASSAYABLE RANGE
		new.peaks.sequence <- new.peaks.sequence[unlist(lapply(new.peaks.sequence, 
			function(sequence) {
				return(lapply(lapply(calcMW(sequence), isAssayable), any))
			}
		))]
		## SEARCH FOR POTENTIAL SNPs IN SEQUENCE, ONE AT A TIME (VERY SLOW STEP) WITH AS FEW AS POSSIBLE
		new.peak.SNPs <- lapply(new.peaks.sequence, identifySNPs, sequence.tagged, rxn)
		print("done")
		if (verbose) cat("FINISHED (", length(new.peak.SNPs), " potential SNPs identified)\n", sep="")
		## MODIFY SNP POSITIONS TO CORRESPOND TO PROPER FWD-STRANDED SEQUENCE
		new.peak.SNPs <- lapply(new.peak.SNPs,
			function(SNPs) {
				if ((length(SNPs)) < 1) return()
				for (i in 1:length(SNPs)) {
					SNPs[[i]]$position <- nchar(sequence.tagged) - SNPs[[i]]$position + 1
					SNPs[[i]]$fragment <- which(fragment.starts <= SNPs[[i]]$position & fragment.ends >= SNPs[[i]]$position)
				}
				return(SNPs)
			}
		)
		best.SNPs <- c()
		## ALL POTENTIAL SNPs HAVE BEEN IDENTIFIED, NOW ITS TIME TO COMPARE ON A SAMPLE-BY-SAMPLE BASIS
		for (i in 1:length(which.rxn)) {
			sample <- x$samples[[which.rxn[i]]]
			## CALCULATE AVERAGE PER-FRAGMENT SNR
			SNR.total <- sum(unlist(lapply(sample$peaks, slot, "SNR")), na.rm=TRUE)
			SNR.avg <- SNR.total / fragments.total
			## IDENTIFY MISSING EXPECTED PEAKS
			missing.peaks <- sample$peaks[which(unlist(lapply(sample$peaks, slot, "missing")) & (unlist(lapply(sample$peaks, slot, "type")) == "Expected"))]
			## MAP EACH MISSING PEAK TO ITS CORRESPONDING EXPECTED FRAGMENT
			missing.fragments <- unlist(lapply(unlist(lapply(unlist(lapply(missing.peaks, slot, "MW.theoretical")), findFragments, fragments)), slot, "ID"))
			## GET NEW PEAKS (AND THEIR PREDICTED SEQUENCES) ON A PER-SAMPLE BASIS
			new.sample.peaks <- sample$peaks[which(unlist(lapply(sample$peaks, slot, "new")))]
			new.sample.sequences <- unlist(lapply(unique(unlist(lapply(new.sample.peaks, slot, "sequence"))), expandSequence))
			## FIND SUBSET OF SNPs THAT MATCH NEW PEAK SEQUENCES FOR THIS SAMPLE
			sample.SNPs <- lapply(new.peak.SNPs,
				function(SNPs) {
					if (length(SNPs) < 1) return()
					if (length(SNPs[[1]]$sequence) < 1) return()
					if (length(new.sample.peaks) < 1) return()
					matched.peak <- which(unlist(lapply(new.sample.peaks, 
						function(peak) {
							return(SNPs[[1]]$sequence %in% unlist(lapply(peak$sequence, expandSequence)))
						}
					)))
#					matched.sequence <- match(SNPs[[1]]$sequence, new.sample.sequences)
					if (length(matched.peak) < 1) return()
					SNR <- max(unlist(lapply(new.sample.peaks[matched.peak], slot, "SNR")), na.rm=TRUE)
					if (is.na(SNR)) return()
					if (SNR <= 0) return()
					for (i in 1:length(SNPs)) {
						## STORE SNR DATA (NOTE THAT THIS ASSUMES THAT A NEW PEAK ARISES FROM ONE SIGNAL, NOT MULTIPLE)
						SNPs[[i]]$SNR <- SNR
						SNPs[[i]]$SNP.present <- min(1, max(0, SNR / SNR.avg))
						if (SNPs[[i]]$fragment %in% missing.fragments) {
							SNPs[[i]]$WT.absent <- 1
						}
						else {
							matched.peaks <- findPeaks(fragments[[SNPs[[i]]$fragment]]$MW, unlist(lapply(sample$peaks, slot, "MW.theoretical")))
#							matched.peaks <- unique(c(matched.peaks, findPeaks(fragments[[SNPs[[i]]$fragment]]$MW, unlist(lapply(sample$peaks, slot, "MW.actual")))))
							matched.peaks <- sample$peaks[matched.peaks[!is.na(matched.peaks)]]
							if (length(matched.peaks) < 1) {
								SNPs[[i]]$WT.absent <- 0
								next
							}
							count.estimate <- sum(unlist(lapply(matched.peaks,
								function(peak) {
									return(min(1, peak$components) - peak$SNR / SNR.avg)
								}
							)), na.rm=TRUE)
							if (is.na(count.estimate)) {
								SNPs[[i]]$WT.absent <- 0
								next
							}
							if (count.estimate > 0) {
								SNPs[[i]]$WT.absent <- 1
							}
							else {
								components <- sum(pmax(unlist(lapply(matched.peaks, slot, "components")), 1), na.rm=TRUE)
								SNR.sum <- sum(unlist(lapply(matched.peaks, slot, "SNR")), na.rm=TRUE)
								SNPs[[i]]$WT.absent <- min(1, max(0, components - SNR.sum / SNR.avg))
							}
						}
					}
					## ALSO MAKE SURE TO RETURN SNR DATA FOR THIS SAMPLE
					return(SNPs)
				}
			)
			## GET LIST OF FRAGMENTS CORRESPONDING TO SNPs
			sample.SNP.fragments <- unlist(lapply(unlist(sample.SNPs, recursive=FALSE), 
				function(SNP) {
					return(SNP$fragment)
				}
			))
			sample.SNP.fragments <- duplicated(sample.SNP.fragments)
			## CALCULATE APPROPRIATE SAMPLE POSITION INFORMATION FOR PLOTTING
			sample.position <- fragment.margin - 0.5 - (3/length(rxns))*(i - 1)/length(which.rxn)
			## PLOT MISSING PEAKS
			if (plot) {
				text(scale(1, 1, N.seq),
					scale(sample.position, 1, 5, usr[3], usr[4]),
					labels=sample$sample, cex=0.25, col="black", pos=4, offset=0.5, srt=90
				)
			}
			for (SNP in unlist(sample.SNPs, recursive=FALSE)) {
				if (!fragments[[SNP$fragment]]$assayable) {
					next
				}
				## MAP SNP POSITION TO THE CENTER OF ITS CONTAINING FRAGMENT
				SNP.position <- fragments[[SNP$fragment]]$position + fragments[[SNP$fragment]]$length/2
				if ((SNP$SNP.present == 1) & (SNP$WT.absent == 1)) {
					color.SNP <- "black"	
					cex.SNP <- 0.75
				}
				else if ((SNP$SNP.present == 1) | (SNP$WT.absent == 1)) {
					color.SNP <- "gray"	
					cex.SNP <- 0.25
				}
				else {
					next
				}
				SNP.map <- match(SNP$sequence, unlist(lapply(best.SNPs, 
						function(SNP) {
							return(SNP$sequence)
						}
					), recursive=FALSE)
				) 
				SNP.map <- SNP.map[which(!is.na(SNP.map))]
				if (length(SNP.map) > 0) {
					best.SNPs[[SNP.map]]$SNP <- c(best.SNPs[[SNP.map]]$SNP, paste(SNP$position, ":", SNP$base, sep=""))
					best.SNPs[[SNP.map]]$SNR <- c(best.SNPs[[SNP.map]]$SNR, SNP$SNR)
					best.SNPs[[SNP.map]]$fragment <- c(best.SNPs[[SNP.map]]$fragment, SNP$fragment)
					best.SNPs[[SNP.map]]$SNP.quality <- c(best.SNPs[[SNP.map]]$SNP.quality, (SNP$SNP.present + SNP$WT.absent))
					best.SNPs[[SNP.map]]$samples <- c(best.SNPs[[SNP.map]]$samples, sample$sample)
					best.SNPs[[SNP.map]]$count <- best.SNPs[[SNP.map]]$count + 1
				}
				else {
					best.SNPs <- c(best.SNPs, list(list(sequence=SNP$sequence, fragment=SNP$fragment, SNR=SNP$SNR, SNP=paste(SNP$position, ":", SNP$base, sep=""), SNP.quality=(SNP$SNP.present + SNP$WT.absent), samples=sample$sample, count=1)))
				}
				if (SNP$base %in% c("A", "T", "C", "G")) {
					pch.SNP <- 2
				}
				else if (SNP$base == "") {
					pch.SNP <- "O"
					cex.SNP <- cex.SNP + 0.25
				}					
				else {
					pch.SNP <- 20
				}	
				if (SNP$fragment %in% missing.fragments) {
					if (plot) {
						text(
							scale(SNP.position, 1, N.seq),
							scale(sample.position, 1, 5, usr[3], usr[4]),
							labels="X", col="red", cex=1
						)
					}
				}
				if (plot) {
					points(scale(SNP.position, 1, N.seq),
						scale(sample.position, 1, 5, usr[3], usr[4]),
						pch=pch.SNP, cex=cex.SNP, col=color.SNP
					)	
				}
				if (SNP$fragment %in% sample.SNP.fragments) {
					if (plot) {
						points(scale(SNP.position, 1, N.seq),
							scale(sample.position, 1, 5, usr[3], usr[4]),
							pch=20, cex=0.5, col="yellow"
						)	
					}
				}
			}
			SNP.list <- c(SNP.list, list(list(sample=sample$sample, SNPs=sample.SNPs)))
		}
	}
	return(best.SNPs)
#	return(SNP.list)
}