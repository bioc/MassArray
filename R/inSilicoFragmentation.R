inSilicoFragmentation <- function (sequence, fwd.tag="", rev.tag="", type=c("T", "C"), lower.threshold=1500, upper.threshold=7000, fwd.primer=0, rev.primer=0, multiple.conversion=FALSE) {
	type <- match.arg(type)
	bisulfite.sequence <- bisConvert(sequence)
	bisulfite.sequence <- paste(fwd.tag, bisulfite.sequence, rev.tag, sep="")
	## PERFORM IN SILICO DIGEST TO GET LIST OF PREDICTED FRAGMENTS
	revcomp.bis.sequence <- revComplement(bisulfite.sequence)
	CG.positions <- unlist(gregexpr("(CG|YG|CR)", bisulfite.sequence))
	fragments <- rnaDigest(revcomp.bis.sequence, type)
	## RE-ASSIGN POSITIONS TO MAP BACK TO REGULAR (I.E. NOT REV-COMPLEMENTED) SEQUENCE
	fragments <- lapply(fragments, function (x) {
			x$position <- as.integer(nchar(revcomp.bis.sequence) - (x$position + x$length - 1) + 1)
			return(x)
		}
	)
	## EXTRACT POSITIONAL INFORMATION
	positions <- unlist(lapply(fragments, slot, "position"))
	lengths <- unlist(lapply(fragments, slot, "length"))
	ends <- positions + lengths - 1
	## STORE PRIMER INFORMATION IN FRAGMENTATION PATTERN
	primer.fragments <- which((ends <= nchar(fwd.tag) + fwd.primer) | (positions >= nchar(sequence) + nchar(fwd.tag) - rev.primer + 1))
	for (i in primer.fragments) {
		fragments[[i]]$primer <- TRUE
	}
	## 5' MODIFICATION OF FIRST FRAGMENT => UPDATE MOLECULAR WEIGHT
	fragments[[1]]$extra <- "5PPP-3P"
	fragments[[1]]$MW <- calcMW(fragments[[1]]$sequence, "5PPP-3P", rxn=type)
	## GET EXPECTED PEAKS
	fragment.peaks <- lapply(fragments, slot, "MW")
	fragment.IDs <- mapply(rep, 1:length(fragments), lapply(fragment.peaks, length))
	peaks.unique <- unique(unlist(fragment.peaks))
	peaks.unique <- peaks.unique[order(peaks.unique)]
	## WHICH PEAKS ARE ASSAYABLE?
	peaks.assayable <- which(isAssayable(peaks.unique, lower.threshold, upper.threshold))
	## MAP PEAKS TO FRAGMENTS (ONE PEAK MAY MAP TO MULTIPLE)
	peak.fragment.map <- mapply(match, rep(list(peaks.unique), length(fragment.peaks)), fragment.peaks)
	peak.fragment.map <- apply(!apply(peak.fragment.map, 2, is.na), 1, which)
	## IDENTIFY ANY PEAK COLLISIONS (UNIQUE PEAKS, BUT MW TOO CLOSE TO DISTINGUISH SEPARATELY)
	peak.collisions <- findCollisions(peaks.unique)
	for (i in which(lapply(peak.collisions, length) > 1)) {
		for (j in peak.collisions[[i]]) {
			peak.fragment.map[[i]] <- unique(c(peak.fragment.map[[i]], peak.fragment.map[[j]]))
		}	
	}
	## DETERMINE NUMBER OF COMPONENT FRAGMENTS FOR EACH PEAK
	peak.counts <- unlist(lapply(peak.fragment.map, length))
	## STORE PER-FRAGMENT COLLISION DATA (INCLUDES MAP TO OTHER FRAGMENTS FOR EACH PEAK COLLISION)
	CG.fragments <- unlist(lapply(CG.positions, 
		function(x) {
			positions <- unlist(lapply(fragments, slot, "position"))
			lengths <- unlist(lapply(fragments, slot, "length"))
			return(which(x >= positions & x < positions + lengths))
		}
	))
	fragments <- lapply(fragments, function (x) {
			x$collisions <- as.integer(peak.counts[match(x$MW, peaks.unique)] - 1)
			x$collision.IDs <- lapply(peak.fragment.map[match(x$MW, peaks.unique)], setdiff, x$ID)
			## IS THIS A CpG-CONTAINING FRAGMENT?
			if (x$CpGs > 0) {
				## WHICH OF THE FRAGMENT COLLISIONS ALSO CONTAIN CpG?
				CpG.collisions <- lapply(x$collision.IDs,
					function(x) {
						if (is.null(x) | identical(x, integer(0))) return()
						return(x[which((unlist(lapply(fragments[x], slot, "CpGs")) > 0)
							& !(unlist(lapply(fragments[x], slot, "conversion.control"))))])
					}
				)
				x$CG.collision.IDs <- lapply(CpG.collisions,
					function(x) {
						return(which(CG.fragments %in% x))
					}
				)
			}
			return(x)
		}
	)
	## ARE THE FRAGMENTS ASSAYABLE?
	fragments.assayable <- unlist(lapply(relist(isAssayable(unlist(fragment.peaks), lower.threshold, upper.threshold), fragment.peaks), all))
	for (i in which(fragments.assayable)) {
		fragments[[i]]$assayable <- TRUE	
	}
	fragments <- convControl(paste(fwd.tag, sequence, rev.tag, sep=""), fragments, multiple=multiple.conversion)

	return(fragments)
}
