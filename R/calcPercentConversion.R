calcPercentAdduct <- function(peaks) {
	sequences <- lapply(peaks, slot, "sequence")
	## IDENTIFY ADDUCTS AND THEIR MATCHING REFERENCE PEAKS
	adduct.peaks <- which((unlist(lapply(peaks, slot, "adduct")) != "") & (lapply(sequences, length) == 1))
	adduct.peaks <- unique(unlist(lapply(peaks[adduct.peaks], slot, "sequence")))
	if (length(adduct.peaks) < 1) { return(numeric(0)) }
	
	## IDENTIFY MATCHING REFERENCE PEAKS
	adduct.ratios <- unlist(lapply(adduct.peaks,
		function (adduct) {
			## WHICH PEAKS MATCH THE SEQUENCE FOR THIS ADDUCT PEAK?
			matched.by.sequence <- which(unlist(lapply(lapply(sequences, "%in%", adduct), any)))
			## IDENTIFY ADDUCTS
			adducts <- matched.by.sequence[which(unlist(lapply(peaks[matched.by.sequence], slot, "adduct")) != "")]
			adducts <- sum(unlist(lapply(peaks[adducts], slot, "SNR")))
			## FILTER OUT ADDUCT, DOUBLY CHARGED, ABORTIVE CYCLING, AND OTHER MODIFIED PEAKS
			non.adducts <- matched.by.sequence[which(unlist(lapply(peaks[matched.by.sequence], slot, "type")) != "Modified")]
			non.adducts <- sum(unlist(lapply(peaks[non.adducts], slot, "SNR")))
			## CALCULATE RATIO OF ADDUCT PEAK HEIGHTS TO THEIR REFERENCES
			return(adducts / (adducts + non.adducts))
		}
	))

	return(median(adduct.ratios))
}


calcPercentConversion <- function(fragments, peaks) {
	controls <- which(unlist(lapply(fragments, slot, "conversion.control")))
	controls.meth <- c()
	for (i in controls) {
		frags.i <- fragments[[i]]$"MW"
		peaks.i <- findPeaks(frags.i, unlist(lapply(peaks, slot, "MW.actual")))
		peaks.i.theory <- findPeaks(frags.i, unlist(lapply(peaks, slot, "MW.theoretical")))
		peaks.i[which(is.na(peaks.i))] <- peaks.i.theory[which(is.na(peaks.i))]
		SNRs.i <- unlist(lapply(peaks[peaks.i], slot, "SNR"))
## EXTRACT SIGNAL-TO-NOISE RATIO (SNR) INFORMATION FOR EACH PEAK
		SNRs.i <- peaks.i
		if (any(!is.na(peaks.i))) {
			SNRs.i[which(!is.na(peaks.i))] <- unlist(lapply(peaks[peaks.i[which(!is.na(peaks.i))]], slot, "SNR"))
		}
## MAP ALL COLLIDING FRAGMENTS TO LIST OF MOLECULAR WEIGHTS
		collisions.i <- lapply(frags.i, 
							   function(x) {
							   temp <- findFragments(x, fragments)
							   IDs <- unlist(lapply(temp, slot, "ID"))
							   return(unique(IDs))
							   }
							   )
## IDENTIFY WHICH FRAGMENTS ARE NON-CONTRIBUTORY TO METHYLATION STATE (E.G. NON-CG-CONTAINING FRAGMENTS)
		non.cg.fragments <- setdiff(unlist(collisions.i), i)		
		controls.meth <- c(controls.meth, calcMeth(SNRs.i, fragments=collisions.i, non.cg.fragments=non.cg.fragments))
	}
	return(as.numeric(controls.meth))
}


calcMeth <- function(SNR.list, fragments=rep(1, length(SNR.list)), non.cg.fragments=numeric(0), method=c("weighted", "proportion"), prune.non.cg.peaks=TRUE, na.rm=FALSE) {
## HANDLE INPUTS TO ENSURE IMPORTANT ASSUMPTIONS (E.G. SNR.LIST REPRESENTS NUMERICAL INPUT)
	method <- match.arg(method)
	SNR.list <- as.numeric(SNR.list)
	num.SNRs <- length(SNR.list)
	non.cg.fragments <- unique(non.cg.fragments)
	if ((num.SNRs > length(fragments)) & (any(SNR.list == 0))) {
		SNR.list <- SNR.list[which(SNR.list != 0)]
		num.SNRs <- length(SNR.list)
	}
	if (num.SNRs != length(fragments)) {
		warning("lengths of input SNR.list and fragments do not match")
		return(NA)
	}
	if (num.SNRs < 1) {
		warning("SNR.list input must contain at least one value")
		return(NA)
	}
## REMOVE NAs FROM INPUT SNR LIST
	if (any(is.na(SNR.list))) {
		if (!na.rm) return(NA)
		if (all(is.na(SNR.list))) return(NA)
		SNR.list[which(is.na(SNR.list))] <- 0
	}
## CALCULATE AVERAGE SNR PER CONTRIBUTING FRAGMENT
	num.fragments <- length(unique(unlist(fragments)))
	SNR.avg <- sum(SNR.list) / num.fragments
## REMOVE NON-CG PEAKS FROM INPUT DATASET VIA ASSUMPTION OF EQUIVALENT AVERAGE PEAK HEIGHTS
	if (prune.non.cg.peaks) {
		for (i in non.cg.fragments) {
			temp <- which(unlist(lapply(lapply(fragments, "%in%", i), any)))
			SNR.list[temp] <- max(SNR.list[temp]-SNR.avg, 0)
		}
		fragments <- lapply(fragments, setdiff, non.cg.fragments)
		num.fragments <- length(unique(unlist(fragments)))
		non.cg.fragments <- numeric(0)
	}
	cgs <- unique(unlist(fragments))
	if (num.fragments < 1) {
		return(NA)	
	}
## SETUP SYSTEM OF LINEAR EQUATIONS
	N <- num.SNRs * num.fragments
	linear.eqs <- matrix(0, nr=0, nc=N)
	solutions <- numeric(0)
## RULE ONE: TOTAL PROPORTIONAL CONTRIBUTIONS FROM EACH CG MUST TOTAL TO 1
	temp <- matrix(0, nr=num.SNRs, nc=N)
	for (i in 1:num.SNRs) {
		cgs.i <- fragments[[i]]	
		index.i <- match(cgs.i, cgs)
		temp[i, (index.i - 1)*num.SNRs + i] <- 1
	}
	linear.eqs <- rbind(linear.eqs, temp)
	solutions <- c(solutions, rep(1, num.SNRs))
## RULE TWO: CGs THAT DO NOT CONTRIBUTE TO A PEAK MUST BE ASSIGNED A COEFFICIENT OF ZERO FOR THAT PEAK
	for (i in 1:num.SNRs) {
		cgs.i <- fragments[[i]]
		missing <- which(!(cgs %in% cgs.i))
		num.missing <- length(missing)
		if (num.missing < 1) next
		temp <- matrix(0, nr=num.missing, nc=N)
		for (j in 1:num.missing) {
			temp[j, (missing[j] - 1)*num.SNRs + i] <- 1
		}
		linear.eqs <- rbind(linear.eqs, temp)
		solutions <- c(solutions, rep(0, num.missing))
	}
## RULE THREE: METHYLATION STATE IS EQUIVALENT FOR TWO (OR MORE) CGs THAT HAVE IDENTICAL PEAK DISTRIBUTIONS
	cg.list <- list()
	for (i in 1:num.fragments) {
		cg.list[[i]] <- which(unlist(lapply(lapply(fragments, match, cgs[i]), any)))
	}
	for (i in which(duplicated(cg.list))) {
		dup.cgs <- which(unlist(lapply(cg.list, identical, cg.list[[i]])))
		dup.cgs <- setdiff(dup.cgs, i)[1]
		temp <- matrix(0, nr=num.SNRs, nc=N)
		for (j in 1:num.SNRs) {
			temp[j, (i - 1)*num.SNRs + j] <- -1
			temp[j, (dup.cgs - 1)*num.SNRs + j] <- 1
		}
		solutions <- c(solutions, rep(0, num.SNRs))
		linear.eqs <- rbind(linear.eqs, temp)
	}
	
## SOLVE SYSTEM OF LINEAR EQUATIONS
	qr.coefs <- qr.coef(qr(linear.eqs), solutions)
## IF SYSTEM OF EQUATIONS IS INCOMPLETE, DETERMINE IDEAL VALUES FOR REMAINING COEFFICIENTS BY RANDOM-VALUE OPTIMIZATION
	if (any(is.na(qr.coefs))) {
		warning("Ambiguous methylation data, calculating optimal values (please be patient as this may take considerable time)", call.=FALSE, immediate.=TRUE)
		na.coefs <- which(is.na(qr.coefs))
		temp <- matrix(0, nr=length(na.coefs), nc=N)
		for (i in 1:length(na.coefs)) {
			temp[i, na.coefs[i]] <- 1
		}
		linear.eqs <- rbind(linear.eqs, temp)
		qr.eqs <- qr(linear.eqs)
		
################################################################################################
## DEFINE OPTIMIZATION FUNCTION FOR SITUATIONS WHERE SYSTEM OF LINEAR EQUATIONS IS INSUFFICIENT
################################################################################################
##
		optimizeCoefficients <- function(coefs=rep(0.5, length(na.coefs)), range=(-500000:500000)/1000000, score.best=Inf, tol=1e-4) {
# DEFINE OPTIMIZATION FUNCTION
			opt.eval <- matrix(0, nr=num.fragments, nc=N)
			for (i in 1:num.fragments) {
				opt.eval[i, (i-1)*num.SNRs + 1:num.SNRs] <- SNR.list
			}
			results <- matrix(NA, nr=length(coefs)+1, nc=length(coefs))
			scores <- rep(NA, length(coefs)+1)
			for (i in 1:(length(coefs)+1)) {
				repeat {
					count <- 100
# IDENTIFY A SET OF COEFFICIENTS THAT SATISFY CRITERIA AS BELOW
					repeat {
						repeat {
# GENERATE RANDOM COEFFICIENTS BY RANDOM-WALK ALGORITHM
							rand <- coefs + sample(range*count/100, size=length(coefs))
# ENSURE THAT STAY WITHIN BOUNDARY LIMITS FOR A COEFFICIENT (0<=X<=1)
							if (all(rand >= 0 & rand <= 1)) break
						}
# ENSURE THAT POINTS DO NOT CONVERGE TOO CLOSE TO PREVIOUS SOLUTIONS						
						if (any(sqrt(apply(as.matrix(results - rand)**2, 1, sum)) < 0.1, na.rm=TRUE)) next
						temp <- c(solutions, rand)
						temp.coefs <- qr.coef(qr.eqs, temp)
# ENSURE THAT STAY WITHIN BOUNDARY LIMITS FOR A COEFFICIENT (0<=X<=1)
						if (any(temp.coefs < 0 & temp.coefs > 1)) next
						score <- sum(abs(rowSums(t(t(opt.eval)*temp.coefs))-SNR.avg))
						if (is.na(score)) next
						if (count <= 0) {
							rand <- coefs
							score <- score.best	
						}
						if (score <= score.best) break
						count <- count-1
					}
					if ((max(abs(range)) < tol | abs(score.best - score) < tol) & count<=90) break
					if (max(abs(range)) > tol & (score == score.best)) break
					coefs <- rand
					range <- range * min(0.9, max(0.5, score/score.best, na.rm=TRUE))
					score.best <- score
					if (max(abs(range)) < tol) range <- range / tol
				}			
				results[i, ] <- coefs
				coefs <- rep(0.5, length(coefs))
				scores[i] <- score.best
				score.best <- Inf
				range <- (-500000:500000)/1000000
			}
			score.worst <- max(scores, na.rm=TRUE)
			slope <- apply(as.matrix(apply(as.matrix(results), 2, diff)), 2, mean)
			results.avg <- apply(as.matrix(results), 2, mean)
			limit.max <- (rep(1, length(rand)) - results.avg)/slope
			which.dim <- which.min(abs(limit.max))
			rand.max <- limit.max[which.dim]*slope + results.avg
			limit.min <- (rep(0, length(rand)) - results.avg)/slope
			which.dim <- which.min(abs(limit.min))
			rand.min <- limit.min[which.dim]*slope + results.avg
			incr <- (rand.max-rand.min)/50
			results <- rep(-1, 51)		
			for (i in 1:51) {
				temp <- c(solutions, rand.min + (i - 1)*incr)
				temp.coefs <- qr.coef(qr.eqs, temp)
				if (any(temp.coefs < 0 & temp.coefs > 1)) next
				score <- sum(abs(rowSums(t(t(opt.eval)*temp.coefs))-SNR.avg))
				results[i] <- score
			}
			results <- which(results <= score.worst + tol*10)
			if (length(results) < 1) {
				return(results.avg)
			}
			rand.max <- rand.min + (results[length(results)] - 1)*incr
			rand.min <- rand.min + (results[1] - 1)*incr
			return((rand.min + rand.max)/2)
		}
##
##
################################################################################################
		
		solutions <- c(solutions, optimizeCoefficients())
		qr.coefs <- qr.coef(qr(linear.eqs), solutions)
	}
	
################################################################################################
## DEFINE METHYLATION CALCULATION SUB-FUNCTION 
################################################################################################
	calcMeth.sub <- function(SNR.list) {
		numerator <- 0
		denominator <- 0
		N <- length(SNR.list)
		SNR.sum <- sum(SNR.list, na.rm=TRUE)	
		for (i in 1:N) {
			SNR.i <- SNR.list[i]
			if (is.na(SNR.i)) next
			denominator <- denominator + SNR.i
			switch(method,
				   "weighted" = numerator <- numerator + (i - 1) * SNR.i / (N - 1),
				   "proportion" = numerator <- numerator + min(1, i - 1) * SNR.i
				   )
		}
		if ((denominator <= 0) | (numerator < 0)) return(NA)
		return(min(1, numerator / denominator))
	}
##
##
################################################################################################
	
	results <- c()
	for (i in 1:num.fragments) {
		results <- c(results, calcMeth.sub(SNR.list[cg.list[[i]]]*qr.coefs[(i - 1)*num.SNRs + cg.list[[i]]]))
	}
	names(results) <- cgs
	return(results)
}


analyzeCpGs <- function(fragments, peaks, method=c("weighted", "proportion")) {
	## HOW MANY CGs TO ANALYZE?
	CpG.num <- sum(unlist(lapply(fragments, slot, "CpGs"))[!unlist(lapply(fragments, slot, "conversion.control"))])
	## STORE CG DATA AS NUMERIC VECTOR
	CpG.data <- rep(NA, CpG.num)
	## IDENTIFY FRAGMENTS THAT HAVE CGs AND ARE NOT CONVERSION CONTROLS
	CpGs <- which((unlist(lapply(fragments, slot, "CpGs")) > 0) & (!unlist(lapply(fragments, slot, "conversion.control"))))
	CpG.fragment.map <- c()
	for (i in CpGs) {
		CpG.fragment.map <- c(CpG.fragment.map, rep(i, fragments[[i]]$CpGs))
	}
	for (i in 1:CpG.num) {
		collisions.i <- CpG.fragment.map[i]
		## FRAGMENT COLLISIONS WITH CpG-CONTAINING FRAGMENTS (IDENTIFY TOTAL SET/STRING OF OVERLAPPING FRAGMENTS)
		repeat {
			temp <- collisions.i
			for (collision in collisions.i) {
				temp <- unique(c(temp, unlist(fragments[[collision]]$collision.IDs)))
			}
			if (identical(collisions.i, temp)) break
			collisions.i <- temp
		}
		## WHAT IS EXPECTED MOLECULAR WEIGHT FOR EACH FRAGMENT?
		MWs <- unique(unlist(lapply(fragments[collisions.i], slot, "MW")))
		MWs <- MWs[order(MWs)]
		## MAP ALL COLLIDING FRAGMENTS TO LIST OF MOLECULAR WEIGHTS
		collisions.i <- lapply(MWs, 
							   function(x) {
								temp <- findFragments(x, fragments)
								IDs <- unlist(lapply(temp, slot, "ID"))
								return(unique(IDs))
							   }
							)
		## IDENTIFY WHICH FRAGMENTS ARE NON-CONTRIBUTORY TO METHYLATION STATE (E.G. NON-CG-CONTAINING FRAGMENTS)
		non.cg.fragments <- setdiff(unlist(collisions.i), CpG.fragment.map)
		## IDENTIFY PEAKS AMONG PEAKLIST THAT MAP TO INPUT MOLECULAR WEIGHTS
		which.peaks.matched <- findPeaks(MWs, unlist(lapply(peaks, slot, "MW.actual")))
		which.peaks.matched.theory <- findPeaks(MWs, unlist(lapply(peaks, slot, "MW.theoretical")))
		which.peaks.matched[which(is.na(which.peaks.matched))] <- which.peaks.matched.theory[which(is.na(which.peaks.matched))]
		## EXTRACT SIGNAL-TO-NOISE RATIO (SNR) INFORMATION FOR EACH PEAK
		SNR.list <- which.peaks.matched
		## HANDLE CASES WHERE UNABLE TO FIND ONE OR MORE PEAKS => INDETERMINATE METHYLATION INFORMATION
		if (any(!is.na(which.peaks.matched))) {
			SNR.list[which(!is.na(which.peaks.matched))] <- unlist(lapply(peaks[which.peaks.matched[which(!is.na(which.peaks.matched))]], slot, "SNR"))
		}
		## CALCULATE PERCENT METHYLATION
		CpG.data[i] <- calcMeth(SNR.list, fragments=collisions.i, non.cg.fragments=non.cg.fragments, method=method, na.rm=TRUE)[as.character(CpG.fragment.map[i])]
	}
	return(CpG.data)
}
