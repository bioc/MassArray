plot.MassArrayData <- function(x, y, ...) UseMethod("plot")


plot.MassArrayData <- function (x, ...,  collapse=TRUE, bars=TRUE, scale=TRUE, sequence=TRUE, labels=TRUE, colors=TRUE, main=position(x), width=1.5) {
	if (!validObject(x)) {
		stop("not a valid object of 'MassArrayData'")
	}
	sequence.tagged <- paste(x$fwd.tag, x$sequence, x$rev.tag, sep="")
	N <- nchar(sequence.tagged)
	CGs <- CG.positions <- gregexpr("(CG|YG|CR)", sequence.tagged)[[1]]
	if (length(CGs) < 1) stop("sequence contains no measurable CpG sites")
	## FOR NOW IGNORE GROUPING/SUBSET INFORMATION!!! => ASSUME ALL SAMPLES WITHIN A MassArrayData OBJECT ARE FROM ONE GROUP
	CGs.num <- dim(x$CpG.data)[2]
	sample.rxns <- unlist(lapply(x$samples, slot, "rxn"))
	sample.names <- samples(x)
	rxns <- unique(sample.rxns)
	num.rxns <- length(rxns)
	if (num.rxns < 1) stop("Spectral data contains no reactions")
	## GROUP METHYLATION DATA INTO ONE OR MORE SUBSETS OF SAMPLES
	if (collapse) {
		if (length(x$groups) < 1) x$groups <- rep("", dim(x$CpG.data)[1])
		x$groups[which(is.na(x$groups))] <- ""
		groups <- length(unique(x$groups))
		error.bars <- data.grouped <- matrix(NA, nrow=groups, ncol=CGs.num)
		rxn <- rep(NA, groups)
		for (i in 1:groups) {
			group.i <- which(x$groups == unique(x$groups)[i])
			if (length(group.i) > 1) {
				error.bars[i, ] <- apply(x$CpG.data[group.i, ], 2, mad, na.rm=TRUE)
				data.grouped[i, ] <- apply(x$CpG.data[group.i, ], 2, mean, na.rm=TRUE)
			}
			else {
				error.bars[i, ] <- NA
				data.grouped[i, ] <- x$CpG.data[group.i, ]
			}
			if (length(unique(sample.rxns[group.i])) > 1) {
				warning("group ('", unique(x$groups)[i], "') contains multiple rxn variants")
				rxn[i] <- "T&C"
			}
			else {
				rxn[i] <- sample.rxns[group.i][1]
			}
		}
		sample.rxns <- rxn
		error.bars[which(is.na(error.bars))] <- 0
		x$CpG.data <- data.grouped
		sample.names <- unique(x$groups)
	}
	else {
		error.bars <- rep(NA, CGs.num)
		bars <- FALSE
	}
	samples.num <- dim(x$CpG.data)[1]

	if (scale) {
		position.units <- "bp"
		CGs <- CGs + 0.5
		if (!is.numeric(width)) width <- 1.5
		if (is.na(width) | is.nan(width)) width <- 1.5
		if (width < 0) width <- abs(width)
		width <- min(width, N)
		ticks <- pretty(1:N, n=ceiling(N/10))
		axis.labels <- pretty(1:N, n=ceiling(N/50))
	}
	else {
		CGs <- 1:length(CGs)
		N <- length(CGs)
		position.units <- "CG number"
		sequence <- FALSE
		width <- 0.25
		ticks <- 1:N
		axis.labels <- 1:N
	}
	## BEGIN PLOTTING
	plot(0, type="n", 
		xlim=c(1, N), 
		ylim=c(0, 10), xlab=paste("Nucleotide position (", position.units, ")", sep=""), ylab="",
		main=main, axes=FALSE)
	axis(1, 
		labels=axis.labels[match(ticks, axis.labels)], 
		at=ticks, 
		cex.lab=0.75, cex.axis=0.75, las=2
	)
	## DEFINE PLOTTING WINDOW AND SCALING FUNCTION
	usr <- par("usr")
	scale <- function (x, min1, max1, min2=usr[1], max2=usr[2]) {
		return(min2 + (x - min1) * (max2 - min2) / (max1 - min1))
	}
	
	if (labels) {
		text(CGs, 
			rep(usr[4], CGs.num), 
			labels=1:CGs.num, 
			cex=0.4, col="black", pos=1
		)
		if (sequence) {
			text(1:N, 
				rep(usr[3], N),
				labels=unlist(strsplit(sequence.tagged, "")),
				cex=0.2, col="darkgray", pos=3
			)
#			text((N + 1)/2, usr[3] + 0.035*(usr[4] - usr[3]), labels=position(x), cex=0.25, col="darkgray")
		}
	}
	else {
		colors <- FALSE	
	}
	for (i in 1:samples.num) {
		if (labels) {
			text(usr[1],
				scale(samples.num - i + 1.5, 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]),
				labels=sample.names[i],
				cex=0.4, col="black", pos=4, srt=90
			)	
			
		}
		switch(sample.rxns[i],
			"T" = fragments <- x$fragments.T,
			"C" = fragments <- x$fragments.C,
			fragments <- c(x$fragments.T, x$fragments.C)
		)
		positions <- unlist(lapply(fragments, slot, "position"))
		lengths <- unlist(lapply(fragments, slot, "length"))
		cg.fragments <- lapply(CG.positions, 
			function (x) {
				return(which(x >= positions & x < positions + lengths))
			}
		)
		for (j in 1:CGs.num) {
			fragment <- fragments[cg.fragments[[j]]]
			shading.density <- NA
			lty <- "solid"
			label.color <- "black"
			label.text <- NA
			if (colors) {
				color.cg <- color.border <- "black"
				color.bg <- "peachpuff"
				if (all(unlist(lapply(fragment, slot, "CpGs")) > 1)) {
					color.bg <- "#FFFFCC"
					color.border <- "peachpuff"
#					lty <- "dotted"
				}
				if (all(unlist(lapply(unlist(lapply(fragment, slot, "CG.collision.IDs"), recursive=FALSE), length)) > 0)) {
					color.bg <- "#FF9900"
				}
			}
			else {
				color.bg <- "#F9F9F9"
				color.cg <- color.border <- "black"
			}
			if (is.na(x$CpG.data[i, j])) {
				color.cg <- color.border <- "#CCCCCC"
				color.bg <- "white"
				shading.density <- NA
#				x$CpG.data[i, j] <- 0
			}
			else if (!any(unlist(lapply(fragment, slot, "assayable")))) {
				color.cg <- color.border <- "#BBBBBB"
				color.bg <- "#EEEEEE"
			}
#			else if (colors) {
#				color.cg <- color.border <- "black"
#				color.bg <- "peachpuff"
#			}
#			else {
#				color.bg <- "#F9F9F9"
#				color.cg <- color.border <- "black"
#			}
			if (all(unlist(lapply(fragment, slot, "required")))) {
#				color.border <- "#DD3333"
				label.text <- "*"
				if (colors) label.color <- "red"
			}
			if (labels & !is.na(label.text)) {
				text(CGs[j], 
					0.98*usr[4], 
					labels=label.text, 
					cex=0.6, col=label.color, pos=1
				)
			}
			rect(CGs[j] - width, 
				scale(samples.num - i + 1.05, 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]), 
				CGs[j] + width, 
				scale(samples.num - i + 1.95, 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]), 
				border=color.border, col=color.bg, lwd=1, lty=lty, density=shading.density
			)
			rect(CGs[j] - width, 
				scale(samples.num - i + 1.05, 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]), 
				CGs[j] + width, 
				scale(samples.num - i + scale(x$CpG.data[i, j], 0, 1, 1.05, 1.95), 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]), 
				border=NA, col=color.cg, lwd=1, lty=lty, density=shading.density
			)
			## DISPLAY ERROR BARS (IF DESIRED)
			if (bars) {
				color.cg <- as.numeric(col2rgb(color.cg))+(as.numeric(col2rgb(color.bg))-as.numeric(col2rgb(color.cg)))/2
				color.cg <- rgb(color.cg[1], color.cg[2], color.cg[3], maxColorValue=255)
				cg <- x$CpG.data[i, j] - c(error.bars[i, j], -error.bars[i, j])
				## SKIP ERROR BAR UNLESS THERE IS NON-ZERO DATA TO PLOT
				if (is.na(x$CpG.data[i, j]) | is.na(error.bars[i, j]) | (error.bars[i, j] <= 0)) next
				## DRAW VERTICAL PORTION OF ERROR BAR
				lines(rep(CGs[j], 2), 
					scale(samples.num - i + scale(c(max(0, cg[1]), min(1, cg[2])), 0, 1, 1.05, 1.95), 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]), 
					col=color.cg, lwd=1)
				## DRAW LOWER DELIMITER OF ERROR BAR
				if (cg[1] > 0) {
					lines(c(CGs[j] - width/3, CGs[j] + width/3), 
						scale(samples.num - i + scale(rep(cg[1], 2), 0, 1, 1.05, 1.95), 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]), 
						col=color.cg, lwd=1)
				}
				## DRAW UPPER DELIMITER OF ERROR BAR
				if (cg[2] < 1) {
					lines(c(CGs[j] - width/3, CGs[j] + width/3), 
						scale(samples.num - i + scale(rep(cg[2], 2), 0, 1, 1.05, 1.95), 1, samples.num + 1, usr[3] + 0.05*(usr[4] - usr[3]), 0.95*usr[4]), 
						col=color.cg, lwd=1)
				}
			}
			
			
		}
				
	}
}