setClass("MassArrayPeak",
	representation(
		ID = "integer",
		MW.theoretical = "numeric",
		MW.actual = "numeric",
		probability = "numeric",
		SNR = "numeric",
		height = "numeric",
		sample.intensity = "numeric",
		ref.intensity = "numeric",
		sequence = "character",
		adduct = "character",
		type = "character",
		charge = "integer",
		collisions = "integer",
		components = "integer",
		missing = "logical",
		new = "logical"
	),
	prototype(
		ID = as.integer(-1),
		MW.theoretical = as.numeric(0),
		MW.actual = as.numeric(0),
		probability = as.numeric(0),
		SNR = as.numeric(0),
		height = as.numeric(0),
		sample.intensity = as.numeric(0),
		ref.intensity = as.numeric(0),
		sequence = "",
		adduct = "",
		type = "",
		charge = as.integer(0),
		collisions = as.integer(0),
		components = as.integer(0),
		missing = FALSE,
		new = FALSE
	)
)


setMethod("initialize",
	"MassArrayPeak",
	function(.Object,
		ID,
		MW.theoretical = as.numeric(0),
		MW.actual = as.numeric(0),
		probability = as.numeric(0),
		SNR = as.numeric(0),
		height = as.numeric(0),
		sample.intensity = as.numeric(0),
		ref.intensity = as.numeric(0),
		sequence = "",
		adduct = c("", "Na", "K"),
		type = c("Expected", "Modified", "Unknown", "Control"),
		charge = as.integer(1),
		collisions = as.integer(0),
		components = as.integer(0),
		missing = FALSE,
		new = FALSE,
		...	
	) {
		.Object@ID <- as.integer(ID)
		.Object@MW.theoretical <- as.numeric(MW.theoretical)
		.Object@MW.actual <- as.numeric(MW.actual)
		.Object@probability <- as.numeric(probability)
		.Object@SNR <- as.numeric(SNR)
		.Object@height <- as.numeric(height)
		.Object@sample.intensity <- as.numeric(sample.intensity)
		.Object@ref.intensity <- as.numeric(ref.intensity)
		.Object@sequence <- as.character(toupper(sequence))
		.Object@adduct <- match.arg(adduct)
		.Object@type <- match.arg(type)
		.Object@charge <- as.integer(charge)
		.Object@collisions <- as.integer(collisions)
		.Object@components <- as.integer(components)
		.Object@missing <- missing
		.Object@new <- new
		return(.Object)
	}
)


setValidity("MassArrayPeak",
	function(object) {
		return(TRUE)
	}
)


setMethod("$", "MassArrayPeak",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)


setMethod("$<-", "MassArrayPeak",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'MassArrayPeak'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)


