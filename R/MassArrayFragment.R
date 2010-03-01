setClass("MassArrayFragment",
	representation(
		ID = "integer",
		assay.name = "character",
		name = "character",
		sequence = "character",
		position = "integer",
		length = "integer",
		CpGs = "integer",
		MW = "numeric",
		collisions = "integer",
		collision.IDs = "list",
		CG.collisions = "integer",
		CG.collision.IDs = "list",
		type = "character",
		direction = "character",
		extra = "character",
		bisulfite.converted = "logical",
		assayable = "logical",
		conversion.control = "logical",
		required = "logical",
		ignored = "logical",
		primer = "logical"
	),
	prototype(
		ID = as.integer(-1),
		assay.name = "",
		name = "",
		sequence = "",
		position = as.integer(-1),
		length = as.integer(0),
		CpGs = as.integer(0),
		MW = as.numeric(0),
		collisions = as.integer(0),
		collision.IDs = as.list(NULL),
		CG.collisions = as.integer(0),
		CG.collision.IDs = as.list(NULL),
		type = "T",
		direction = "+",
		extra = "5OH",
		bisulfite.converted = FALSE,
		assayable = FALSE,
		conversion.control = FALSE,
		required = FALSE,
		ignored = FALSE,
		primer = FALSE
	)
)


setMethod("initialize",
	"MassArrayFragment",
	function(.Object,
		ID,
		sequence,
		assay.name = "",
		name = "",
		position = -1,
		type = c("T", "C"),
		direction = c("+", "-"),
		extra = c("5OH", "5PPP-3P", "5PPP-3OH"),
		bisulfite.converted = FALSE,
		assayable = FALSE,
		primer = FALSE,
		...	
	) {
		.Object@ID <- as.integer(ID)
		.Object@assay.name <- as.character(assay.name)
		.Object@name <- as.character(name)
		.Object@sequence <- gsub("[^ATCGYRN()]", "", toupper(sequence))
		.Object@position <- as.integer(position)
		.Object@type <- match.arg(type)
		.Object@direction <- match.arg(direction)
		.Object@extra <- match.arg(extra)
		.Object@bisulfite.converted <- bisulfite.converted
		.Object@assayable <- assayable
		.Object@primer <- primer
		.Object@length <- nchar(.Object@sequence)
		.Object@CpGs <- as.integer(countCGs(.Object@sequence))
		.Object@MW <- calcMW(.Object@sequence, rxn=as.character(.Object@type))
		return(.Object)
	}
)


setValidity("MassArrayFragment",
	function(object) {
		return(TRUE)
	}
)


setMethod("$", "MassArrayFragment",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)


setMethod("$<-", "MassArrayFragment",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'MassArrayFragment'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)


