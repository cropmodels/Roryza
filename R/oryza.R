

if (!isGeneric("run")) { setGeneric("run", function(x, ...) standardGeneric("run")) }
if (!isGeneric("crop<-")) { setGeneric("crop<-", function(x, value) standardGeneric("crop<-")) }
if (!isGeneric("soil<-")) { setGeneric("soil<-", function(x, value) standardGeneric("soil<-")) }
if (!isGeneric("control<-")) { setGeneric("control<-", function(x, value) standardGeneric("control<-")) }
if (!isGeneric("weather<-")) { setGeneric("weather<-", function(x, value) standardGeneric("weather<-")) }


oryza <- function(crop, weather, soil, control) {
	control$modelstart <- as.integer(as.Date(control$modelstart))
	if (!is.null(control$TMCTB) && is.matrix(control$TMCTB)) {
		control$TMCTB <- as.vector(control$TMCTB)
	}
	d <- .oryza(crop, weather, soil, control)
	date <- as.Date(control$modelstart, origin = "1970-01-01") + (d[, "step"] - 1)
	data.frame(date = date, d)
}


oryza_model <- function(crop, weather, soil, control) {
	m <- OryzaModel$new()
	if (!missing(crop)) { crop(m) <- crop }
	if (!missing(soil)) { soil(m) <- soil }
	if (!missing(control)) { control(m) <- control }
	if (!missing(weather)) { weather(m) <- weather }
	m
}


setMethod("run", signature("Rcpp_OryzaModel"),
	function(x, ...) {
		x$run()
		stopError <- isTRUE(list(...)$stopError)
		msgs <- x$messages
		nm <- length(msgs)
		if (nm > 0) {
			x$messages <- character(0)
			if (x$fatalError) {
				x$fatalError <- FALSE
				errm <- msgs[nm]
				if (nm > 1) {
					warning(paste(msgs[-nm], collapse = "\n"))
				}
				if (stopError) stop(errm) else warning(paste("Error :", errm))
			} else {
				warning(paste(msgs, collapse = "\n"))
			}
		}
		out <- matrix(x$output$values, ncol = length(x$output$names), byrow = TRUE)
		colnames(out) <- x$output$names
		out <- data.frame(out)
		date <- as.Date(x$control$modelstart, origin = "1970-01-01") + (out$step - 1)
		data.frame(date, out)
	}
)


.crop_pars <- c(
	"TBD", "TBLV", "TMD", "TOD", "DVRJ", "DVRI", "DVRP", "DVRR", "MOPP", "PPSE",
	"SHCKD", "COLDMIN", "COLDEAD", "RGRLMX", "RGRLMN", "SHCKL", "SWISLA",
	"ASLA", "BSLA", "CSLA", "DSLA", "SLAMAX", "SLATB", "SSGATB",
	"FRPAR", "SCP", "CO2REF", "CO2", "KDFTB", "KNFTB", "EFFTB", "REDFTT", "NFLVTB",
	"MAINLV", "MAINST", "MAINSO", "MAINRT", "TREF", "Q10",
	"CRGLV", "CRGST", "CRGSO", "CRGRT", "CRGSTR", "LRSTR",
	"FSTR", "TCLSTR", "SPGF", "WGRMX",
	"FSHTB", "FLVTB", "FSTTB", "FSOTB", "DRLVT",
	"FCLV", "FCST", "FCSO", "FCRT", "FCSTR",
	"GZRT", "ZRTMCW", "ZRTMCD",
	"NFLVI", "FNLVI", "NMAXUP", "NMAXSO", "RFNLV", "RFNST", "TCNTRF", "FNTRT",
	"NMAXLT", "NMINLT", "NMINSOT", "NSLLVT",
	"LAPE", "DVSI", "WLVGI", "WSTI", "WRTI", "WSOI", "ZRTI", "ZRTTR",
	"NH", "NPLH", "NPLSB", "NPLDS",
	"ULLS", "LLLS", "ULDL", "LLDL", "ULLE", "LLLE", "ULRT", "LLRT", "SWIRTR"
)

.crop_char_pars <- c("SWISLA", "SWIRTR")

.set_crop_list <- function(x, value) {
	nms <- names(value)
	miss <- setdiff(.crop_pars, nms)
	if (length(miss)) stop(paste("parameters missing:", paste(miss, collapse = ", ")))
	value <- value[.crop_pars]
	for (nm in names(value)) {
		v <- value[[nm]]
		if (nm %in% .crop_char_pars) {
			eval(parse(text = paste0("x$crop$", nm, " <- ", deparse(as.character(v)))))
		} else if (is.matrix(v)) {
			eval(parse(text = paste0("x$crop$", nm, " <- ", paste("c(", paste(as.vector(v), collapse = ","), ")"))))
		} else if (length(v) > 1) {
			eval(parse(text = paste0("x$crop$", nm, " <- ", paste("c(", paste(as.numeric(v), collapse = ","), ")"))))
		} else {
			eval(parse(text = paste0("x$crop$", nm, " <- ", as.numeric(v))))
		}
	}
	x
}

setMethod("crop<-", signature("Rcpp_OryzaModel", "list"),
	function(x, value) .set_crop_list(x, value)
)


.soil_pars <- c("SCODE", "NL", "TKL", "ZRTMS", "WL0MX", "WL0I", "SWITPD", "SWITGW",
	"SWITPF", "SWITVP", "SWITKH", "FIXPERC", "KST", "WCST", "VGA", "VGL", "VGN", "VGR",
	"WCLI", "ZWTB", "RIWCLI")

setMethod("soil<-", signature("Rcpp_OryzaModel", "list"),
	function(x, value) {
		x$setSoil(value)
		x
	}
)


setMethod("weather<-", signature("Rcpp_OryzaModel", "data.frame"),
	function(x, value) {
		parameters <- c("date", "srad", "tmin", "tmax", "prec", "wind", "vapr")
		nms <- colnames(value)
		if (!all(parameters %in% nms)) {
			stop(paste("weather variables missing:", paste(parameters[!(parameters %in% nms)], collapse = ", ")))
		}
		w <- new("Rcpp_Weather")
		w$date <- as.integer(value$date)
		w$srad <- value$srad
		w$tmin <- value$tmin
		w$tmax <- value$tmax
		w$prec <- value$prec
		w$wind <- value$wind
		w$vapr <- value$vapr
		x$wth <- w
		x
	}
)


.req_ctr_pars <- c("modelstart", "cropstart", "max_duration", "water_limited", "latitude", "CO2")
.opt_ctr_pars <- c("output", "ANGSTA", "ANGSTB", "FAOF", "ESTAB", "ETMOD", "SBDUR", "TMPSB", "TMCTB",
                   "RICETYPE", "elevation", "nitrogen_limited", "WATBAL", "SWITIR", "DVSIMAX",
                   "IRRI", "SLMIN", "KPAMIN", "WCMIN", "WL0DAY", "WL0MIN", "RIRRIT", "ISTAGET",
                   "FERTIL", "RECNIT", "SOILSP")

setMethod("control<-", signature("Rcpp_OryzaModel", "list"),
	function(x, value) {
		nms <- names(value)
		if (!all(.req_ctr_pars %in% nms)) {
			stop(paste("parameters missing:", paste(.req_ctr_pars[!(.req_ctr_pars %in% nms)], collapse = ", ")))
		}
		x$control$modelstart <- as.integer(as.Date(value$modelstart))
		x$control$cropstart <- as.integer(value$cropstart)
		x$control$max_duration <- as.integer(value$max_duration)
		x$control$water_limited <- as.logical(value$water_limited)
		x$control$latitude <- as.numeric(value$latitude)
		x$control$CO2 <- as.numeric(value$CO2)
		if (!is.null(value$elevation)) x$control$elevation <- as.numeric(value$elevation)
		if (!is.null(value$output)) x$control$output_option <- as.character(value$output)
		if (!is.null(value$ANGSTA)) x$control$ANGSTA <- as.numeric(value$ANGSTA)
		if (!is.null(value$ANGSTB)) x$control$ANGSTB <- as.numeric(value$ANGSTB)
		if (!is.null(value$FAOF)) x$control$FAOF <- as.numeric(value$FAOF)
		if (!is.null(value$ESTAB)) x$control$ESTAB <- as.character(value$ESTAB)
		if (!is.null(value$ETMOD)) x$control$ETMOD <- as.character(value$ETMOD)
		if (!is.null(value$SBDUR)) x$control$SBDUR <- as.integer(value$SBDUR)
		if (!is.null(value$TMPSB)) x$control$TMPSB <- as.numeric(value$TMPSB)
		if (!is.null(value$RICETYPE)) x$control$RICETYPE <- as.character(value$RICETYPE)
		if (!is.null(value$nitrogen_limited)) x$control$nitrogen_limited <- as.logical(value$nitrogen_limited)
		if (!is.null(value$WATBAL)) x$control$WATBAL <- as.character(value$WATBAL)
		if (!is.null(value$SWITIR)) x$control$SWITIR <- as.integer(value$SWITIR)
		if (!is.null(value$DVSIMAX)) x$control$DVSIMAX <- as.numeric(value$DVSIMAX)
		if (!is.null(value$IRRI)) x$control$IRRI <- as.numeric(value$IRRI)
		if (!is.null(value$SLMIN)) x$control$SLMIN <- as.integer(value$SLMIN)
		if (!is.null(value$KPAMIN)) x$control$KPAMIN <- as.numeric(value$KPAMIN)
		if (!is.null(value$WCMIN)) x$control$WCMIN <- as.numeric(value$WCMIN)
		if (!is.null(value$WL0DAY)) x$control$WL0DAY <- as.integer(value$WL0DAY)
		if (!is.null(value$WL0MIN)) x$control$WL0MIN <- as.numeric(value$WL0MIN)
		if (!is.null(value$TMCTB)) {
			v <- value$TMCTB
			if (is.matrix(v)) v <- as.vector(v)
			x$control$TMCTB <- as.numeric(v)
		}
		if (!is.null(value$RIRRIT)) {
			v <- value$RIRRIT
			if (is.matrix(v)) v <- as.vector(v)
			x$control$RIRRIT <- as.numeric(v)
		}
		if (!is.null(value$ISTAGET)) {
			v <- value$ISTAGET
			if (is.matrix(v)) v <- as.vector(v)
			x$control$ISTAGET <- as.numeric(v)
		}
		if (!is.null(value$FERTIL)) {
			v <- value$FERTIL
			if (is.matrix(v)) v <- as.vector(v)
			x$control$FERTIL <- as.numeric(v)
		}
		if (!is.null(value$RECNIT)) {
			v <- value$RECNIT
			if (is.matrix(v)) v <- as.vector(v)
			x$control$RECNIT <- as.numeric(v)
		}
		if (!is.null(value$SOILSP)) x$control$SOILSP <- as.numeric(value$SOILSP)
		x
	}
)
