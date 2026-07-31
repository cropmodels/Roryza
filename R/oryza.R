

if (!isGeneric("run")) { setGeneric("run", function(x, ...) standardGeneric("run")) }
if (!isGeneric("crop<-")) { setGeneric("crop<-", function(x, value) standardGeneric("crop<-")) }
if (!isGeneric("soil<-")) { setGeneric("soil<-", function(x, value) standardGeneric("soil<-")) }
if (!isGeneric("control<-")) { setGeneric("control<-", function(x, value) standardGeneric("control<-")) }
if (!isGeneric("weather<-")) { setGeneric("weather<-", function(x, value) standardGeneric("weather<-")) }


oryza <- function(crop, weather, soil, control, engine = c("fse", "cpp")) {
	engine <- match.arg(engine)
	if (engine == "fse") {
		return(.oryza_via_fse(crop, weather, soil, control))
	}
	control$modelstart <- as.integer(as.Date(control$modelstart))
	if (!is.null(control$TMCTB) && is.matrix(control$TMCTB)) {
		control$TMCTB <- as.vector(control$TMCTB)
	}
	d <- .oryza(crop, weather, soil, control)
	date <- as.Date(control$modelstart, origin = "1970-01-01") + (d[, "step"] - 1)
	data.frame(date = date, d)
}


.oryza_via_fse <- function(crop, weather, soil, control) {
	wl <- isTRUE(as.logical(control$water_limited))
	nl <- isTRUE(as.logical(control$nitrogen_limited))
	watbal <- if (!is.null(control$WATBAL)) toupper(as.character(control$WATBAL)) else "PADDY"
	# TTUTIL distinguishes INTEGER vs REAL by token form ("2" vs "2.")
	.fse_real <- function(x) sprintf("%.6f", as.numeric(x))
	overrides <- list()
	if (!is.null(control$ESTAB)) overrides$ESTAB <- sprintf("'%s'", control$ESTAB)
	if (!is.null(control$ETMOD)) overrides$ETMOD <- sprintf("'%s'", control$ETMOD)
	if (!is.null(control$RICETYPE)) overrides$RICETYPE <- sprintf("'%s'", control$RICETYPE)
	if (!is.null(control$SWITIR)) overrides$SWITIR <- as.integer(control$SWITIR)
	if (!is.null(control$IRRI)) overrides$IRRI <- .fse_real(control$IRRI)
	if (!is.null(control$KPAMIN)) overrides$KPAMIN <- .fse_real(control$KPAMIN)
	if (!is.null(control$DVSIMAX)) overrides$DVSIMAX <- .fse_real(control$DVSIMAX)
	if (!is.null(control$WCMIN)) overrides$WCMIN <- .fse_real(control$WCMIN)
	if (!is.null(control$WL0MIN)) overrides$WL0MIN <- .fse_real(control$WL0MIN)
	if (!is.null(control$SBDUR)) overrides$SBDUR <- as.integer(control$SBDUR)
	# When N-limited and fertilizer not supplied, use the classic IRRI 225 kg N schedule
	if (nl && is.null(control$FERTIL)) {
		overrides$FERTIL <- c(
			"0.", "0.", "1.", "0.", "11.", "0.", "12.", "60.", "13.", "0.",
			"29.", "0.", "30.", "60.", "31.", "0.", "66.", "0.", "67.", "60.",
			"68.", "0.", "94.", "0.", "95.", "45.", "96.", "0.", "366.", "0."
		)
	} else if (!is.null(control$FERTIL)) {
		v <- control$FERTIL
		if (is.matrix(v)) v <- as.vector(v)
		overrides$FERTIL <- as.character(as.numeric(v))
	}
	if (!is.null(control$SOILSP)) overrides$SOILSP <- .fse_real(control$SOILSP)

	prdel <- if (!is.null(control$PRDEL)) as.numeric(control$PRDEL) else 1
	workdir <- oryza_fse_prepare(watbal = watbal, nitrogen_limited = nl,
		water_limited = wl, prdel = prdel, overrides = overrides)
	on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

	if (!missing(weather) && !is.null(weather)) {
		ms <- as.Date(control$modelstart)
		yr <- as.integer(format(ms, "%Y"))
		wfile <- file.path(workdir, sprintf("PHIL1.%03d", yr %% 1000L))
		lon <- if (!is.null(control$longitude)) control$longitude else 121.25
		lat <- if (!is.null(control$latitude)) control$latitude else 14.18
		elev <- if (!is.null(control$elevation)) control$elevation else 21
		oryza_write_weather(weather, wfile, longitude = lon, latitude = lat, elevation = elev)
		exp <- readLines(file.path(workdir, "experiment.dat"), warn = FALSE)
		exp <- .fse_set_raw(exp, "IYEAR", yr)
		exp <- .fse_set_raw(exp, "STTIME", as.numeric(format(ms, "%j")))
		writeLines(exp, file.path(workdir, "experiment.dat"), useBytes = TRUE)
	}

	# crop.dat: keep packaged IR72 template (matches oryza_crop("IR72"))
	# soil: selected via WATBAL in oryza_fse_prepare; soil list reserved for future DAT writer
	invisible(crop); invisible(soil)

	res <- oryza_fse(workdir, quiet = TRUE)
	.oryza_res_to_api(res, control)
}


.oryza_res_to_api <- function(res, control) {
	ms <- as.Date(control$modelstart)
	# TIME in ORYZA is days since STTIME on the timer; for template STTIME=1, TIME==DOY when IYEAR matches
	if ("DOY" %in% names(res) && "YEAR" %in% names(res)) {
		date <- as.Date(paste(as.integer(res$YEAR), as.integer(res$DOY)), format = "%Y %j")
	} else {
		date <- ms + (res$TIME - 1)
	}
	keep <- c("DVS", "LAI", "WAGT", "WST", "WLVG", "WLVD", "WLV", "WSO", "WRR14",
		"TRC", "TRW", "EVSC", "WL0", "IR", "MSKPA1", "CROPSTA", "PCEW", "LESTRS",
		"LRSTRS", "NFLV", "ZRT", "RAIN", "RAINCU")
	keep <- intersect(keep, names(res))
	out <- data.frame(date = date, step = as.integer(res$TIME), res[keep], check.names = FALSE)
	out
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
	"NFLVI", "NMAXLT",
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
                   "IRRI", "SLMIN", "KPAMIN", "WCMIN", "WL0DAY", "WL0MIN", "RIRRIT", "ISTAGET")

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
		x
	}
)
