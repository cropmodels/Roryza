
setMethod("predict", signature("Rcpp_OryzaModel"),
function(object, weather, mstart, soils = NULL, soiltypes = NULL,
         filename = "", overwrite = FALSE, ...) {

	if (!requireNamespace("terra", quietly = TRUE)) {
		stop("package 'terra' is required for spatial predict")
	}
	stopifnot(inherits(weather, "SpatRasterDataset"))
	stopifnot(inherits(mstart, "Date"))
	if (any(is.na(mstart))) {
		stop("mstart cannot be NA")
	}

	watlim <- isTRUE(object$control$water_limited)

	# latitude / elevation are always needed for radiation and ET
	need_soil_layers <- c("latitude", "elevation")
	if (watlim) {
		need_soil_layers <- c("soil", "soildepth", need_soil_layers)
		needed <- c("tmin", "tmax", "srad", "prec", "vapr", "wind")
		if (is.null(soils) || is.null(soiltypes)) {
			stop("water-limited predict requires 'soils' SpatRaster and 'soiltypes' list")
		}
		stopifnot(inherits(soils, "SpatRaster"))
		stopifnot(inherits(soiltypes, "list"), length(soiltypes) > 0)
		miss <- setdiff(need_soil_layers, names(soils))
		if (length(miss)) {
			stop(paste("soils missing layer(s):", paste(miss, collapse = ", ")))
		}
		terra::compareGeom(weather[1], soils, lyrs = FALSE)
		scol <- .makeSoilCollection(soiltypes)
	} else {
		needed <- c("tmin", "tmax", "srad")
		if (is.null(soils)) {
			stop("predict requires a 'soils' SpatRaster with latitude and elevation layers")
		}
		stopifnot(inherits(soils, "SpatRaster"))
		miss <- setdiff(need_soil_layers, names(soils))
		if (length(miss)) {
			stop(paste("soils missing layer(s):", paste(miss, collapse = ", ")))
		}
		terra::compareGeom(weather[1], soils, lyrs = FALSE)
		if (is.null(soiltypes)) {
			soiltypes <- list(oryza_soil("potential"))
		}
		scol <- .makeSoilCollection(soiltypes)
	}

	nms <- names(weather)
	if (!all(needed %in% nms)) {
		stop(paste("missing these weather variables:",
			paste(needed[!(needed %in% nms)], collapse = ", ")))
	}
	weather <- weather[needed]
	nss <- terra::nlyr(weather)
	if (!all(nss == nss[1])) {
		stop("all weather subdatasets must have the same number of layers")
	}

	if (!terra::has.time(weather$tmin)) {
		stop("weather$tmin has no time stamps (dates)")
	}
	dates <- as.Date(terra::time(weather$tmin))
	stopifnot(length(dates) == terra::nlyr(weather$tmin))
	if (any(is.na(dates))) {
		stop("NA in weather dates not allowed")
	}
	dates_i <- as.integer(dates)
	mstart_i <- as.integer(mstart)

	rout <- terra::rast(weather)
	terra::nlyr(rout) <- length(mstart)
	terra::time(rout) <- mstart

	use_raster <- FALSE
	src1 <- tryCatch(unlist(terra::sources(weather[1])[1]), error = function(e) "")
	if (length(src1) && substr(src1[1], 1, 6) == "NETCDF") {
		if (!requireNamespace("raster", quietly = TRUE)) {
			stop("package 'raster' is required to read NETCDF weather sources")
		}
		use_raster <- TRUE
		p <- src1[1]
		f <- unlist(strsplit(gsub("NETCDF:\"", "", p), "\""))[1]
		weather <- lapply(needed, function(i) raster::brick(f, varname = i))
		names(weather) <- needed
	}

	nc <- ncol(rout)
	nr <- nrow(rout)
	if (!use_raster) terra::readStart(weather)
	terra::readStart(soils)

	wopt <- list(...)
	if (is.null(wopt$names)) wopt$names <- as.character(mstart)
	terra::writeStart(rout, filename, overwrite, wopt = wopt)

	# row-by-row (same strategy as Rwofost)
	b <- list(row = seq_len(nr), nrows = rep(1L, nr), n = nr)

	for (i in seq_len(b$n)) {
		if (use_raster) {
			tmin <- as.vector(t(raster::getValues(weather$tmin, b$row[i], b$nrows[i])))
			tmax <- as.vector(t(raster::getValues(weather$tmax, b$row[i], b$nrows[i])))
			srad <- as.vector(t(raster::getValues(weather$srad, b$row[i], b$nrows[i])))
			if (watlim) {
				prec <- as.vector(t(raster::getValues(weather$prec, b$row[i], b$nrows[i])))
				vapr <- as.vector(t(raster::getValues(weather$vapr, b$row[i], b$nrows[i])))
				wind <- as.vector(t(raster::getValues(weather$wind, b$row[i], b$nrows[i])))
			}
		} else {
			tmin <- as.vector(t(terra::readValues(weather$tmin, b$row[i], b$nrows[i], 1, nc, mat = TRUE)))
			tmax <- as.vector(t(terra::readValues(weather$tmax, b$row[i], b$nrows[i], 1, nc, mat = TRUE)))
			srad <- as.vector(t(terra::readValues(weather$srad, b$row[i], b$nrows[i], 1, nc, mat = TRUE)))
			if (watlim) {
				prec <- as.vector(t(terra::readValues(weather$prec, b$row[i], b$nrows[i], 1, nc, mat = TRUE)))
				vapr <- as.vector(t(terra::readValues(weather$vapr, b$row[i], b$nrows[i], 1, nc, mat = TRUE)))
				wind <- as.vector(t(terra::readValues(weather$wind, b$row[i], b$nrows[i], 1, nc, mat = TRUE)))
			}
		}

		ncell_row <- length(tmin) / length(dates_i)
		if (!watlim) {
			zeros <- rep(0, length(tmin))
			prec <- zeros
			vapr <- zeros
			wind <- zeros
		}

		elv <- as.vector(terra::readValues(soils$elevation, b$row[i], b$nrows[i], 1, nc))
		lat <- as.vector(terra::readValues(soils$latitude, b$row[i], b$nrows[i], 1, nc))

		if (watlim) {
			sidx <- as.vector(terra::readValues(soils$soil, b$row[i], b$nrows[i], 1, nc))
			sidx[is.na(sidx)] <- -99
			sidx <- as.integer(sidx)
			depth <- as.vector(terra::readValues(soils$soildepth, b$row[i], b$nrows[i], 1, nc))
			depth[is.na(depth)] <- -99
		} else {
			sidx <- as.integer(rep(1L, ncell_row))
			depth <- rep(-99, ncell_row)
		}

		y <- object$run_batch(tmin, tmax, srad, prec, vapr, wind,
			dates_i, mstart_i, sidx, scol, depth, elv, lat)
		terra::writeValues(rout, round(y), b$row[i], b$nrows[i])
	}

	if (!use_raster) terra::readStop(weather)
	terra::readStop(soils)
	terra::writeStop(rout)
}
)
