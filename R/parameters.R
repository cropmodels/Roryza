

.getNumLst <- function(ini, make_matrix = TRUE) {
	v <- ini[, 3]
	vv <- sapply(v, function(i) strsplit(i, ","), USE.NAMES = FALSE)
	# character parameters kept as character
	char_pars <- c("SWISLA", "SWIRTR", "ESTAB", "ETMOD", "RICETYPE", "PRODENV",
	               "NITROENV", "WATBAL", "RUNMODE", "SCODE", "RIWCLI")
	# 1-d layer vectors (not x-y AFGEN tables)
	vector_pars <- c("TKL", "KST", "WCST", "VGA", "VGL", "VGN", "VGR", "PN",
	                 "WCLI", "WCFC", "WCWP", "WCAD", "WCSTRP", "ISTAGET")
	out <- lapply(seq_along(vv), function(i) {
		nm <- ini[i, 2]
		if (nm %in% char_pars) {
			gsub("'", "", trimws(vv[[i]][1]))
		} else {
			x <- as.numeric(vv[[i]])
			if (make_matrix && length(x) > 1 && !(nm %in% vector_pars)) {
				matrix(x, nrow = 2)
			} else {
				x
			}
		}
	})
	names(out) <- ini[, 2]
	out
}


.notavailable <- function(group, error = TRUE) {
	if (group == "crop") {
		f <- list.files(system.file("oryza/crop", package = "Roryza"), pattern = "\\.ini$", full.names = TRUE)
	} else if (group == "soil") {
		f <- list.files(system.file("oryza/soil", package = "Roryza"), pattern = "\\.ini$", full.names = TRUE)
	} else {
		f <- character(0)
	}
	x <- gsub("\\.ini$", "", basename(f))
	if (error) {
		stop(paste(group, "not available. Choose one of:\n"), paste(x, collapse = ", "), "\n")
	}
	x
}


oryza_control <- function(filename = "") {
	x <- list()
	filename <- trimws(filename)
	if (filename == "") {
		filename <- system.file("oryza/control.ini", package = "Roryza")
	}
	ini <- .readIniFile(filename)

	s <- which(ini[, 2] == "modelstart")
	if (length(s) > 0) {
		x$modelstart <- as.Date(unname(ini[s[1], 3]))
		ini <- ini[-s, , drop = FALSE]
	}
	lst <- .getNumLst(ini, make_matrix = FALSE)
	# logicals
	for (nm in c("water_limited", "nitrogen_limited")) {
		if (!is.null(lst[[nm]])) lst[[nm]] <- as.logical(as.numeric(lst[[nm]]))
	}
	append(x, lst)
}


oryza_soil <- function(name = "") {
	if (missing(name) || trimws(name) == "") {
		return(.notavailable("soil", FALSE))
	}
	name <- trimws(name)
	if (file.exists(name)) {
		return(.getNumLst(.readIniFile(name), make_matrix = FALSE))
	}
	f <- list.files(system.file("oryza/soil", package = "Roryza"), pattern = "\\.ini$", full.names = TRUE)
	soils <- gsub("\\.ini$", "", basename(f))
	if (name %in% soils) {
		return(.getNumLst(.readIniFile(f[which(name == soils)[1]]), make_matrix = FALSE))
	}
	.notavailable("soil")
}


oryza_crop <- function(name = "") {
	if (missing(name) || trimws(name) == "") {
		return(.notavailable("crop", FALSE))
	}
	name <- trimws(name)
	if (file.exists(name)) {
		return(.getNumLst(.readIniFile(name)))
	}
	f <- list.files(system.file("oryza/crop", package = "Roryza"), pattern = "\\.ini$", full.names = TRUE)
	crops <- gsub("\\.ini$", "", basename(f))
	if (name %in% crops) {
		ini <- .readIniFile(f[which(name == crops)[1]])
		ini <- ini[ini[, 2] %in% .crop_pars, , drop = FALSE]
		j <- .crop_pars %in% ini[, 2]
		if (!all(j)) {
			warning(paste("missing parameter(s):", paste(.crop_pars[!j], collapse = ", ")))
		}
		return(.getNumLst(ini))
	}
	.notavailable("crop")
}
