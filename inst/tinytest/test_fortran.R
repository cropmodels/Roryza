# FORTRAN/FSE reference engine (dev/ only). Skipped unless at_home and binary present.

source("helper.R")

if (!isTRUE(tinytest::at_home())) {
	exit_file("FORTRAN FSE tests run only at_home")
}

exe <- .oryza3_dev_exe()
tmpl <- .oryza3_templates()
if (is.na(exe) || is.na(tmpl)) {
	exit_file("dev/inst/oryza3.exe or templates not found — run Rscript dev/tools/build_oryza3.R")
}

.parse_res <- function(f) {
	lines <- readLines(f, warn = FALSE)
	i <- grep("^TIME", lines)[1]
	hdr <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
	body <- lines[(i + 1L):length(lines)]
	body <- body[nzchar(trimws(body))]
	dat <- utils::read.delim(textConnection(body), header = FALSE, sep = "\t",
		strip.white = TRUE, na.strings = c("-", ""), check.names = FALSE)
	names(dat) <- hdr
	for (n in names(dat)) dat[[n]] <- suppressWarnings(as.numeric(dat[[n]]))
	dat
}

.set_char <- function(lines, name, value) {
	pat <- paste0("^\\s*", name, "\\s*=")
	hit <- which(!grepl("^\\s*\\*", lines) & grepl(pat, lines))
	asn <- sprintf("%s = '%s'", name, value)
	if (length(hit)) {
		lines[hit[1]] <- asn
		if (length(hit) > 1) for (i in hit[-1]) lines[i] <- paste0("*", lines[i])
		lines
	} else {
		c(lines, asn)
	}
}

.prepare <- function(watbal, water_limited = TRUE, nitrogen_limited = FALSE, prdel = 1) {
	wd <- tempfile("oryza3_")
	dir.create(wd)
	ok <- file.copy(list.files(tmpl, full.names = TRUE), wd, overwrite = TRUE)
	expect_true(all(ok))
	soil <- switch(toupper(watbal),
		PADDY = "PADDYIN.DAT", SAWAH = "SAWAIN.DAT", SAHEL = "SAHELIN.DAT",
		LOWBAL = "LOWBALIN.DAT", SOILPF = "SOILPFIN.DAT",
		stop("bad watbal"))
	ctrl <- readLines(file.path(wd, "control.dat"), warn = FALSE)
	ctrl <- sub("PRDEL\\s*=\\s*[0-9.]+", sprintf("PRDEL  = %.6f", prdel), ctrl)
	ctrl <- vapply(ctrl, function(ln) {
		if (!grepl("^\\s*\\*", ln) && grepl("^\\s*FILEI2\\s*=", ln)) paste0("*", ln) else ln
	}, character(1), USE.NAMES = FALSE)
	ctrl <- c(ctrl, sprintf("   FILEI2 = '%s'", soil))
	writeLines(ctrl, file.path(wd, "control.dat"), useBytes = TRUE)

	exp <- readLines(file.path(wd, "experiment.dat"), warn = FALSE)
	if (water_limited) {
		exp <- .set_char(exp, "PRODENV", "WATER BALANCE")
		exp <- .set_char(exp, "WATBAL", toupper(watbal))
	} else {
		exp <- .set_char(exp, "PRODENV", "POTENTIAL")
	}
	exp <- .set_char(exp, "NITROENV", if (nitrogen_limited) "NITROGEN BALANCE" else "POTENTIAL")
	writeLines(exp, file.path(wd, "experiment.dat"), useBytes = TRUE)
	wd
}

.run <- function(wd) {
	owd <- setwd(wd)
	on.exit(setwd(owd), add = TRUE)
	out <- system2(exe, stdout = TRUE, stderr = TRUE)
	status <- attr(out, "status")
	expect_true(is.null(status) || identical(as.integer(status), 0L))
	expect_true(file.exists("RES.DAT") && file.info("RES.DAT")$size > 0)
	.parse_res("RES.DAT")
}

.final_wso <- function(res) {
	utils::tail(res$WSO[is.finite(res$WSO)], 1)
}

# Anchors from rebuilt oryza3.exe on packaged IRRI 1992 templates (PRDEL=1)
cases <- list(
	list(name = "POTENTIAL", watbal = "PADDY", wl = FALSE, nl = FALSE, wso = 9648.8),
	list(name = "PADDY", watbal = "PADDY", wl = TRUE, nl = FALSE, wso = 7456.9),
	list(name = "SAWAH", watbal = "SAWAH", wl = TRUE, nl = FALSE, wso = 9645.3),
	list(name = "SAHEL", watbal = "SAHEL", wl = TRUE, nl = FALSE, wso = 6916.3),
	list(name = "LOWBAL", watbal = "LOWBAL", wl = TRUE, nl = FALSE, wso = 9648.8),
	list(name = "SOILPF", watbal = "SOILPF", wl = TRUE, nl = FALSE, wso = 2654.3)
)

for (sc in cases) {
	wd <- .prepare(sc$watbal, sc$wl, sc$nl, prdel = 1)
	res <- .run(wd)
	expect_true(nrow(res) > 50, info = sc$name)
	expect_equal(.final_wso(res), sc$wso, tolerance = 0.05, info = sc$name)
	expect_equal(max(res$DVS, na.rm = TRUE), 2.012, tolerance = 0.001, info = sc$name)
	unlink(wd, recursive = TRUE)
}

# Repeatability: identical RES for two PADDY runs
wd1 <- .prepare("PADDY", TRUE, FALSE, 1)
r1 <- .run(wd1)
unlink(wd1, recursive = TRUE)
wd2 <- .prepare("PADDY", TRUE, FALSE, 1)
r2 <- .run(wd2)
unlink(wd2, recursive = TRUE)
vars <- intersect(c("DVS", "LAI", "WSO", "WAGT", "TRW", "WL0"), names(r1))
n <- min(nrow(r1), nrow(r2))
maxabs <- max(vapply(vars, function(v) {
	ok <- is.finite(r1[[v]][seq_len(n)]) & is.finite(r2[[v]][seq_len(n)])
	if (!any(ok)) 0 else max(abs(r1[[v]][seq_len(n)][ok] - r2[[v]][seq_len(n)][ok]))
}, numeric(1)))
expect_equal(maxabs, 0)
