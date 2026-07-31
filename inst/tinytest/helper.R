# Shared helpers for Roryza tinytests (sourced from each test_*.R)

library(Roryza)

.oryza_synth_weather <- function(n, start = as.Date("1992-01-01"), seed = 1L) {
	dates <- start + seq_len(n) - 1L
	set.seed(seed)
	data.frame(
		date = dates,
		srad = 18000 + stats::rnorm(n, 0, 1000),
		tmin = 22 + 2 * sin(2 * pi * seq_len(n) / 365),
		tmax = 31 + 2 * sin(2 * pi * seq_len(n) / 365),
		prec = pmax(0, stats::rnorm(n, 5, 8)),
		wind = 2 + stats::runif(n),
		vapr = 2.5 + 0.2 * stats::runif(n)
	)
}

.oryza3_dev_exe <- function() {
	exe <- if (.Platform$OS.type == "windows") "oryza3.exe" else "oryza3"
	# tinytest cwd is inst/tinytest (source or installed)
	cands <- c(
		file.path("..", "..", "dev", "inst", exe),
		file.path("..", "..", "dev", "oryza", "InputData", exe)
	)
	for (p in cands) {
		if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
	}
	NA_character_
}

.oryza3_templates <- function() {
	cands <- c(
		file.path("..", "..", "dev", "oryza2000", "templates"),
		file.path("..", "..", "dev", "inst", "input")
	)
	for (p in cands) {
		if (dir.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
	}
	NA_character_
}
