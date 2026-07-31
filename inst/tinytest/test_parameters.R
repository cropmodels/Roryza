# Parameter loaders

crops <- oryza_crop()
expect_true("IR72" %in% crops)

crop <- oryza_crop("IR72")
expect_true(is.list(crop))
expect_true(all(c("TBD", "DVRJ", "DVRR", "RGRLMX", "SWISLA") %in% names(crop)))
expect_equal(crop$TBD, 8)
expect_equal(crop$SWISLA, "FUNCTION")

soils <- oryza_soil()
expect_true(all(c("paddy", "potential") %in% soils))

paddy <- oryza_soil("paddy")
expect_true(is.list(paddy))
expect_true("NL" %in% names(paddy) || "TKL" %in% names(paddy))

ctrl <- oryza_control()
expect_true(is.list(ctrl))
expect_true(all(c("modelstart", "cropstart", "max_duration", "water_limited", "latitude") %in% names(ctrl)))
expect_false(isTRUE(ctrl$water_limited))
expect_true(inherits(ctrl$modelstart, "Date") || is.numeric(ctrl$modelstart) || is.character(ctrl$modelstart))
