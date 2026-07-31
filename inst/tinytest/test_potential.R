# Potential production (C++ engine), fixed synthetic weather

source("helper.R")

crop <- oryza_crop("IR72")
soil <- oryza_soil("potential")
control <- oryza_control()
control$modelstart <- as.Date("1992-01-01")
control$cropstart <- 0L
control$max_duration <- 160L
control$water_limited <- FALSE

wth <- .oryza_synth_weather(180L)
out <- oryza(crop, wth, soil, control)

expect_true(is.data.frame(out))
expect_true(nrow(out) > 50)
expect_true(all(c("date", "step", "DVS", "LAI", "WSO", "WRR14", "DAE", "CROPSTA") %in% names(out)))
expect_true(inherits(out$date, "Date"))
expect_equal(out$date[1], as.Date("1992-01-01"))

# Anchors from fixed seed=1 synthetic weather (see capture_cpp_anchors.R)
expect_equal(max(out$DVS), 2.018796, tolerance = 1e-5)
expect_equal(max(out$LAI), 6.956648, tolerance = 1e-4)
expect_equal(utils::tail(out$WSO, 1), 2921.962, tolerance = 0.1)
expect_true(max(out$CROPSTA) >= 4)
expect_true(all(is.finite(out$WSO)))
