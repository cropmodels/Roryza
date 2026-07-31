# Water-limited PADDY: rainfed and irrigated (C++ engine)

source("helper.R")

crop <- oryza_crop("IR72")
soil <- oryza_soil("paddy")
control <- oryza_control()
control$modelstart <- as.Date("1992-01-01")
control$cropstart <- 0L
control$max_duration <- 150L
control$water_limited <- TRUE
control$SWITIR <- 0L

# Same RNG stream as potential anchors: generate 180 days, use first 160
wth <- .oryza_synth_weather(180L)[seq_len(160L), , drop = FALSE]
out <- oryza(crop, wth, soil, control)

expect_true(is.data.frame(out))
expect_true(nrow(out) > 50)
expect_true(all(c("TRW", "WL0", "IR", "MSKPA1", "PCEW", "LESTRS") %in% names(out)))

expect_equal(max(out$DVS), 2.018796, tolerance = 1e-5)
expect_equal(max(out$LAI), 3.718733, tolerance = 1e-4)
expect_equal(utils::tail(out$WSO, 1), 2235.800, tolerance = 0.1)
expect_true(max(out$TRW) > 0)
expect_true(all(out$PCEW >= 0 & out$PCEW <= 1.0000001))

# Rainfed yield should be below potential (~2922) on same weather
expect_true(utils::tail(out$WSO, 1) < 2922)

# Irrigated SWITIR=2
control$SWITIR <- 2L
control$WL0MIN <- 10
control$IRRI <- 50
out2 <- oryza(crop, wth, soil, control)
expect_true(nrow(out2) > 50)
expect_equal(utils::tail(out2$WSO, 1), 2670.130, tolerance = 0.1)
expect_true(sum(out2$IR > 0) > 0)
