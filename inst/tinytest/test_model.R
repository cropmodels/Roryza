# Object API: oryza_model + run()

source("helper.R")

crop <- oryza_crop("IR72")
soil <- oryza_soil("potential")
control <- oryza_control()
control$modelstart <- as.Date("1992-01-01")
control$cropstart <- 0L
control$max_duration <- 160L

wth <- .oryza_synth_weather(180L)
m <- oryza_model(crop, wth, soil, control)
expect_true(inherits(m, "Rcpp_OryzaModel") || inherits(m, "C++Object") || !is.null(m))

out <- run(m)
expect_true(is.data.frame(out))
expect_true(nrow(out) > 50)
expect_equal(utils::tail(out$WSO, 1), 2923.095, tolerance = 0.1)

# one-shot oryza() should match object API on same inputs
out2 <- oryza(crop, wth, soil, control)
expect_equal(out$WSO, out2$WSO, tolerance = 1e-8)
expect_equal(out$DVS, out2$DVS, tolerance = 1e-8)
