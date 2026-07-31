# Roryza

R package for the **ORYZA2000** rice crop growth simulation model.

The computational engine is the original ORYZA2000 FORTRAN (FSE 2.1) shipped in the package and run as a subprocess — the same model used for reference validation. The R API mirrors [RWofost](https://github.com/cropmodels/Rwofost): daily weather (`meteor` units), crop/soil/control parameter lists, and either a one-shot `oryza()` call or an `oryza_model()` object with `run()` (C++ path).

## Features

| Mode | Status |
|------|--------|
| Potential production | Full (FORTRAN) |
| Water-limited PADDY | Full (FORTRAN) |
| Water-limited SAHEL / LOWBAL / SOILPF | Full (FORTRAN) |
| Water-limited SAWAH | Sources shipped; shipped Windows binary crashes in DRSAWA |
| Irrigation (SWITIR 0–6) | Full (FORTRAN) |
| Nitrogen balance (NCROP2 + NSOIL) | Full (FORTRAN) |
| Numerical parity with ORYZA2000 | Engine **is** ORYZA2000 |

## Install

```r
# Windows: Rtools recommended for the optional C++ path
remotes::install_github("cropmodels/Roryza")
```

A prebuilt `oryza3` binary is included under `inst/bin/`. To rebuild from `src/fse` sources:

```bash
# gfortran on PATH (Rtools on Windows)
Rscript tools/build_oryza3.R
```

## Example — FORTRAN engine (default)

```r
library(Roryza)

crop <- oryza_crop("IR72")
control <- oryza_control()
soil <- oryza_soil("paddy")
control$water_limited <- TRUE
control$WATBAL <- "PADDY"
control$SWITIR <- 6

# weather: date, srad (kJ m-2 d-1), tmin, tmax, prec, wind, vapr
out <- oryza(crop, weather, soil, control)           # engine = "fse"
```

## Example — nitrogen-limited

```r
control$nitrogen_limited <- TRUE
# optional: control$FERTIL <- c(dae1, kgN1, dae2, kgN2, ...)
out <- oryza(crop, weather, soil, control)
```

## Direct FSE rundir API

```r
wd <- oryza_fse_prepare(watbal = "SAWAH", water_limited = TRUE, prdel = 1)
res <- oryza_fse(wd)          # full RES.DAT table
```

## Weather units

| Column | Unit |
|--------|------|
| `date` | `Date` |
| `srad` | kJ m⁻² day⁻¹ |
| `tmin`, `tmax` | °C |
| `vapr` | kPa |
| `wind` | m s⁻¹ |
| `prec` | mm day⁻¹ |

## Validation

```r
# From package root after INSTALL:
Rscript tests/validate_fortran.R
```

This runs potential, all five water balances, and N-limited cases, and checks repeatability to absolute tolerance 0 (identical subprocess outputs).

## Object API (experimental C++ path)

```r
m <- oryza_model(crop, weather, soil, control)
out <- run(m)
# or: oryza(crop, weather, soil, control, engine = "cpp")
```

The C++ path is incomplete relative to ORYZA2000; use `engine = "fse"` for full physics and FORTRAN equality.

## References

Bouman et al. (2001) *ORYZA2000: modeling lowland rice*. IRRI / Wageningen University.
