# Roryza

R package for the **ORYZA2000** rice crop growth simulation model.

The API mirrors [RWofost](https://github.com/cropmodels/Rwofost): daily weather (`meteor` units), crop/soil/control parameter lists (INI files under `inst/oryza/`), and either a one-shot `oryza()` call or an `oryza_model()` object with `run()`.

Reference FORTRAN/FSE sources, rebuild scripts, and numerical tests live under `dev/` (not part of the installed package).

## Install

```r
# Windows: Rtools; macOS/Linux: a C++ toolchain
remotes::install_github("cropmodels/Roryza")
```

## Example

```r
library(Roryza)

crop <- oryza_crop("IR72")
control <- oryza_control()
soil <- oryza_soil("paddy")
control$water_limited <- TRUE

# weather: data.frame with date, srad, tmin, tmax, prec, wind, vapr
out <- oryza(crop, weather, soil, control)
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

## Tests

```r
# C++ API smoke + numeric anchors (always)
tinytest::test_all()

# Also runs FORTRAN/FSE scenarios when at_home and
# dev/inst/oryza3.exe exists (rebuild first):
#   Rscript dev/tools/build_oryza3.R
```

## FORTRAN reference (`dev/`)

```bash
Rscript dev/tools/build_oryza3.R
Rscript dev/tests/test_oryza3_rebuild.R
```

Binary is written to `dev/inst/oryza3.exe`.

## References

Bouman et al. (2001) *ORYZA2000: modeling lowland rice*. IRRI / Wageningen University.
