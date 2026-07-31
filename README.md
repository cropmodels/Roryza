# Roryza

R package for the **ORYZA2000** rice crop growth simulation model.

Either a one-shot `oryza()` call or an `oryza_model()` object with `run()`.

## Install

```r
install.packages("Roryza", repos = "https://rspatial.r-universe.dev")
```

## Example

```r
library(Roryza)

# weather data (IRRI, Los Banos)
f <- system.file("extdata/Philippines_IRRI.csv", package = "meteor")
wth <- read.csv(f)
wth$date <- as.Date(wth$date)

crop <- oryza_crop("IR72")
soil <- oryza_soil("potential")
control <- oryza_control()
control$modelstart <- as.Date("1991-01-15")
control$latitude <- 14.18
control$elevation <- 21

# potential production
out <- oryza(crop, wth, soil, control)
plot(out$date, out$LAI, type = "l")

# water-limited (rainfed PADDY)
soil <- oryza_soil("paddy")
control$water_limited <- TRUE
control$SWITIR <- 0L
outw <- oryza(crop, wth, soil, control)
plot(outw$date, outw$WSO, type = "l")
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

## References

Bouman, B. A. M., Kropff, M. J., Tuong, T. P., Wopereis, M. C. S., ten Berge, H. F. M., & van Laar, H. H. (2001). ORYZA2000 : modeling lowland rice. International Rice Research Institute / Wageningen University. https://books.irri.org/9712201716_content.pdf
