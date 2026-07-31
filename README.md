# Roryza

R package for the **ORYZA2000** rice crop growth simulation model.

Either a one-shot `oryza()` call or an `oryza_model()` object with `run()`.

The API mirrors [RWofost](https://github.com/cropmodels/Rwofost)

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


## References

Bouman, B. A. M., Kropff, M. J., Tuong, T. P., Wopereis, M. C. S., ten Berge, H. F. M., & van Laar, H. H. (2001). ORYZA2000 : modeling lowland rice. ORYZA2000: modeling lowland rice. International Rice Research Institute/ Wageningen University. https://books.irri.org/9712201716_content.pdf
