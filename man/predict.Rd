\name{predict}

\alias{predict}
\alias{predict,Rcpp_OryzaModel-method}

\title{Spatial ORYZA model predictions}

\description{
Make spatial predictions with an ORYZA model. First create a model object with crop and control settings, then call \code{predict} with a weather \code{SpatRasterDataset} and a soils \code{SpatRaster}.
}

\usage{
\S4method{predict}{Rcpp_OryzaModel}(object, weather, mstart, soils = NULL,
     soiltypes = NULL, filename = "", overwrite = FALSE, ...)
}

\arguments{
  \item{object}{ORYZA model (\code{Rcpp_OryzaModel})}
  \item{weather}{SpatRasterDataset with daily weather. Required sub-datasets:  \code{tmin}, \code{tmax}, \code{srad}. For water-limited runs also \code{prec}, \code{vapr}, and \code{wind}. All must share the same number of  layers; \code{tmin} must have Date time stamps}
  \item{mstart}{Date. Simulation start dates (one output layer per date)}
  \item{soils}{SpatRaster. Must include layer \code{elevation}. For water-limited runs also \code{soil} (1-based index into \code{soiltypes}) and \code{soildepth} (maximum rootable depth in m; used as \code{ZRTMS})}
  \item{soiltypes}{list of soil parameter lists (see \code{\link{oryza_soil}}).  Required for water-limited runs; optional for potential (defaults to \code{oryza_soil("potential")})}
  \item{filename}{character. Optional output filename}
  \item{overwrite}{logical. If \code{TRUE}, overwrite \code{filename}}
  \item{...}{additional arguments passed to \code{\link[terra]{writeRaster}}}
}

\value{
SpatRaster with one layer per element of \code{mstart}. Cell values are final grain biomass at 14\% humidty (\code{WRR14}, kg ha\eqn{^{-1}}, rounded).
}

\note{
Weather variables and units are the same as for the point simulation: 
\code{srad} in kJ m\eqn{^{-2}} d\eqn{^{-1}}, \code{vapr} in kPa, \code{prec} in mm d\eqn{^{-1}}, temperatures in \eqn{^\circ}C, wind in m s\eqn{^{-1}}.
}

\seealso{
\code{\link{oryza_model}}, \code{\link{oryza}}
}


\examples{
f <- system.file("extdata/Philippines_IRRI.csv", package="meteor")
w <- read.csv(f)[1:600, ]
w$date <- as.Date(w$date)
	
r <- terra::rast(nrow=2, ncol=2, nlyr=nrow(w))
v <- c("srad", "tmin", "tmax", "vapr", "wind", "prec")
wthfun <- function(i) {
	a <- rep(0, 4)		
	if ((i == "tmin") || (i == "tmax")) {
		a <- c(-2, 0, 2, 4)
	}
	out <- terra::setValues(r, rbind(w[,i]+a[1], w[,i]+a[2], w[,i]+a[3], w[,i]+a[4]))
	terra::time(out) <- as.Date(w$date)
	out 
}

wth <- terra::sds(lapply(v, wthfun))
names(wth) <- v

crop <- oryza_crop("IR72")
soil <- oryza_soil("paddy")
control <- oryza_control()
control$water_limited <- FALSE
control$SWITIR <- 0L
mod <- oryza_model(crop, w, soil, control)

soils <- terra::rast(nrow=2, ncol=2, names="elevation", vals=100)

p <- predict(mod, wth, mstart=as.Date(c("1979-01-05", "1979-05-05")), soils=soils)

}
