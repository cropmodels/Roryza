\name{predict}

\alias{predict}
\alias{predict,Rcpp_OryzaModel-method}

\title{Spatial ORYZA model predictions}

\description{
Make spatial predictions with an ORYZA2000 model. First create a model object with
crop and control settings, then call \code{predict} with a weather
\code{SpatRasterDataset} and a soils \code{SpatRaster}.
}

\usage{
\S4method{predict}{Rcpp_OryzaModel}(object, weather, mstart, soils = NULL,
     soiltypes = NULL, filename = "", overwrite = FALSE, ...)
}

\arguments{
  \item{object}{ORYZA model (\code{Rcpp_OryzaModel})}
  \item{weather}{SpatRasterDataset with daily weather. Required sub-datasets:
    \code{tmin}, \code{tmax}, \code{srad}. For water-limited runs also
    \code{prec}, \code{vapr}, and \code{wind}. All must share the same number of
    layers; \code{tmin} must have Date time stamps}
  \item{mstart}{Date. Simulation start dates (one output layer per date)}
  \item{soils}{SpatRaster. Must include layers \code{latitude} and \code{elevation}.
    For water-limited runs also \code{soil} (1-based index into \code{soiltypes})
    and \code{soildepth} (maximum rootable depth in m; used as \code{ZRTMS})}
  \item{soiltypes}{list of soil parameter lists (see \code{\link{oryza_soil}}).
    Required for water-limited runs; optional for potential (defaults to
    \code{oryza_soil("potential")})}
  \item{filename}{character. Optional output filename}
  \item{overwrite}{logical. If \code{TRUE}, overwrite \code{filename}}
  \item{...}{additional arguments passed to \code{\link[terra]{writeRaster}} /
    \code{writeStart} via \code{wopt}}
}

\value{
SpatRaster with one layer per element of \code{mstart}. Cell values are final
storage organ biomass (\code{WSO}, kg ha\eqn{^{-1}}, rounded).
}

\note{
Requires the \pkg{terra} package (\code{Suggests}). Weather units match the point
API: \code{srad} in kJ m\eqn{^{-2}} d\eqn{^{-1}}, \code{vapr} in kPa,
\code{prec} in mm d\eqn{^{-1}}, temperatures in \eqn{^\circ}C, wind in m s\eqn{^{-1}}.
}

\seealso{
\code{\link{oryza_model}}, \code{\link{oryza}}
}
