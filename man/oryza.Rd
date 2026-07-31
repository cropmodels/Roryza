\name{oryza}
\alias{oryza}
\title{ORYZA2000 rice crop growth model}
\description{
Run the ORYZA2000 rice crop growth model (potential or water-limited).
}
\usage{
oryza(crop, weather, soil, control)
}
\arguments{
\item{crop}{list. Crop parameters (see \code{\link{oryza_crop}})}
\item{weather}{data.frame with weather data}
\item{soil}{list. Soil parameters (see \code{\link{oryza_soil}}); used when \code{water_limited=TRUE}}
\item{control}{list. Model control options (see \code{\link{oryza_control}})}
}
\value{
data.frame with daily simulation results. Water-limited runs also return
\code{TRW}, \code{WL0}, \code{IR}, \code{MSKPA1}, \code{PCEW}, \code{LESTRS}.
}
\seealso{
\code{\link{oryza_model}}, \code{\link{oryza_crop}}, \code{\link{oryza_soil}}, \code{\link{oryza_control}}
}
\details{
Weather columns (same units as RWofost / meteor):

\tabular{lll}{
\code{date} \tab Date \tab - \cr
\code{srad} \tab Solar radiation \tab kJ m-2 day-1 \cr
\code{tmin}, \code{tmax} \tab Temperature \tab degrees C \cr
\code{vapr} \tab Vapor pressure \tab kPa \cr
\code{wind} \tab Wind speed \tab m s-1 \cr
\code{prec} \tab Precipitation \tab mm day-1 \cr
}

Set \code{control$water_limited <- TRUE} and supply a PADDY soil for water-limited production.
Irrigation is controlled by \code{control$SWITIR} (0=rainfed, 2=pond depth, 3=soil KPa, 6=stage thresholds).
}
\examples{
\dontrun{
crop <- oryza_crop("IR72")
control <- oryza_control()
soil <- oryza_soil("potential")
# out <- oryza(crop, weather, soil, control)
}
}
