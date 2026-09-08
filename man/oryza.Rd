\name{oryza}
\alias{oryza}
\title{ORYZA rice crop growth model}

\description{
Run the ORYZA rice crop growth model (potential or water-limited).
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
data.frame with daily simulation results.
}

\seealso{
\code{\link{oryza_model}}, \code{\link{oryza_crop}}, \code{\link{oryza_soil}}, \code{\link{oryza_control}}
}
\details{
Weather variables

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
# weather data (IRRI, Los Banos)
f <- system.file("extdata/Philippines_IRRI.csv", package="meteor")
wth <- read.csv(f)
wth$date <- as.Date(wth$date)
head(wth)

crop <- oryza_crop("IR72")
soil <- oryza_soil("potential")
control <- oryza_control()
control$modelstart <- as.Date("1991-01-15")
control$latitude <- 14.18
control$elevation <- 21

# potential production
out <- oryza(crop, wth, soil, control)
tail(out)
plot(out$date, out$LAI, type="l", xlab="Date", ylab="LAI")

# water-limited (rainfed PADDY)
soil <- oryza_soil("paddy")
control$water_limited <- TRUE
control$SWITIR <- 0L
outw <- oryza(crop, wth, soil, control)
plot(outw$date, outw$WSO, type="l", xlab="Date", ylab="Storage organs (kg/ha)")
}
