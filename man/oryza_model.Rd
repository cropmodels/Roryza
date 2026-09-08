\name{oryza_model}

\alias{oryza_model}
\alias{run}
\alias{crop<-}
\alias{soil<-}
\alias{control<-}
\alias{weather<-}
\alias{crop<-,Rcpp_OryzaModel,list-method}
\alias{soil<-,Rcpp_OryzaModel,list-method}
\alias{control<-,Rcpp_OryzaModel,list-method}
\alias{weather<-,Rcpp_OryzaModel,data.frame-method}
\alias{run,Rcpp_OryzaModel-method}

\title{ORYZA model object interface}

\description{
Create an ORYZA model object and run it. Through this interface you first create a model object and then you run it. You can change single parameters and run the model again.
}

\usage{
oryza_model(crop, weather, soil, control)
run(x, ...)
crop(x) <- value
soil(x) <- value
control(x) <- value
weather(x) <- value
}

\arguments{
\item{crop}{list. Crop parameters (see \code{\link{oryza_crop}})}
\item{weather}{data.frame with weather data. See Details}
\item{soil}{list. Soil parameters (see \code{\link{oryza_soil}})}
\item{control}{list. Model control options (see \code{\link{oryza_control}})}
\item{x}{ORYZA model object (\code{Rcpp_OryzaModel})}
\item{value}{crop, weather, soil, or control data, as above}
\item{\dots}{additional arguments. \code{stopError} (logical): if an error occurs and
\code{stopError} is \code{TRUE}, an error is thrown; otherwise a warning is given and
partial simulation data may be returned}
}

\value{
\code{oryza_model} returns an \code{OryzaModel} object.
\code{run} returns a data.frame of daily simulation results.
}

\details{
The weather data must be passed as a data.frame with the following variables and units.

\tabular{llll}{
	\tab \bold{variable} \tab \bold{description} \tab \bold{class / unit} \cr
	\tab \code{date} \tab "Date" class variable \tab - \cr
	\tab \code{srad} \tab Solar radiation  \tab kJ m-2 day-1 \cr
	\tab \code{tmin} \tab Minimum temperature \tab degrees C \cr
	\tab \code{tmax} \tab Maximum temperature \tab degrees C \cr
	\tab \code{vapr} \tab Vapor pressure \tab kPa \cr
	\tab \code{wind} \tab Wind speed \tab m s-1 \cr
	\tab \code{prec} \tab Precipitation \tab mm day-1 \cr
}

There should be no time gaps between the days in the data.frame.
}

\seealso{
\code{\link{oryza}}, \code{\link{oryza_crop}}, \code{\link{oryza_soil}},
\code{\link{oryza_control}}, \code{\link{predict}}
}


\references{
Bouman, B. A. M., Kropff, M. J., Tuong, T. P., Wopereis, M. C. S., ten Berge, H. F. M., & van Laar, H. H. (2001). ORYZA2000 : modeling lowland rice. International Rice Research Institute / Wageningen University. https://books.irri.org/9712201716_content.pdf
}


\examples{
# weather data (IRRI, Los Banos)
f <- system.file("extdata/Philippines_IRRI.csv", package="meteor")
wth <- read.csv(f)
wth$date <- as.Date(wth$date)

crop <- oryza_crop("IR72")
soil <- oryza_soil("paddy")
control <- oryza_control()
control$modelstart <- as.Date("1991-01-15")
control$latitude <- 14.18
control$elevation <- 21
control$water_limited <- TRUE
control$SWITIR <- 0L

m <- oryza_model(crop, wth, soil, control)
out <- run(m)
plot(out$date, out$LAI, type="l", xlab="Date", ylab="LAI")

# change soil / control and re-run
soil(m) <- oryza_soil("potential")
control$water_limited <- FALSE
control(m) <- control
out2 <- run(m)
plot(out2$date, out2$WSO, type="l", xlab="Date", ylab="Storage organs (kg/ha)")
}
