\name{oryza_model}
\alias{oryza_model}
\alias{oryza_crop}
\alias{oryza_soil}
\alias{oryza_control}
\alias{run}
\title{ORYZA2000 model object interface}
\description{
Create and run an ORYZA2000 model object, and load crop, soil and control parameters.
}
\usage{
oryza_model(crop, weather, soil, control)
oryza_crop(name = "")
oryza_soil(name = "")
oryza_control(filename = "")
}
\arguments{
\item{crop}{crop parameter list}
\item{weather}{weather data.frame}
\item{soil}{soil parameter list}
\item{control}{control parameter list}
\item{name}{name of built-in crop or soil (\code{"IR72"}, \code{"paddy"}, \code{"potential"}), or path to an INI file}
\item{filename}{path to control INI file (default: package control.ini)}
}
\value{
\code{oryza_model} returns an \code{OryzaModel} object. Parameter helpers return named lists.
Use \code{run(m)} to simulate. Replacement methods \code{crop<-}, \code{soil<-},
\code{control<-}, and \code{weather<-} update the model object.
}
\examples{
\dontrun{
crop <- oryza_crop("IR72")
soil <- oryza_soil("paddy")
control <- oryza_control()
control$water_limited <- TRUE
m <- oryza_model(crop, weather, soil, control)
out <- run(m)
}
}
