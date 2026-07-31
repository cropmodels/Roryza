\name{oryza_model}

\alias{oryza_model}
\alias{oryza_crop}
\alias{oryza_soil}
\alias{oryza_control}
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

\title{ORYZA2000 model object interface}

\description{
Create and run an ORYZA2000 model object, and load crop, soil and control parameters.
Through this interface you first create a model object and then you run it. You can
change single parameters and run the model again.
}

\usage{
oryza_model(crop, weather, soil, control)
oryza_crop(name = "")
oryza_soil(name = "")
oryza_control(filename = "")
run(x, ...)
crop(x) <- value
soil(x) <- value
control(x) <- value
weather(x) <- value
}

\arguments{
\item{crop}{list. Crop parameters (see \code{\link{oryza_crop}})}
\item{weather}{data.frame with weather data. See \code{\link{oryza}}}
\item{soil}{list. Soil parameters (see \code{\link{oryza_soil}})}
\item{control}{list. Model control options (see \code{\link{oryza_control}})}
\item{name}{name of built-in crop or soil (\code{"IR72"}, \code{"paddy"}, \code{"potential"}), or path to an INI file}
\item{filename}{path to control INI file (default: package \code{control.ini})}
\item{x}{ORYZA model object (\code{Rcpp_OryzaModel})}
\item{value}{crop, weather, soil, or control data, as above}
\item{\dots}{additional arguments. \code{stopError} (logical): if an error occurs and
\code{stopError} is \code{TRUE}, an error is thrown; otherwise a warning is given and
partial simulation data may be returned}
}

\value{
\code{oryza_model} returns an \code{OryzaModel} object.
\code{oryza_crop}, \code{oryza_soil}, and \code{oryza_control} return named lists.
\code{run} returns a data.frame of daily simulation results.
}

\seealso{
\code{\link{oryza}}
}

\examples{
\dontrun{
crop <- oryza_crop("IR72")
soil <- oryza_soil("paddy")
control <- oryza_control()
control$water_limited <- TRUE
m <- oryza_model(crop, weather, soil, control)
out <- run(m)

crop(m) <- oryza_crop("IR72")
out2 <- run(m)
}
}
