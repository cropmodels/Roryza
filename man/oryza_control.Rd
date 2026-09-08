\name{oryza_control}

\alias{oryza_control}

\title{
ORYZA control parameters
}

\description{
Returns a list of control and experiment parameters to run the ORYZA model. Important fields include \code{modelstart}, \code{cropstart}, \code{max_duration}, \code{latitude}, \code{CO2}, \code{water_limited}, and \code{WATBAL}. See Bouman et al. (2001) .
}

\usage{
oryza_control(filename = "")
}

\arguments{
\item{filename}{character. Path to a control INI file. If empty, default values are returned}
}

\note{
In the FORTRAN versions of ORYZA these parameters in the "experiment" file, not in the "control" file.
}

\value{list}

\seealso{\code{\link{oryza}}, \code{\link{oryza_model}}}

\examples{
contr <- oryza_control()
str(contr)
}
