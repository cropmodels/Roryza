\name{oryza_soil}

\alias{oryza_soil}

\title{
ORYZA soil parameters
}

\description{
Returns a list of soil parameters for the ORYZA model. Built-in soils are
provided as INI files under \code{inst/oryza/soil}. Use \code{"potential"} when
running potential production (water balance not used). See Bouman et al. (2001)
for parameter interpretation.
}

\usage{
oryza_soil(name = "")
}

\arguments{
\item{name}{character. Either the name of a soil that comes with the package
(see examples) or a path to a similarly formatted INI file. If empty, the
available soil names are returned}
}

\value{list}

\seealso{\code{\link{oryza}}, \code{\link{oryza_model}}}

\examples{
oryza_soil("")
soil <- oryza_soil("paddy")
str(soil)
}
