\name{oryza_crop}

\alias{oryza_crop}

\title{
ORYZA crop parameters
}

\description{
Returns a list of crop parameters for the ORYZA model. Built-in varieties are
provided as INI files under \code{inst/oryza/crop}. See Bouman et al. (2001) for
parameter interpretation.
}

\usage{
oryza_crop(name = "")
}

\arguments{
\item{name}{character. Either the name of a crop that comes with the package
(see examples) or a path to a similarly formatted INI file. If empty, the
available crop names are returned}
}

\value{list}

\seealso{\code{\link{oryza}}, \code{\link{oryza_model}}}

\examples{
oryza_crop("")
crop <- oryza_crop("IR72")
str(crop)
}
