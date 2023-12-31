% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.R
\name{map.habitat}
\alias{map.habitat}
\title{Map species distribution of habitat specialist.}
\usage{
map.habitat(
  longlat,
  layer,
  move = TRUE,
  mcp = FALSE,
  points = FALSE,
  eval = TRUE
)
}
\arguments{
\item{longlat}{Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.}

\item{layer}{RasterLayer object representing the presence/absence (1/0) of a single habitat type.}

\item{move}{If TRUE, identifies and moves presence records to closest cells with suitable habitat. Use when spatial error might put records outside the correct patch.}

\item{mcp}{If TRUE, all habitat patches inside the minimum convex hull polygon encompassing all occurrence records are converted to presence.}

\item{points}{If TRUE, force map to include cells with presence records even if suitable habitat was not identified.}

\item{eval}{If TRUE, build a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).}
}
\value{
One raster object and, if eval = TRUE, a matrix with EOO (from raw data), EOO (from model), AOO (from raw data) and AOO (from model).
}
\description{
Mapping of all habitat patches where the species is known to occur.
}
\details{
In many cases a species has a very restricted habitat and we generally know where it occurs. In such cases using the distribution of the known habitat patches may be enough to map the species.
}
