\name{getPolygons}
\alias{getPolygons}
\title{Extract coordinates of all polygons from a Seurat image slot.}
\usage{
getPolygons(obj, fov)
}
\description{
Extract cell polygon coordinates and meta.data from a specific image within a \code{SeuratObject}.}
\arguments{
\item{obj}{\code{SeuratObject}.}
\item{fov}{\code{character}, name of image to use.}
}
\value{
Returns a data frame containing polygon info and meta data for each cell in the defined image.
}
\examples{
# coordDF <- getPolygons(Seurat, "FOV_1")
}


