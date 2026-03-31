\name{findNeighbors}
\alias{findNeighbors}
\title{Identify closest neighbors of seed cells.}
\usage{
findNeighbors(seuratObj, seed.col=NULL, nNeigh=8, split.by=NULL, 
              xname=NULL, yname=NULL, neighborCol="isNeighbor")
}
\description{
Extract the k-closest neighbors of a set of pre-defined cells in a \code{SeuratObject}.}
\arguments{
\item{seuratObj}{\code{SeuratObject}.}
\item{seed.col}{\code{character}, name of meta.data column containing seed designations. Must be boolean.}
\item{nNeigh}{\code{integer}, number of k-nearest neighbors to select.}
\item{split.by}{\code{character}, name of meta.data column containing groups to split by \(e.g., 
if there are multiple areas or FOVs in the dataset\). Default=NULL.}
\item{xname}{\code{character}, name of meta.data column for x-coordinate if not simply 'x'. Default=NULL.}
\item{yname}{\code{character}, name of meta.data column for x-coordinate if not simply 'y'. Default=NULL.}
\item{neighborCol}{\code{character}, name for output column. Default="isNeighbor"}
}
\value{
Returns a Seurat object with an additional meta.data column indicating neighborhood status.
}
\examples{
# Obj <- findNeighbors(Obj, seed.col="isTCell", nNeigh = 8, neighborCol = "TCell_Neighbor")
}


