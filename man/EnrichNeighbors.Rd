\name{EnrichNeighbors}
\alias{EnrichNeighbors}
\title{Check if any meta data variables are enriched in pre-defined neighborhoods.}
\usage{
# EnrichNeighbors(seuratObj, neighborCol, group.by, exclude.seeds=F, seedCol=NULL, split.by=NULL)
}
\description{
Identifies if any levels of a user-defined meta data variable are enriched in 
pre-defined neighborhoods.}
\arguments{
\item{seuratObj}{\code{SeuratObject}.}
\item{neighborCol}{\code{character}, column name of neighborhood designation in obj meta.data; 
must be boolean.}
\item{group.by}{\code{character}, column name of meta.data variable to query}
\item{exclude.seeds}{\code{boolean}, whether or not to remove seed cells when querying enrichment.
This does NOT necessarily include all seeds; only those neighboring another seed. Default=FALSE}
\item{seedCol}{\code{character}, if exclude.seeds == TRUE, column name of seed designation in 
obj meta.data; must be boolean. Default=NULL.}
\item{split.by}{\code{character}, optional column name of meta.data variable to split enrichment analysis.
E.g., for identifying enrichment differences between samples, conditions, etc. Default=NULL.}

}
\value{
Returns a \code{data.frame} containing enrichment stats for each group in neighborhoods vs. background.
}
\examples{
# Obj <- findNeighbors(Obj, seed.col="isTCell", nNeigh = 8, neighborCol = "TCell_Neighbor")
# enrichDF <- EnrichNeighbors(Obj, "TCell_Neighbor", "CellType")
# enrichDF_split <- EnrichNeighbors(Obj, neighborCol="TCell_Neighbor", group.by="CellType", split.by="Genotype")
}


