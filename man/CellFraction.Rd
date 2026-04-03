\name{CellFraction}
\alias{CellFraction}
\title{Extract coordinates of all polygons from a Seurat image slot.}
\usage{
CellFraction(obj, group.by, split.by, label=T, label.size=3, min.perc=5, 
             tilt=T, colors=NULL, return.data=F)
}
\description{
Extract cell polygon coordinates and meta.data from a specific image within a \code{SeuratObject}.}
\arguments{
\item{obj}{\code{SeuratObject}.}
\item{group.by}{\code{character}, name of meta.data column to use for fill variable.}
\item{split.by}{\code{character}, name of meta.data column to use for x-axis variable.}
\item{label}{\code{boolean}, whether or not to draw \% labels on columns, default=TRUE.}
\item{label.size}{\code{integer}, fontsize for labels, if drawn, default=3.}
\item{min.perc}{\code{numeric}, minimum percentage for a label to be drawn, default=5.}
\item{tilt}{\code{boolean}, whether or not to angle the x-axis text, default=TRUE.}
\item{colors}{\code{character}, fill colours, default=NULL.}
\item{return.data}{\code{boolean}, if TRUE, returns data instead of a plot. Default=FALSE.}
}
\value{
Returns a grob or a dataframe with cell fraction information
}
\examples{
# draw a barplot with cell type fractions by condition:
# CellFraction(merfish, group.by="CellType", split.by="Condition")
# get a data.frame with cell type fractions by condition:
# CellFraction(merfish, group.by="CellType", split.by="Condition", return.data=T)
}


