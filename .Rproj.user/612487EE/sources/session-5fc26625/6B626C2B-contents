\name{PlotPolygons}
\alias{PlotPolygons}
\title{Extract coordinates of all polygons from a Seurat image slot.}
\usage{
PlotPolygons(coordDF, fillVar, border_color="black", fillColor=NULL, fill_lims=NULL, 
             legend="right", background=NULL, legend.name=NULL)
}
\description{
Uses cell polygon coordinates and meta.data from a specific image within a \code{SeuratObject} to plot a polygon image.}
\arguments{
\item{coordDF}{\code{data.frame}, output from getPolygons().}
\item{fillVar}{\code{character}, name of meta.data variable to for colouration. Ignored if background is non-null.}
\item{border_color}{\code{character}, colour for polygon outlines.}
\item{fillColor}{\code{vector} of colours to use for the fill gradient, defaults to reds.}
\item{fill_lims}{\code{vector}, numeric vector of length 2, manual limits for fill colour.}
\item{legend}{\code{character}, legend position, inherited from ggplot2::theme(legend.position).}
\item{background}{\code{character}, if set to a colour instead of NULL, all plotted cells will be this colour. Intended use is for creating a background layer.}
\item{legend.name}{\code{character}, legend title, defaults to fillVar.}
}
\value{
Returns a grob containing spatially plotted polygons.
}
\examples{
# coordDF <- getPolygons(Seurat, "FOV_1")

# PlotPolygons(coordDF, "Trp53_Expr") # colour polygons by meta.data column
# PlotPolygons(coordDF, background="grey") # makes all polygons grey
}


