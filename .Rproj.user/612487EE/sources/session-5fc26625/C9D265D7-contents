\name{prepStack}
\alias{prepStack}
\title{Setup one layer of a layered polygon plot for stacking.}
\usage{
prepStack(p, coordDF, first=F)
}
\description{Sets up all parameters for grob stacking.}
\arguments{
\item{p}{\code{grob}, output from PlotPolygons().}
\item{coordDF}{\code{data.frame}, output from getPolygons(), must include all cells within the FOV of interest.}
\item{first}{\code{boolean}, whether or not the image is the base layer of the image stack.}
}
\value{
Returns several grobs containing spatially plotted polygons from the same FOV overlaid on top of each other.
}
\examples{
# coordDF <- getPolygons(Seurat, "FOV_1")

# p1 <- PlotPolygons(subset(coordDF, CellType == "T Cells"),
#                    "Trp53_Expr") # colour T cells by Trp53 expression

# p2 <- PlotPolygons(subset(coordDF, CellType != "T Cells"), 
#                    background="grey") # colour other cells grey

# prepStack(p2, coordDF, first=T)+ # use non-T Cells as base layer
#   prepStack(p1, coordDF) # overlay with coloured T Cells

}


