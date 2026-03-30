\name{MarkerPlot}
\alias{MarkerPlot}
\title{Plots a set of user defined markers using a seurat object.}
\usage{
# MarkerPlot(Obj, genes_df)
}
\description{
Plots an annotated and optionally clustered dotplot used a set of user-defined annotated genes and a \code{SeuratObject}.}
\arguments{
\item{obj}{\code{SeuratObject}.}
\item{genes}{\code{data.frame} with gene names in column 1 and an annotation in column 2.}
\item{margin_factor}{\code{numeric}, margin size modulator, increase or decrease to optimize annotation spacing, default=0.5}
\item{maxsize}{\code{numeric}, maximum size of dots, default=4}
\item{label.fontsize}{\code{numeric}, fontsize of annotations, default=3}
\item{assay}{Assay of the Seurat object to use for values, default="RNA"}
\item{show.annotations}{\code{boolean}, whether or not to label annotations, default=\code{TRUE}}
\item{cluster}{\code{boolean}, whether or not to cluster the x axis by correlation default=\code{TRUE}}
}
\value{
Returns a grob containing a dotplot of input genes.
}
\examples{
# MarkerPlot(Obj, genes_df)
# MarkerPlot(Obj, genes_df, cluster=F) # no clustering
# MarkerPlot(Obj, genes_df, show.annotations=F, margin_factor=0.2) # no annotation labels 
# MarkerPlot(Obj, genes_df, assay="SCT") # use the SCT assay to inform gene info
}


